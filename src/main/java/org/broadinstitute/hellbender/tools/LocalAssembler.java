package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.util.SequenceUtil;
import org.apache.hadoop.fs.DF;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.MultiplePassReadWalker;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.SetSizeUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.LongFunction;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "experiment",
        oneLineSummary = "experiment",
        usageExample = "gatk LocalAssembler",
        programGroup = CoverageAnalysisProgramGroup.class
)
@BetaFeature
public class LocalAssembler extends MultiplePassReadWalker {
    public static final byte QMIN = 22;
    public static final int MIN_THIN_OBS = 4;

    @Override public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.PRIMARY_LINE);
    }

    @Override public void traverseReads() {
        final KmerSet<KmerAdjacency> kmerAdjacencySet = new KmerSet<>(1000000);
        final int nReads = kmerizeReadsPass(kmerAdjacencySet);
        final List<ContigImpl> contigs = buildContigs(kmerAdjacencySet);
        connectContigs(contigs);
        while ( removeThinContigs(contigs, kmerAdjacencySet) ) {}
        weldPipes(contigs);
        final List<Path> readPaths = new ArrayList<>(nReads);
        final Map<GapFill, List<List<PathPart>>> gapFillCountMap = new HashMap<>();
        pathReadsPass(kmerAdjacencySet, readPaths, gapFillCountMap);
        fillGaps(contigs, gapFillCountMap, kmerAdjacencySet);
        weldPipesAndPatchPaths(contigs, readPaths);
        final int nComponents = markComponents(contigs);
        final List<List<Contig>> cycles = new ArrayList<>();
        markCycles(contigs, cycles);

        contigs.sort((tig1, tig2) -> Integer.compare(tig1.getId(), tig2.getId()));
        writeDOT(contigs, "assembly.dot");
        writeContigs(contigs);
        writePaths(readPaths);
        writeCycles(cycles);
        System.out.println("There are " + nComponents + " assembly graph components.");
    }

    private int kmerizeReadsPass( final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        final int[] nReads = new int[1];
        forEachRead( (read, ref, feature, nReadsProcessed) -> {
            final byte[] calls = read.getBasesNoCopy();
            final byte[] quals = read.getBaseQualitiesNoCopy();
            KmerAdjacencyImpl.kmerize(calls, quals, QMIN, kmerAdjacencySet);
            nReads[0] += 1;
        });
        return nReads[0];
    }

    private static List<ContigImpl> buildContigs( final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        final List<ContigImpl> contigs = new ArrayList<>();
        for ( final KmerAdjacency kmerAdjacency : kmerAdjacencySet ) {
            if ( kmerAdjacency.getContig() == null ) {
                ContigImpl contig = null;
                final KmerAdjacency predecessor = kmerAdjacency.getSolePredecessor();
                if ( predecessor == null || predecessor.getSuccessorCount() > 1 ) {
                    contig = new ContigImpl(kmerAdjacency);
                } else {
                    final KmerAdjacency successor = kmerAdjacency.getSoleSuccessor();
                    if ( successor == null || successor.getPredecessorCount() > 1 ) {
                        contig = new ContigImpl(kmerAdjacency.rc());
                    }
                }
                if ( contig != null ) {
                    contigs.add(contig);
                }
            }
        }
        return contigs;
    }

    private static void connectContigs( final List<ContigImpl> contigs ) {
        final int nContigs = contigs.size();
        final KmerSet<ContigEndKmer> contigEnds = new KmerSet<>(2*nContigs);
        for ( int contigId = 0; contigId != nContigs; ++contigId ) {
            final ContigImpl contig = contigs.get(contigId);
            final KmerAdjacency fwdKmer = contig.getFirstKmer();
            final KmerAdjacency revKmer = contig.getLastKmer().rc();
            if ( fwdKmer == revKmer ) {
                contigEnds.findOrAdd(fwdKmer.getKVal(), kVal -> new ContigEndKmer(kVal, contig, ContigOrientation.BOTH));
            } else {
                contigEnds.findOrAdd(fwdKmer.getKVal(), kVal -> new ContigEndKmer(kVal, contig, ContigOrientation.FWD));
                contigEnds.findOrAdd(revKmer.getKVal(), kVal -> new ContigEndKmer(kVal, contig, ContigOrientation.REV));
            }
        }

        for ( int contigId = 0; contigId != nContigs; ++contigId ) {
            final Contig contig = contigs.get(contigId);

            final KmerAdjacency start = contig.getFirstKmer();
            final int predecessorCount = start.getPredecessorCount();
            if ( predecessorCount > 0 ) {
                final List<Contig> predecessors = contig.getPredecessors();
                final int mask = start.getPredecessorMask();
                for ( int call = 0; call != 4; ++call ) {
                    if ( (mask & (1 << call)) != 0 ) {
                        final ContigEndKmer contigEndKmer =
                                contigEnds.find(KmerAdjacency.reverseComplement(start.getPredecessorVal(call)));
                        switch ( contigEndKmer.getContigOrientation() ) {
                            case FWD:
                                predecessors.add(contigEndKmer.getContig().rc());
                                break;
                            case REV:
                                predecessors.add(contigEndKmer.getContig());
                                break;
                            case BOTH:
                                predecessors.add(contigEndKmer.getContig());
                                predecessors.add(contigEndKmer.getContig().rc());
                                break;
                        }
                    }
                }
            }

            final KmerAdjacency end = contig.getLastKmer();
            final int successorCount = end.getSuccessorCount();
            if ( successorCount > 0 ) {
                final List<Contig> successors = contig.getSuccessors();
                final int mask = end.getSuccessorMask();
                for ( int call = 0; call != 4; ++call ) {
                    if ( (mask & (1 << call)) != 0 ) {
                        final ContigEndKmer contigEndKmer = contigEnds.find(end.getSuccessorVal(call));
                        switch ( contigEndKmer.getContigOrientation() ) {
                            case FWD:
                                successors.add(contigEndKmer.getContig());
                                break;
                            case REV:
                                successors.add(contigEndKmer.getContig().rc());
                                break;
                            case BOTH:
                                successors.add(contigEndKmer.getContig());
                                successors.add(contigEndKmer.getContig().rc());
                                break;
                        }
                    }
                }
            }
        }
    }

    private static boolean removeThinContigs( final List<ContigImpl> contigs,
                                           final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        for ( final Contig contig : contigs ) {
            contig.setCut(false);
            contig.setAuxData(null);
            contig.rc().setAuxData(null);
        }

        for ( final Contig contig : contigs ) {
            if ( contig.getAuxData() != null ) continue;
            contig.setAuxData(new CutData());
            int children = 0;
            for ( final Contig nextContig : contig.getSuccessors() ) {
                if ( nextContig.getAuxData() == null ) {
                    findCuts(nextContig, contig);
                    children += 1;
                }
            }
            for ( final Contig nextContig : contig.getPredecessors() ) {
                if ( nextContig.getAuxData() == null ) {
                    findCuts(nextContig, contig);
                    children += 1;
                }
            }
            if ( children >= 2 ) {
                contig.setCut(true);
            }
        }

        return contigs.removeIf( tig -> {
            if ( tig.getMaxObservations() < MIN_THIN_OBS && !tig.isCut() ) {
                unlinkContig(tig, kmerAdjacencySet);
                return true;
            }
            return false;
        } );
    }

    private static CutData findCuts( final Contig contig, final Contig parent ) {
        final CutData cutData = new CutData();
        contig.setAuxData(cutData);
        for ( final Contig nextContig : contig.getSuccessors() ) {
            if ( nextContig == parent ) continue;
            CutData nextCutData = (CutData)nextContig.getAuxData();
            if ( nextCutData != null ) {
                cutData.minVisitNum = Math.min(cutData.minVisitNum, nextCutData.visitNum);
            } else {
                nextCutData = findCuts(nextContig, contig);
                cutData.minVisitNum = Math.min(cutData.minVisitNum, nextCutData.minVisitNum);
                if ( nextCutData.minVisitNum >= cutData.visitNum ) {
                    contig.setCut(true);
                }
            }
        }
        for ( final Contig nextContig : contig.getPredecessors() ) {
            if ( nextContig == parent ) continue;
            CutData nextCutData = (CutData)nextContig.getAuxData();
            if ( nextCutData != null ) {
                cutData.minVisitNum = Math.min(cutData.minVisitNum, nextCutData.visitNum);
            } else {
                nextCutData = findCuts(nextContig, contig);
                cutData.minVisitNum = Math.min(cutData.minVisitNum, nextCutData.minVisitNum);
                if ( nextCutData.minVisitNum >= cutData.visitNum ) {
                    contig.setCut(true);
                }
            }
        }
        return cutData;
    }

    private static void unlinkContig( final Contig contig, final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        final KmerAdjacency firstKmer = contig.getFirstKmer();
        firstKmer.removeAllPredecessors();
        final int firstKmerFinalCall = firstKmer.getFinalCall();
        for ( final Contig predecessor : contig.getPredecessors() ) {
            if ( predecessor != contig && predecessor != contig.rc() ) {
                predecessor.getLastKmer().removeSuccessor(firstKmerFinalCall, kmerAdjacencySet);
                if ( !predecessor.getSuccessors().remove(contig) ) {
                    throw new GATKException("failed to find predecessor link");
                }
            }
        }

        final KmerAdjacency lastKmer = contig.getLastKmer();
        lastKmer.removeAllSuccessors();
        final int lastKmerInitialCall = lastKmer.getInitialCall();
        for ( final Contig successor : contig.getSuccessors() ) {
            if ( successor != contig && successor != contig.rc() ) {
                successor.getFirstKmer().removePredecessor(lastKmerInitialCall, kmerAdjacencySet);
                if ( !successor.getPredecessors().remove(contig) ) {
                    throw new GATKException("failed to find successor link");
                }
            }
        }

        updateKmerContig( contig, null, 0 );
    }

    private static int updateKmerContig( final Contig oldContig, final Contig newContig, final int initialOffset ) {
        int offset = initialOffset;
        KmerAdjacency kmer = oldContig.getFirstKmer();
        do {
            kmer.setContig(newContig, offset);
            if ( newContig != null ) offset += 1;
        } while ( (kmer = kmer.getSoleSuccessor()) != null && kmer.getContig() == oldContig );
        return offset;
    }

    private static void weldPipes( final List<ContigImpl> contigs ) {
        for ( int contigId = 0; contigId != contigs.size(); ++contigId ) {
            final ContigImpl contig = contigs.get(contigId);
            if ( contig.getSuccessors().size() == 1 ) {
                final Contig successor = contig.getSuccessors().get(0);
                if ( successor != contig && successor != contig.rc() && successor.getPredecessors().size() == 1 ) {
                    final ContigImpl joinedContig = join(contig, successor);
                    final List<Contig> successors = joinedContig.getSuccessors();
                    if ( successors.remove(successor) ) {
                        successors.add(joinedContig);
                    }
                    if ( successors.remove(successor.rc()) ) {
                        successors.add(joinedContig.rc());
                    }
                    contigs.set(contigId, joinedContig);
                    if ( !contigs.remove(successor.canonical()) ) {
                        throw new GATKException("successor linkage is messed up");
                    }
                    contigId -= 1; // reconsider the new contig -- there might be more joining possible
                    continue;
                }
            }
            if ( contig.getPredecessors().size() == 1 ) {
                final Contig predecessor = contig.getPredecessors().get(0);
                if ( predecessor != contig && predecessor != contig.rc() && predecessor.getSuccessors().size() == 1 ) {
                    final ContigImpl joinedContig = join(predecessor, contig);
                    final List<Contig> predecessors = joinedContig.getPredecessors();
                    if ( predecessors.remove(predecessor) ) {
                        predecessors.add(joinedContig);
                    }
                    if ( predecessors.remove(predecessor.rc()) ) {
                        predecessors.add(joinedContig.rc());
                    }
                    contigs.set(contigId, joinedContig);
                    if ( !contigs.remove(predecessor.canonical()) ) {
                        throw new GATKException("predecessor linkage is messed up");
                    }
                    contigId -= 1; // reconsider
                }
            }
        }
    }

    private static void weldPipesAndPatchPaths( final List<ContigImpl> contigs, final List<Path> readPaths ) {
        for ( int contigId = 0; contigId != contigs.size(); ++contigId ) {
            final ContigImpl contig = contigs.get(contigId);
            if ( contig.getSuccessors().size() == 1 ) {
                final Contig successor = contig.getSuccessors().get(0);
                if ( successor != contig && successor != contig.rc() && successor.getPredecessors().size() == 1 ) {
                    final ContigImpl newContig = join(contig, successor);
                    contigs.set(contigId, newContig);
                    if ( !contigs.remove(successor.canonical()) ) {
                        throw new GATKException("linkage is messed up");
                    }
                    contigId -= 1; // reconsider the new contig -- there might be more joining possible
                    patchPaths(contig, successor, newContig, readPaths);
                }
            }
        }
    }

    private static ContigImpl join( final Contig predecessor, final Contig successor ) {
        final ContigImpl joinedContig = new ContigImpl(predecessor, successor);
        for ( final Contig contig : joinedContig.getPredecessors() ) {
            if ( contig != joinedContig ) {
                final List<Contig> successorContigs = contig.getSuccessors();
                successorContigs.set(successorContigs.indexOf(predecessor), joinedContig);
            }
        }
        for ( final Contig contig : joinedContig.getSuccessors() ) {
            if ( contig != joinedContig ) {
                final List<Contig> predecessorContigs = contig.getPredecessors();
                predecessorContigs.set(predecessorContigs.indexOf(successor), joinedContig);
            }
        }
        updateKmerContig(successor, joinedContig, updateKmerContig(predecessor, joinedContig, 0));
        return joinedContig;
    }

    private static void patchPaths( final Contig predecessor, final Contig successor, final Contig joinedContig,
                                    final List<Path> readPaths ) {
        final int predecessorMaxStop = predecessor.size() - Kmer.KSIZE + 1;
        final int successorMaxStop = successor.size() - Kmer.KSIZE + 1;
        for ( final Path path : readPaths ) {
            final List<PathPart> parts = path.getParts();
            int nParts = parts.size();
            for ( int partId = 0; partId != nParts; ++partId ) {
                final PathPart part = parts.get(partId);
                final Contig contig = part.getContig();
                if ( contig == predecessor ) {
                    final PathPart replacementPart = new PathPart(joinedContig, part.getStart(), part.getStop());
                    final int nextId = partId + 1;
                    if ( part.getStop() == predecessorMaxStop && nextId < nParts ) {
                        final PathPart nextPart = parts.get(nextId);
                        if ( nextPart.getContig() == successor && nextPart.getStart() == 0 ) {
                            replacementPart.setStop(predecessorMaxStop + nextPart.getStop());
                            parts.remove(nextId);
                            nParts -= 1;
                        }
                    }
                    parts.set(partId, replacementPart);
                } else if ( contig == successor.rc() ) {
                    final PathPart replacementPart = new PathPart(joinedContig.rc(), part.getStart(), part.getStop());
                    final int nextId = partId + 1;
                    if ( part.getStop() == successorMaxStop && nextId < nParts ) {
                        final PathPart nextPart = parts.get(nextId);
                        if ( nextPart.getContig() == predecessor.rc() && nextPart.getStart() == 0 ) {
                            replacementPart.setStop(successorMaxStop + nextPart.getStop());
                            parts.remove(nextId);
                            nParts -= 1;
                        }
                    }
                    parts.set(partId, replacementPart);
                } else if ( contig == successor ) {
                    parts.set(partId, new PathPart(joinedContig,
                                                    part.getStart() + predecessorMaxStop,
                                                    part.getStop() + predecessorMaxStop));
                } else if ( contig == predecessor.rc() ) {
                    parts.set(partId, new PathPart(joinedContig.rc(),
                                                    part.getStart() + successorMaxStop,
                                                    part.getStop() + successorMaxStop));
                }
            }
        }
    }
/*
    private void extendSinks( final List<ContigImpl> contigs,
                              final Map<Contig, String> contigNames,
                              final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        final Map<Contig, List<int[]>> extensions = new HashMap<>(contigs.size() * 3);
        for ( final Contig contig : contigs ) {
            if ( contig.getSuccessors().size() == 0 ) {
                extensions.put(contig, new ArrayList<>());
            }
            if ( contig.rc().getSuccessors().size() == 0 ) {
                extensions.put(contig.rc(), new ArrayList<>());
            }
        }
        forEachRead( (read, ref, feature, nReadsProcessed) -> {
            final byte[] calls = read.getBasesNoCopy();
            buildExtensions(calls, kmerAdjacencySet, extensions);
            SequenceUtil.reverseComplement(calls);
            buildExtensions(calls, kmerAdjacencySet, extensions);
        });
        for ( final Map.Entry<Contig, List<int[]>> entry : extensions.entrySet() ) {
            final List<int[]> callCounts = entry.getValue();
            if ( callCounts.size() > 0 ) {
                final Contig contig = entry.getKey();
                long kVal = contig.getLastKmer().getKVal();
                final StringBuilder sb = new StringBuilder(callCounts.size());
                for ( final int[] counts : callCounts ) {
                    int max = -1;
                    int argMax = -1;
                    for ( int idx = 0; idx < 4; ++idx ) {
                        final int count = counts[idx];
                        if ( count > max ) {
                            max = count;
                            argMax = idx;
                        } else if ( count == max ) {
                            argMax = -1;
                        }
                    }
                    if ( argMax == -1 ) {
                        break;
                    }
                    kVal = ((kVal << 2) | argMax) & KMASK;
                    final KmerAdjacency kmer = KmerAdjacency.find(kVal, kmerAdjacencySet);
                    if ( kmer != null && kmer.getContig() != null ) {
                        System.out.println(contigNames.get(contig) + " + " + sb + " -> " +
                                contigNames.get(kmer.getContig()) + ":" + kmer.getContigOffset());
                        sb.setLength(0);
                        break;
                    }
                    sb.append("ACGT".charAt(argMax));
                }
                if ( sb.length() > 0 ) {
                    System.out.println(contigNames.get(contig) + " + " + sb);
                }
            }
        }
    }

    private static void buildExtensions( final byte[] calls,
                                         final KmerSet<KmerAdjacency> kmerAdjacencySet,
                                         final Map<Contig, List<int[]>> extensions ) {
        long kVal = 0;
        int readOffset = 0;
        for ( final byte call : calls ) {
            kVal <<= 2;
            switch ( call ) {
                case 'C': case 'c': kVal += 1; break;
                case 'G': case 'g': kVal += 2; break;
                case 'T': case 't': kVal += 3; break;
            }
            if ( ++readOffset >= KSIZE ) {
                final KmerAdjacency kmer = KmerAdjacency.find(kVal & KMASK, kmerAdjacencySet);
                if ( kmer != null ) {
                    final Contig contig = kmer.getContig();
                    if ( contig != null ) {
                        int extensionLength = readOffset - KSIZE - kmer.getContigOffset();
                        if ( extensionLength > 0 ) {
                            // if contig.rc() is not a sink, the lookup will return null
                            final List<int[]> extension = extensions.get(contig.rc());
                            if ( extension != null ) {
                                int extensionOffset = 0;
                                while ( extensionLength > 0 ) {
                                    final int rcCall;
                                    switch ( calls[--extensionLength] ) {
                                        case 'A': case 'a': rcCall = 3; break;
                                        case 'C': case 'c': rcCall = 2; break;
                                        case 'G': case 'g': rcCall = 1; break;
                                        case 'T': case 't': rcCall = 0; break;
                                        default: rcCall = -1; break;
                                    }
                                    if ( rcCall >= 0 ) {
                                        while ( extensionOffset >= extension.size() ) {
                                            extension.add(new int[4]);
                                        }
                                        extension.get(extensionOffset)[rcCall] += 1;
                                    }
                                    extensionOffset += 1;
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }
    }
*/

    private static int markComponents( final List<ContigImpl> contigs ) {
        for ( final ContigImpl contig : contigs ) {
            contig.setComponentId(0);
        }

        int componentId = 0;
        for ( final ContigImpl contig : contigs ) {
            if ( contig.getComponentId() == 0 ) {
                contig.setComponentId(++componentId);
                markSuccessorComponents(contig);
                markSuccessorComponents(contig.rc());
            }
        }
        return componentId;
    }

    private static void markSuccessorComponents( final Contig contig ) {
        final int componentId = contig.getComponentId();
        for ( final Contig successor : contig.getSuccessors() ) {
            if ( successor.getComponentId() == 0 ) {
                successor.canonical().setComponentId(componentId);
                markSuccessorComponents(successor);
                markSuccessorComponents(successor.rc());
            }
        }
    }

    private static void markCycles( final List<ContigImpl> contigs, final List<List<Contig>> cycles ) {
        for ( final ContigImpl contig : contigs ) {
            contig.setCyclic(false);
            contig.setAuxData(DFSearchStatus.UNVISITED);
            contig.rc().setCyclic(false);
            contig.rc().setAuxData(DFSearchStatus.UNVISITED);
        }

        final List<Contig> visiting = new ArrayList<>(contigs.size());
        for ( final ContigImpl contig : contigs ) {
            if ( contig.getAuxData() == DFSearchStatus.UNVISITED ) {
                markSuccessorCycles(contig, visiting, cycles);
            }
            if ( contig.rc().getAuxData() == DFSearchStatus.UNVISITED ) {
                markSuccessorCycles(contig.rc(), visiting, cycles);
            }
        }
    }

    private static void markSuccessorCycles( final Contig contig,
                                             final List<Contig> visiting,
                                             final List<List<Contig>> cycles ) {
        visiting.add(contig);
        contig.setAuxData(DFSearchStatus.VISITING);
        for ( final Contig successor : contig.getSuccessors() ) {
            final Object successorState = successor.getAuxData();
            if ( successorState == DFSearchStatus.VISITING ) {
                markCycle(visiting, successor, cycles);
            } else if ( successorState == DFSearchStatus.UNVISITED ) {
                markSuccessorCycles(successor, visiting, cycles);
            }
        }
        contig.setAuxData(DFSearchStatus.VISITED);
        visiting.remove(visiting.size() - 1);
    }

    private static void markCycle( final List<Contig> visiting,
                                   final Contig contig,
                                   final List<List<Contig>> cycles ) {
        for ( int idx = visiting.size() - 1; idx >= 0; --idx ) {
            final Contig cyclicContig = visiting.get(idx);
            cyclicContig.setCyclic(true);
            if ( cyclicContig == contig ) {
                cycles.add(new ArrayList<>(visiting.subList(idx, visiting.size())));
                return;
            }
        }
        throw new GATKException("shouldn't be able to get here -- cycle-starting contig not found");
    }

/*
    private static void phaseBubbles( final List<ContigImpl> contigs ) {
        for ( final Contig contig : contigs ) {
            final List<Contig> predecessors = contig.getPredecessors();
            if ( predecessors.size() > 1 ) {
                final List<Contig> successors = contig.getSuccessors();
                if ( successors.size() > 1 ) {
                    final List<List<Long>> predecessorsObservations = new ArrayList<>(predecessors.size());
                    for ( final Contig predecessor : predecessors ) {
                        predecessorsObservations.add(predecessor.getLastKmer().getObservations());
                    }
                    final List<List<Long>> successorsObservations = new ArrayList<>(successors.size());
                    for ( final Contig successor : successors ) {
                        successorsObservations.add(successor.getFirstKmer().getObservations());
                    }
                    System.out.print(contig);
                    for ( final List<Long> list1 : predecessorsObservations ) {
                        for ( final List<Long> list2 : successorsObservations ) {
                            System.out.print('\t');
                            System.out.print(commonCount(list1, list2));
                        }
                        System.out.println();
                    }
                }
            }
        }
    }

    private static int commonCount( final List<Long> list1, final List<Long> list2 ) {
        final Iterator<Long> itr2 = list2.iterator();
        if ( !itr2.hasNext() ) return 0;
        long ele2 = itr2.next();
        int result = 0;
        for ( final long ele1 : list1 ) {
            while ( ele2 <= ele1 ) {
                if ( ele1 == ele2 ) result += 1;
                if ( !itr2.hasNext() ) return result;
                ele2 = itr2.next();
            }
        }
        return result;
    }
*/

    private static Map<Contig, String> nameContigs( final List<ContigImpl> contigs ) {
        final Map<Contig, String> contigNames = new HashMap<>(contigs.size() * 3);
        int id = 0;
        for ( final ContigImpl contig : contigs ) {
            final String contigName = "c" + ++id;
            contigNames.put( contig, contigName);
            contigNames.put( contig.rc(), contigName + "RC");
        }
        return contigNames;
    }

    private void pathReadsPass( final KmerSet<KmerAdjacency> kmerAdjacencySet,
                                final List<Path> paths,
                                final Map<GapFill, List<List<PathPart>>> gapFillCountMap ) {
        forEachRead( (read, ref, feature, nReadsProcessed) -> {
            final Path path = new Path(read.getBasesNoCopy(), read.getBaseQualitiesNoCopy(), kmerAdjacencySet);
            paths.add(path);
            final List<PathPart> parts = path.getParts();
            final int nParts = parts.size();
            for ( int idx = 1; idx < nParts - 1; ++idx ) {
                final PathPart pathPart = parts.get(idx);
                if ( pathPart.isGap() ) {
                    final Contig start = parts.get(idx-1).getContig();
                    final Contig end = parts.get(idx+1).getContig();
                    final String seq1 = start.canonical().getSequence().toString();
                    final String seq2 = end.canonical().getSequence().toString();
                    if ( seq1.compareTo(seq2) <= 0 ) {
                        final GapFill gapFill = new GapFill(start, end, pathPart.getLength());
                        final List<PathPart> gapParts = parts.subList(idx - 1, idx + 1);
                        gapFillCountMap.computeIfAbsent(gapFill, k -> new ArrayList<>()).add(gapParts);
                    } else {
                        final GapFill gapFill = new GapFill(end.rc(), start.rc(), pathPart.getLength());
                        final List<PathPart> gapParts = parts.subList(idx - 1, idx + 1);
                        gapFillCountMap.computeIfAbsent(gapFill, k -> new ArrayList<>()).add(gapParts);
                    }
                }
            }
        });
    }

    private static void fillGaps( final List<ContigImpl> contigs,
                                  final Map<GapFill, List<List<PathPart>>> gapFillCountMap,
                                  final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
        for ( final Map.Entry<GapFill, List<List<PathPart>>> entry : gapFillCountMap.entrySet() ) {
            final GapFill gapFill = entry.getKey();
            final List<List<PathPart>> gapParts = entry.getValue();
            final int gapSize = gapFill.getDistance();
            final Contig start = gapFill.getStart();
            final Contig end = gapFill.getEnd();
            final int count = gapParts.size();
            if ( count >= 3 && gapSize <= Kmer.KSIZE/2 ) {
                final CharSequence seq1 = start.getSequence();
                final int seq1Size = seq1.length();
                final CharSequence seq2 = end.getSequence();
                final int seq2Start = Kmer.KSIZE - gapSize - 1;
                final String sequence = seq1.subSequence(seq1Size - Kmer.KSIZE + 1, seq1Size).toString() +
                        seq2.subSequence(seq2Start, seq2Start + gapSize);
                final int seqLen = sequence.length();
                KmerAdjacency firstAdjacency = null;
                KmerAdjacency prevAdjacency = null;
                KmerAdjacency curAdjacency = start.getLastKmer();
                KmerAdjacency nextAdjacency;
                int callCount = 0;
                long kVal = 0;
                for ( int idx = 0; idx != seqLen; ++idx ) {
                    kVal <<= 2;
                    switch ( sequence.charAt(idx) ) {
                        case 'C': kVal += 1; break;
                        case 'G': kVal += 2; break;
                        case 'T': kVal += 3; break;
                    }
                    if ( ++callCount >= Kmer.KSIZE ) {
                        nextAdjacency = KmerAdjacency.findOrAdd(kVal, kmerAdjacencySet);
                        if ( callCount == Kmer.KSIZE ) firstAdjacency = nextAdjacency;
                        curAdjacency.observe(prevAdjacency, nextAdjacency, count);
                        prevAdjacency = curAdjacency;
                        curAdjacency = nextAdjacency;
                    }
                }
                curAdjacency.observe(prevAdjacency, end.getFirstKmer(), count);
                final ContigImpl gapTig = new ContigImpl(sequence, count, start, end, firstAdjacency, curAdjacency);
                contigs.add(gapTig);
                start.getSuccessors().add(gapTig);
                end.getPredecessors().add(gapTig);
                final byte[] quals = new byte[sequence.length()];
                Arrays.fill(quals, QMIN);
                KmerAdjacencyImpl.kmerize(sequence.getBytes(), quals, QMIN, kmerAdjacencySet);
                for ( final List<PathPart> parts : gapParts ) {
                    if ( parts.get(0).getContig() == start ) {
                        parts.set(1, new PathPart(gapTig, 0, gapSize));
                    } else {
                        parts.set(1, new PathPart(gapTig.rc(), 0, gapSize));
                    }
                }
            }
        }
    }

    private static void writeDOT( final List<ContigImpl> contigs,
                                  final String fileName ) {
        try ( final BufferedWriter writer = new BufferedWriter(new FileWriter(fileName)) ) {
            writer.write("digraph {\n");
            for ( final Contig contig : contigs ) {
                final double width = contig.getSequence().length() / 100.;
                writer.write(contig + " [width=" + width + "]\n");
                writer.write( contig.rc() + " [width=" + width + "]\n");
            }
            for ( final Contig contig : contigs ) {
                for ( final Contig predecessor : contig.getPredecessors() ) {
                    final String predecessorName = predecessor.rc().toString();
                    writer.write(contig.rc() + " -> " + predecessorName + "\n");
                }
                for ( final Contig successor : contig.getSuccessors() ) {
                    final String successorName = successor.toString();
                    writer.write(contig + " -> " + successorName + "\n");
                }
            }
            writer.write("}\n");
        } catch ( final IOException ioe ) {
            throw new GATKException("Failed to write assembly DOT file.", ioe);
        }
    }

    private static void writeContigs( final List<ContigImpl> contigs ) {
        for ( final Contig contig : contigs ) {
            final List<Contig> predecessors = contig.getPredecessors();
            final String predecessorDescription;
            if ( predecessors.size() == 0 ) {
                predecessorDescription = "\tnone";
            } else {
                final StringBuilder sb = new StringBuilder();
                char prefix = '\t';
                for ( final Contig predecessor : predecessors ) {
                    sb.append(prefix);
                    prefix = ',';
                    sb.append(predecessor);
                }
                predecessorDescription = sb.toString();
            }

            final List<Contig> successors = contig.getSuccessors();
            final String successorDescription;
            if ( successors.size() == 0 ) {
                successorDescription = "\tnone";
            } else {
                final StringBuilder sb = new StringBuilder();
                char prefix = '\t';
                for ( final Contig successor : successors ) {
                    sb.append(prefix);
                    prefix = ',';
                    sb.append(successor);
                }
                successorDescription = sb.toString();
            }

            final String contigName = contig.toString();
            final String component = (contig.isCyclic() ? "(C)\t" : "\t") + contig.getComponentId();
            System.out.println(
                    contigName + component + predecessorDescription + successorDescription + "\t" +
                            contig.getMaxObservations() + "\t" +
                            contig.getFirstKmer().getNObservations() + "\t" +
                            contig.getLastKmer().getNObservations() + "\t" +
                            contig.getSequence().length() + "\t" +
                            contig.getSequence());
        }
    }

    private static void writePaths( final List<Path> readPaths ) {
        final int nReads = readPaths.size();
        for ( int readId = 0; readId != nReads; ++readId ) {
            final Path path = readPaths.get(readId);
            final String pathDesc = path.toString();
            final int nErrors = path.getErrors().size();
            if ( nErrors == 0 ) {
                System.out.println((readId + 1) + ": " + pathDesc);
            } else {
                System.out.println((readId + 1) + ": " + pathDesc + " with " + nErrors + " errors");
            }
        }
    }

    private static void writeCycles( final List<List<Contig>> cycles ) {
        for ( final List<Contig> cycle : cycles ) {
            final StringBuilder sb = new StringBuilder();
            String prefix = "Cycle: ";
            for ( final Contig contig : cycle ) {
                sb.append(prefix).append(contig);
                prefix = ", ";
            }
            System.out.println(sb);
        }
    }

    public static class Kmer {
        public static final int KSIZE = 31;
        public static final long KMASK = (1L << 2*KSIZE) - 1L;
        private final long kVal;

        public Kmer( final long kVal ) { this.kVal = kVal; }

        public long getKVal() { return kVal; }
        public boolean isCanonical() { return isCanonical(kVal); }
        public int getInitialCall() { return (int)(kVal >> (KSIZE*2 - 2)) & 3; }
        public int getFinalCall() { return (int)kVal & 3; }

        public long getPredecessorVal( final int call ) { return (kVal >> 2) | ((long)call << (2 * (KSIZE - 1))); }
        public long getSuccessorVal( final int call ) { return ((kVal << 2) & KMASK) | call; }

        public static boolean isCanonical( final long val ) {
            return (val & (1L << KSIZE)) == 0L;
        }
    }

    public static final class KmerSet<KMER extends Kmer> implements Iterable<KMER> {
        private int capacity;
        private int size;
        // unused buckets contain null.  (this data structure does not support null entries.)
        // if the bucket is unused, the corresponding status byte is irrelevant, but is always set to 0.
        private KMER[] buckets;
        // format of the status bytes:
        // high bit set indicates that the bucket contains a "chain head" (i.e., an entry that naturally belongs in the
        // corresponding bucket).  high bit not set indicates a "squatter" (i.e., an entry that got placed here through the
        // collision resolution methodology).  we use Byte.MIN_VALUE (i.e., 0x80) to pick off this bit.
        // low 7 bits give the (unsigned) offset from the current entry to the next entry in the collision resolution chain.
        // if the low 7 bits are 0, then we'd be pointing at ourselves, which is nonsense, so that particular value marks
        // "end of chain" instead.  we use Byte.MAX_VALUE (i.e., 0x7f) to pick off these bits.
        private byte[] status;

        private static final double LOAD_FACTOR = .85;
        private static final int SPREADER = 241;

        public KmerSet( final int capacity ) {
            this.capacity = computeCapacity(capacity);
            this.size = 0;
            this.buckets = makeBuckets(this.capacity);
            this.status = new byte[this.capacity];
        }

        public final void clear() {
            Arrays.fill(buckets, null);
            Arrays.fill(status, (byte)0);
            size = 0;
        }

        public Iterator<KMER> iterator() { return new Itr(); }

        public KMER find( final long kVal ) {
            return findOrAdd(kVal, k -> null);
        }

        public KMER findOrAdd( final long kVal, final LongFunction<KMER> producer  ) {
            try {
                return findOrAddInternal(kVal, producer);
            } catch ( final HopscotchException he ) {
                resize();
                return findOrAddInternal(kVal, producer);
            }
        }

        private KMER findOrAddInternal( final long kVal, final LongFunction<KMER> producer ) {
            final int bucketIndex = bucketForKVal(kVal);
            if ( !isChainHead(bucketIndex) ) {
                final KMER entry = producer.apply(kVal);
                if ( entry != null ) {
                    insert(entry, bucketIndex);
                }
                return entry;
            }
            KMER entry = buckets[bucketIndex];
            if ( kVal == entry.getKVal() ) {
                return entry;
            }
            int offset;
            int chainIndex = bucketIndex;
            while ( (offset = getOffset(chainIndex)) != 0 ) {
                chainIndex = getIndex(chainIndex, offset);
                entry = buckets[chainIndex];
                if ( kVal == entry.getKVal() ) {
                    return entry;
                }
            }
            entry = producer.apply(kVal);
            if ( entry != null ) {
                append(entry, bucketIndex, chainIndex);
            }
            return entry;
        }

        private void insert( final KMER entry, final int bucketIndex ) {
            if ( buckets[bucketIndex] != null ) evict(bucketIndex);
            buckets[bucketIndex] = entry;
            status[bucketIndex] = Byte.MIN_VALUE;
            size += 1;
        }

        private void append( final KMER entry, final int bucketIndex, final int endOfChainIndex ) {
            final int offsetToEndOfChain = getIndexDiff(bucketIndex, endOfChainIndex);

            // find an empty bucket for the new entry
            int emptyBucketIndex = findEmptyBucket(bucketIndex);

            // if the distance to the empty bucket is larger than this, we'll have to hopscotch
            final int maxOffset = offsetToEndOfChain + Byte.MAX_VALUE;

            // hopscotch the empty bucket into range if it's too far away
            int offsetToEmpty;
            while ( (offsetToEmpty = getIndexDiff(bucketIndex, emptyBucketIndex)) > maxOffset ) {
                emptyBucketIndex = hopscotch(bucketIndex, emptyBucketIndex);
            }

            // if the new entry lies downstream of the current chain end, just link it in
            if ( offsetToEmpty > offsetToEndOfChain ) {
                status[endOfChainIndex] += offsetToEmpty - offsetToEndOfChain;
            } else {
                linkIntoChain(bucketIndex, emptyBucketIndex);
            }
            buckets[emptyBucketIndex] = entry;
            size += 1;
        }

        private void evict( final int bucketToEvictIndex ) {
            final int bucketIndex = bucketForKVal(buckets[bucketToEvictIndex].getKVal());
            final int offsetToEvictee = getIndexDiff(bucketIndex, bucketToEvictIndex);
            int emptyBucketIndex = findEmptyBucket(bucketIndex);
            int fromIndex = bucketIndex;
            while ( true ) {
                while ( getIndexDiff(bucketIndex, emptyBucketIndex) > offsetToEvictee ) {
                    emptyBucketIndex = hopscotch(fromIndex, emptyBucketIndex);
                }
                if ( emptyBucketIndex == bucketToEvictIndex ) return;
                fromIndex = emptyBucketIndex;
                linkIntoChain(bucketIndex, emptyBucketIndex);
                int prevIndex = bucketIndex;
                int offsetToNext = getOffset(prevIndex);
                int nextIndex = getIndex(prevIndex, offsetToNext);
                while ( (offsetToNext = getOffset(nextIndex)) != 0 ) {
                    prevIndex = nextIndex;
                    nextIndex = getIndex(nextIndex, offsetToNext);
                }
                buckets[emptyBucketIndex] = buckets[nextIndex];
                buckets[nextIndex] = null;
                status[nextIndex] = 0;
                status[prevIndex] -= getOffset(prevIndex);
                emptyBucketIndex = nextIndex;
            }
        }

        private int bucketForKVal( final long kVal ) {
            int bucketIndex = (int)((kVal * SPREADER) % capacity);
            if ( bucketIndex < 0 ) {
                bucketIndex += capacity;
            }
            return bucketIndex;
        }

        private int findEmptyBucket( int bucketIndex ) {
            do {
                bucketIndex = getIndex(bucketIndex, 1);
            }
            while ( buckets[bucketIndex] != null );
            return bucketIndex;
        }

        // walk the chain until we find where the new slot gets linked in
        private void linkIntoChain( final int bucketIndex, final int emptyBucketIndex ) {
            int offsetToEmpty = getIndexDiff(bucketIndex, emptyBucketIndex);
            int tmpIndex = bucketIndex;
            int offset;
            while ( (offset = getOffset(tmpIndex)) < offsetToEmpty ) {
                tmpIndex = getIndex(tmpIndex, offset);
                offsetToEmpty -= offset;
            }
            offset -= offsetToEmpty;
            status[tmpIndex] -= offset;
            status[emptyBucketIndex] = (byte) offset;
        }

        private boolean isChainHead( final int bucketIndex ) {
            return (status[bucketIndex] & Byte.MIN_VALUE) != 0;
        }

        private int getOffset( final int bucketIndex ) {
            return status[bucketIndex] & Byte.MAX_VALUE;
        }

        private int getIndex( final int bucketIndex, final int offset ) {
            int result = bucketIndex + offset;
            if ( result >= capacity ) result -= capacity;
            else if ( result < 0 ) result += capacity;
            return result;
        }

        // bucket1 is assumed to be upstream of bucket2 (even if bucket2's index has wrapped)
        // i.e., the result is always positive
        private int getIndexDiff( final int bucketIndex1, final int bucketIndex2 ) {
            int result = bucketIndex2 - bucketIndex1;
            if ( result < 0 ) result += capacity;
            return result;
        }

        private int hopscotch( final int fromIndex, final int emptyBucketIndex ) {
            final int fromToEmptyDistance = getIndexDiff(fromIndex, emptyBucketIndex);
            int offsetToEmpty = Byte.MAX_VALUE;
            while ( offsetToEmpty > 1 ) {
                final int bucketIndex = getIndex(emptyBucketIndex, -offsetToEmpty);
                final int offsetInBucket = getOffset(bucketIndex);
                if ( offsetInBucket != 0 &&
                        offsetInBucket < offsetToEmpty &&
                        offsetToEmpty-offsetInBucket < fromToEmptyDistance ) {
                    final int bucketToMoveIndex = getIndex(bucketIndex, offsetInBucket);
                    move(bucketIndex, bucketToMoveIndex, emptyBucketIndex);
                    return bucketToMoveIndex;
                }
                offsetToEmpty -= 1;
            }
            // this happens now and then, but is usually caught and remedied by a resize
            throw new HopscotchException("Hopscotching failed at load factor "+(1.*size/capacity));
        }


        private void move( int predecessorBucketIndex, final int bucketToMoveIndex, final int emptyBucketIndex ) {
            int toEmptyDistance = getIndexDiff(bucketToMoveIndex, emptyBucketIndex);
            int nextOffset = getOffset(bucketToMoveIndex);
            if ( nextOffset == 0 || nextOffset > toEmptyDistance ) {
                status[predecessorBucketIndex] += toEmptyDistance;
            } else {
                status[predecessorBucketIndex] += nextOffset;
                toEmptyDistance -= nextOffset;
                predecessorBucketIndex = getIndex(bucketToMoveIndex, nextOffset);
                while ( (nextOffset = getOffset(predecessorBucketIndex)) != 0 && nextOffset < toEmptyDistance ) {
                    toEmptyDistance -= nextOffset;
                    predecessorBucketIndex = getIndex(predecessorBucketIndex, nextOffset);
                }
                status[predecessorBucketIndex] = (byte) toEmptyDistance;
            }
            if ( nextOffset != 0 ) {
                status[emptyBucketIndex] = (byte) (nextOffset - toEmptyDistance);
            }
            buckets[emptyBucketIndex] = buckets[bucketToMoveIndex];
            buckets[bucketToMoveIndex] = null;
            status[bucketToMoveIndex] = 0;
        }

        private void resize() {
            final int oldCapacity = capacity;
            final int oldSize = size;
            final KMER[] oldBuckets = buckets;
            final byte[] oldStatus = status;

            capacity = SetSizeUtils.getLegalSizeAbove(capacity);
            size = 0;
            buckets = makeBuckets(capacity);
            status = new byte[capacity];

            try {
                int idx = 0;
                do {
                    final KMER entry = oldBuckets[idx];
                    if ( entry != null ) add(entry);
                }
                while ( (idx = (idx+127)%oldCapacity) != 0 );
            } catch ( final IllegalStateException ise ) {
                capacity = oldCapacity;
                size = oldSize;
                buckets = oldBuckets;
                status = oldStatus;
                // this shouldn't happen except in the case of really bad hashCode implementations
                throw new IllegalStateException("Hopscotching failed at load factor "+1.*size/capacity+", and resizing didn't help.");
            }

            if ( size != oldSize ) {
                // this should never happen, period.
                throw new IllegalStateException("Lost some elements during resizing.");
            }
        }

        private void add( final KMER entry ) {
            final int bucketIndex = bucketForKVal(entry.getKVal());

            // if there's a squatter where the new entry should go, move it elsewhere and put the entry there
            if ( buckets[bucketIndex] != null && !isChainHead(bucketIndex) ) evict(bucketIndex);

            // if the place where it should go is empty, just put the new entry there
            if ( buckets[bucketIndex] == null ) {
                buckets[bucketIndex] = entry;
                status[bucketIndex] = Byte.MIN_VALUE;
                size += 1;
                return;
            }

            // walk to end of chain
            // along the way, make sure the entry isn't already present if necessary
            int endOfChainIndex = bucketIndex;
            int offset;
            while ( (offset = getOffset(endOfChainIndex)) != 0 ) {
                endOfChainIndex = getIndex(endOfChainIndex, offset);
            }

            append(entry, bucketIndex, endOfChainIndex);
        }

        @SuppressWarnings("unchecked")
        private KMER[] makeBuckets( final int size ) {
            return (KMER[])new Kmer[size];
        }

        private static int computeCapacity( final int size ) {
            if ( size < LOAD_FACTOR*Integer.MAX_VALUE ) {
                final int augmentedSize = (int) (size / LOAD_FACTOR);
                for ( final int legalSize : SetSizeUtils.legalSizes ) {
                    if ( legalSize >= augmentedSize ) return legalSize;
                }
            }
            return SetSizeUtils.legalSizes[SetSizeUtils.legalSizes.length - 1];
        }

        private static final class HopscotchException extends IllegalStateException {
            private static final long serialVersionUID = 1L;
            public HopscotchException( final String message ) { super(message); }
        }

        private final class Itr implements Iterator<KMER> {
            private int bucketsIndex = -1;
            private KMER next = null;

            public Itr() { advance(); }

            public boolean hasNext() { return next != null; }

            public KMER next() {
                final KMER result = next;
                if ( result == null ) throw new NoSuchElementException("iterator is exhausted");
                advance();
                return result;
            }

            private void advance() {
                next = null;
                while ( next == null && ++bucketsIndex < buckets.length ) {
                    next = buckets[bucketsIndex];
                }
            }
        }
    }

    public static abstract class KmerAdjacency extends Kmer {
        public KmerAdjacency( final long kVal ) { super(kVal); }

        public abstract KmerAdjacency getSolePredecessor();
        public abstract int getPredecessorMask();
        public abstract int getPredecessorCount();
        public abstract void removeAllPredecessors();
        public abstract void removePredecessor( final int callToRemove, final KmerSet<KmerAdjacency> kmerAdjacencySet );

        public abstract KmerAdjacency getSoleSuccessor();
        public abstract int getSuccessorMask();
        public abstract int getSuccessorCount();
        public abstract void removeAllSuccessors();
        public abstract void removeSuccessor( final int callToRemove, final KmerSet<KmerAdjacency> kmerAdjacencySet );

        public abstract Contig getContig();
        public abstract int getContigOffset();
        public abstract void setContig( final Contig contig, final int contigOffset );

        public abstract int getNObservations();
        public abstract KmerAdjacency rc();

        public void observe( final KmerAdjacency predecessor, final KmerAdjacency successor ) {
            observe(predecessor, successor, 1);
        }

        public abstract void observe( final KmerAdjacency predecessor, final KmerAdjacency successor, final int count );

        @Override public String toString() {
            final StringBuilder sb = new StringBuilder(KSIZE);
            long currentVal = getKVal();
            for ( int idx = 0; idx != KSIZE; ++idx ) {
                sb.append("ACGT".charAt((int)currentVal & 3));
                currentVal >>= 2;
            }
            sb.reverse();
            return sb.toString();
        }

        // Lookup table for reverse-complementing each possible byte value.
        // Each pair of bits represents a base, so you have to reverse bits pairwise and then invert all bits.
        // This is most quickly and easily done with a lookup table.
        private static final long[] BYTEWISE_REVERSE_COMPLEMENT;
        static {
            BYTEWISE_REVERSE_COMPLEMENT = new long[256];
            for ( int bIn = 0; bIn != 256; ++bIn ) {
                BYTEWISE_REVERSE_COMPLEMENT[bIn] =
                        ~(((bIn & 3) << 6) | (((bIn >> 2) & 3) << 4) | (((bIn >> 4) & 3) << 2) | ((bIn >> 6) & 3)) & 0xffL;
            }
        }

        public static long reverseComplement( long val ) {
            // process val one byte at a time
            long result = BYTEWISE_REVERSE_COMPLEMENT[(int)val & 0xFF]; // handle the low-order byte
            int nBytes = 8;
            while ( --nBytes != 0 ) { // pre-decrementing:  we'll go through the loop 7 times
                // rotate down by a byte
                val >>= 8;
                // rotate up by a byte and OR in the reverse complement of the next byte
                result = (result << 8) | BYTEWISE_REVERSE_COMPLEMENT[(int)val & 0xFF];
            }
            return result >>> (Long.SIZE - 2*KSIZE);
        }

        public static KmerAdjacency find( final long kVal, final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            if ( isCanonical(kVal) ) return kmerAdjacencySet.find(kVal & KMASK);
            final KmerAdjacency result = kmerAdjacencySet.find(reverseComplement(kVal));
            return result == null ? null : result.rc();
        }

        public static KmerAdjacency findOrAdd( final long kVal, final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            if ( isCanonical(kVal) ) return kmerAdjacencySet.findOrAdd(kVal & KMASK, KmerAdjacencyImpl::new);
            return kmerAdjacencySet.findOrAdd(reverseComplement(kVal), KmerAdjacencyImpl::new).rc();
        }
    }

    public static final class KmerAdjacencyRC extends KmerAdjacency {
        private final KmerAdjacencyImpl rc;
        private static final int[] NIBREV =
        // 0000,  0001,  0010,  0011,  0100,  0101,  0110,  0111,  1000,  1001,  1010,  1011,  1100,  1101,  1110,  1111
        {0b0000,0b1000,0b0100,0b1100,0b0010,0b1010,0b0110,0b1110,0b0001,0b1001,0b0101,0b1101,0b0011,0b1011,0b0111,0b1111};

        public KmerAdjacencyRC( final KmerAdjacencyImpl rc ) {
            super(reverseComplement(rc.getKVal()));
            this.rc = rc;
        }

        @Override public KmerAdjacency getSolePredecessor() {
            final KmerAdjacency successor = rc.getSoleSuccessor();
            return successor == null ? null : successor.rc();
        }
        @Override public int getPredecessorMask() { return NIBREV[rc.getSuccessorMask()]; }
        @Override public int getPredecessorCount() { return rc.getSuccessorCount(); }
        @Override public void removeAllPredecessors() { rc.removeAllSuccessors(); }
        @Override
        public void removePredecessor( final int callToRemove, final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            rc.removeSuccessor(3 - callToRemove, kmerAdjacencySet);
        }

        @Override public KmerAdjacency getSoleSuccessor() {
            final KmerAdjacency predecessor = rc.getSolePredecessor();
            return predecessor == null ? null : predecessor.rc();
        }
        @Override public int getSuccessorMask() { return NIBREV[rc.getPredecessorMask()]; }
        @Override public int getSuccessorCount() { return rc.getPredecessorCount(); }
        @Override public void removeAllSuccessors() { rc.removeAllPredecessors(); }
        @Override
        public void removeSuccessor( final int callToRemove, final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            rc.removePredecessor(3 - callToRemove, kmerAdjacencySet);
        }

        @Override public Contig getContig() {
            final Contig contig = rc.getContig();
            return contig == null ? null : contig.rc();
        }
        @Override public int getContigOffset() {
            final Contig contig = rc.getContig();
            return contig == null ? 0 : contig.size() - rc.getContigOffset() - KSIZE;
        }
        @Override public void setContig( final Contig contig, final int contigOffset ) {
            if ( contig == null ) rc.setContig(null, 0);
            else rc.setContig(contig.rc(), contig.size() - contigOffset - KSIZE);
        }

        @Override public int getNObservations() { return rc.getNObservations(); }
        @Override public KmerAdjacency rc() { return rc; }

        @Override public void observe( final KmerAdjacency predecessor, final KmerAdjacency successor, final int count ) {
            rc.observe(successor == null ? null : successor.rc(), predecessor == null ? null : predecessor.rc(), count);
        }
    }

    public static final class KmerAdjacencyImpl extends KmerAdjacency {
        private KmerAdjacency solePredecessor = null; // set to null if there are no predecessors, or multiple predecessors
        private KmerAdjacency soleSuccessor = null; // set to null if there are no successors, or multiple successors
        private int predecessorMask = 0; // bit mask of observed kmers preceding this one
        private int successorMask = 0; // bit mask observed kmers following this one
        private Contig contig = null; // the contig that contains this Kmer
        private int contigOffset;
        private int nObservations = 0; // the reads in which this kmer was observed
        private final KmerAdjacencyRC rc; // the reverse-complement of this kmer
        private static final int[] COUNT_FOR_MASK =
                //side sum for binary values from 0 -> 15
                //0000  0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111
                {    0,    1,   1,   2,   1,   2,   2,   3,   1,   2,   2,   3,   2,   3,   3,   4 };

        public KmerAdjacencyImpl( final long kVal ) {
            super(kVal);
            this.rc = new KmerAdjacencyRC(this);
        }

        @Override public KmerAdjacency getSolePredecessor() { return solePredecessor; } // may return null
        @Override public int getPredecessorMask() { return predecessorMask; }
        @Override public int getPredecessorCount() { return COUNT_FOR_MASK[predecessorMask]; }
        @Override public void removeAllPredecessors() { predecessorMask = 0; solePredecessor = null; }
        @Override
        public void removePredecessor( final int callToRemove, final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            predecessorMask &= ~(1 << callToRemove);
            solePredecessor = null;
            if ( getPredecessorCount() == 1 ) {
                for ( int call = 0; call != 4; ++call ) {
                    if ( ((1 << call) & predecessorMask) != 0 ) {
                        solePredecessor = find(getPredecessorVal(call), kmerAdjacencySet);
                        break;
                    }
                }
            }
        }

        @Override public KmerAdjacency getSoleSuccessor() { return soleSuccessor; } // may return null
        @Override public int getSuccessorMask() { return successorMask; }
        @Override public int getSuccessorCount() { return COUNT_FOR_MASK[successorMask]; }
        @Override public void removeAllSuccessors() { successorMask = 0; soleSuccessor = null; }
        @Override
        public void removeSuccessor( final int callToRemove, final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            successorMask &= ~(1 << callToRemove);
            soleSuccessor = null;
            if ( getSuccessorCount() == 1 ) {
                for ( int call = 0; call != 4; ++call ) {
                    if ( ((1 << call) & successorMask) != 0 ) {
                        soleSuccessor = find(getSuccessorVal(call), kmerAdjacencySet);
                        break;
                    }
                }
            }
        }

        @Override public Contig getContig() { return contig; }
        @Override public int getContigOffset() { return contigOffset; }
        @Override public void setContig( final Contig contig, final int contigOffset ) {
            this.contig = contig;
            this.contigOffset = contigOffset;
        }

        @Override public int getNObservations() { return nObservations; }
        @Override public KmerAdjacency rc() { return rc; }

        @Override public void observe( final KmerAdjacency predecessor,
                                       final KmerAdjacency successor,
                                       final int count ) {
            if ( predecessor != null ) {
                final int initialCall = predecessor.getInitialCall();
                final int newPredecessorMask = 1 << initialCall;
                if ( (newPredecessorMask & predecessorMask) == 0 ) {
                    if ( predecessorMask == 0 ) {
                        solePredecessor = predecessor;
                        predecessorMask = newPredecessorMask;
                    } else {
                        solePredecessor = null;
                        predecessorMask |= newPredecessorMask;
                    }
                }
            }
            if ( successor != null ) {
                final int finalCall = successor.getFinalCall();
                final int newSuccessorMask = 1 << finalCall;
                if ( (newSuccessorMask & successorMask) == 0 ) {
                    if ( successorMask == 0 ) {
                        soleSuccessor = successor;
                        successorMask = newSuccessorMask;
                    } else {
                        soleSuccessor = null;
                        successorMask |= newSuccessorMask;
                    }
                }
            }
            nObservations += count;
        }

        public static void kmerize( final byte[] calls,
                                    final byte[] quals,
                                    final byte qMin,
                                    final KmerSet<KmerAdjacency> kmerSet ) {
            int currentCount = 0;
            long currentKVal = 0;
            KmerAdjacency prevAdjacency = null;
            KmerAdjacency currentAdjacency = null;
            for ( int idx = 0; idx < calls.length; ++idx ) {
                if ( quals[idx] <  qMin ) {
                    if ( currentAdjacency != null ) {
                        currentAdjacency.observe(prevAdjacency, null);
                    }
                    currentCount = 0;
                    currentAdjacency = prevAdjacency = null;
                    continue;
                }
                currentKVal <<= 2;
                switch ( calls[idx] ) {
                    case 'A': case 'a': break;
                    case 'C': case 'c': currentKVal += 1; break;
                    case 'G': case 'g': currentKVal += 2; break;
                    case 'T': case 't': currentKVal += 3; break;
                    default:
                        if ( currentAdjacency != null ) {
                            currentAdjacency.observe(prevAdjacency, null);
                        }
                        currentCount = 0;
                        currentAdjacency = prevAdjacency = null;
                        continue;
                }
                if ( ++currentCount >= KSIZE ) {
                    final KmerAdjacency nextAdjacency = findOrAdd(currentKVal, kmerSet);
                    if ( currentAdjacency != null ) {
                        currentAdjacency.observe(prevAdjacency, nextAdjacency);
                    }
                    prevAdjacency = currentAdjacency;
                    currentAdjacency = nextAdjacency;
                }
            }
            if ( currentAdjacency != null ) {
                currentAdjacency.observe(prevAdjacency, null);
            }
        }
    }

    public enum ContigOrientation {
        FWD, // k-mer appears at the 5' end of the contig
        REV, // k-mer appears at the 5' end of the reverse-complemented contig
        BOTH // k-mer occurs on 5' end of the contig and its RC (can happen when the contig is a palindrome)
    }

    public enum DFSearchStatus {
        UNVISITED,
        VISITING,
        VISITED
    }

    public static final class ContigEndKmer extends Kmer {
        private final Contig contig;
        private final ContigOrientation contigOrientation;

        public ContigEndKmer( final long kVal, final Contig contig, final ContigOrientation contigEnd ) {
            super(kVal);
            this.contig = contig;
            this.contigOrientation = contigEnd;
        }

        public Contig getContig() { return contig; }
        public ContigOrientation getContigOrientation() { return contigOrientation; }
    }

    public interface Contig {
        CharSequence getSequence();
        int getMaxObservations();
        KmerAdjacency getFirstKmer();
        KmerAdjacency getLastKmer();
        List<Contig> getPredecessors();
        List<Contig> getSuccessors();
        int getComponentId();
        int size();
        Contig rc();
        boolean isCyclic();
        void setCyclic( final boolean cyclic );
        boolean isCut();
        void setCut( final boolean cut );
        boolean isCanonical();
        ContigImpl canonical();
        Object getAuxData();
        void setAuxData( final Object obj );
    }

    public static final class ContigImpl implements Contig {
        private static int nContigs;
        private final int id;
        private final CharSequence sequence;
        private final int maxObservations;
        private final KmerAdjacency firstKmer;
        private final KmerAdjacency lastKmer;
        private final List<Contig> predecessors;
        private final List<Contig> successors;
        private int componentId;
        private boolean cyclic;
        private boolean cut;
        private final Contig rc;
        private Object auxData;

        public ContigImpl( final CharSequence sequence, final int maxObservations,
                           final Contig predecessor, final Contig successor,
                           final KmerAdjacency firstKmer, final KmerAdjacency lastKmer ) {
            this.id = nContigs++;
            this.sequence = sequence;
            this.maxObservations = maxObservations;
            this.firstKmer = firstKmer;
            this.lastKmer = lastKmer;
            this.predecessors = new ArrayList<>(1);
            predecessors.add(predecessor);
            this.successors = new ArrayList<>(1);
            successors.add(successor);
            this.rc = new ContigRCImpl(this);
        }

        public ContigImpl( final KmerAdjacency firstKmerAdjacency ) {
            this.id = nContigs++;
            final StringBuilder sb = new StringBuilder(firstKmerAdjacency.toString());
            int maxObservations = firstKmerAdjacency.getNObservations();
            KmerAdjacency lastKmerAdjacency = firstKmerAdjacency;
            for ( KmerAdjacency kmerAdjacency = firstKmerAdjacency.getSoleSuccessor();
                  kmerAdjacency != null;
                  kmerAdjacency = kmerAdjacency.getSoleSuccessor() ) {
                // if we've gone around a circle, or if we're branching backwards, or if we hit a palindrome u-turn
                if ( firstKmerAdjacency == kmerAdjacency ||
                        kmerAdjacency.getPredecessorCount() != 1 ||
                        kmerAdjacency == lastKmerAdjacency.rc() ) {
                    break;
                }
                sb.append("ACGT".charAt(kmerAdjacency.getFinalCall()));
                maxObservations = Math.max(maxObservations, kmerAdjacency.getNObservations());
                lastKmerAdjacency = kmerAdjacency;
            }
            this.sequence = sb.toString();
            this.maxObservations = maxObservations;
            this.firstKmer = firstKmerAdjacency;
            this.lastKmer = lastKmerAdjacency;
            this.predecessors = new ArrayList<>(firstKmer.getPredecessorCount());
            this.successors = new ArrayList<>(lastKmer.getSuccessorCount());
            this.rc = new ContigRCImpl(this);

            int offset = 0;
            for ( KmerAdjacency kmerAdjacency = firstKmerAdjacency;
                  kmerAdjacency != lastKmerAdjacency;
                  kmerAdjacency = kmerAdjacency.getSoleSuccessor() ) {
                kmerAdjacency.setContig(this, offset++);
            }
            lastKmerAdjacency.setContig(this, offset);
        }

        // create a new contig by joining two contigs
        public ContigImpl( final Contig predecessor, final Contig successor ) {
            this.id = nContigs++;
            final StringBuilder sb = new StringBuilder(predecessor.getSequence());
            final CharSequence successorSequence = successor.getSequence();
            sb.append(successorSequence.subSequence(Kmer.KSIZE - 1, successorSequence.length()));
            this.sequence = sb.toString();
            this.maxObservations = Math.max(predecessor.getMaxObservations(), successor.getMaxObservations());
            this.firstKmer = predecessor.getFirstKmer();
            this.lastKmer = successor.getLastKmer();
            this.predecessors = new ArrayList<>(predecessor.getPredecessors());
            predecessors.replaceAll( contig -> contig == successor ? this : contig );
            this.successors = new ArrayList<>(successor.getSuccessors());
            successors.replaceAll( contig -> contig == predecessor ? this : contig );
            this.rc = new ContigRCImpl(this);
        }

        public int getId() { return id; }
        @Override public CharSequence getSequence() { return sequence; }
        @Override public int getMaxObservations() { return maxObservations; }
        @Override public KmerAdjacency getFirstKmer() { return firstKmer; }
        @Override public KmerAdjacency getLastKmer() { return lastKmer; }
        @Override public List<Contig> getPredecessors() { return predecessors; }
        @Override public List<Contig> getSuccessors() { return successors; }
        @Override public int getComponentId() { return componentId; }
        public void setComponentId( final int id ) { this.componentId = id; }
        @Override public int size() { return sequence.length(); }
        @Override public Contig rc() { return rc; }
        @Override public boolean isCyclic() { return cyclic; }
        @Override public void setCyclic( final boolean cyclic ) { this.cyclic = cyclic; }
        @Override public boolean isCut() { return cut; }
        @Override public void setCut( final boolean cut ) { this.cut = cut; }
        @Override public boolean isCanonical() { return true; }
        @Override public ContigImpl canonical() { return this; }
        @Override public Object getAuxData() { return auxData; }
        @Override public void setAuxData( final Object auxData ) { this.auxData = auxData; }
        @Override public String toString() { return "c" + Integer.toString(id); }
    }

    public static final class ContigRCImpl implements Contig {
        private final CharSequence sequence;
        private final List<Contig> predecessors;
        private final List<Contig> successors;
        private final ContigImpl rc;
        private Object auxData;

        public ContigRCImpl( final ContigImpl contig ) {
            this.sequence = new SequenceRC(contig.getSequence());
            this.predecessors = new ListRC(contig.getSuccessors());
            this.successors = new ListRC(contig.getPredecessors());
            this.rc = contig;
        }

        @Override public CharSequence getSequence() { return sequence; }
        @Override public int getMaxObservations() { return rc.getMaxObservations(); }
        @Override public KmerAdjacency getFirstKmer() { return rc.getLastKmer().rc(); }
        @Override public KmerAdjacency getLastKmer() { return rc.getFirstKmer().rc(); }
        @Override public List<Contig> getPredecessors() { return predecessors; }
        @Override public List<Contig> getSuccessors() { return successors; }
        @Override public int getComponentId() { return rc.getComponentId(); }
        @Override public int size() { return sequence.length(); }
        @Override public Contig rc() { return rc; }
        @Override public boolean isCyclic() { return rc.isCyclic(); }
        @Override public void setCyclic( final boolean cyclic ) { rc.setCyclic(cyclic); }
        @Override public boolean isCut() { return rc.isCut(); }
        @Override public void setCut( final boolean cut ) { rc.setCut(cut); }
        @Override public boolean isCanonical() { return false; }
        @Override public ContigImpl canonical() { return rc; }
        @Override public Object getAuxData() { return auxData; }
        @Override public void setAuxData( final Object auxData ) { this.auxData = auxData; }
        @Override public String toString() { return rc.toString() + "RC"; }

        public static final class SequenceRC implements CharSequence {
            private final int lenLess1;
            private final CharSequence sequence;

            public SequenceRC( final CharSequence sequence ) {
                this.lenLess1 = sequence.length() - 1;
                this.sequence = sequence;
            }

            @Override public int length() { return sequence.length(); }
            @Override public char charAt( final int index ) {
                final char result;
                switch ( sequence.charAt(lenLess1 - index) ) {
                    case 'A': result = 'T'; break;
                    case 'C': result = 'G'; break;
                    case 'G': result = 'C'; break;
                    case 'T': result = 'A'; break;
                    default: result = 'N'; break;
                }
                return result;
            }
            @Override public CharSequence subSequence( final int start, final int end ) {
                return new StringBuilder(end - start).append(this, start, end);
            }
            @Override public String toString() { return new StringBuilder(this).toString(); }
        }

        public static final class ListRC extends AbstractList<Contig> {
            private final List<Contig> contigList;

            public ListRC( final List<Contig> contigList ) {
                this.contigList = contigList;
            }

            @Override public Contig get( final int index ) { return contigList.get(index).rc(); }
            @Override public int size() { return contigList.size(); }
            @Override public Contig set( final int index, final Contig contig ) {
                return contigList.set(index, contig.rc()).rc();
            }
            @Override public void add( final int index, final Contig contig ) { contigList.add(index, contig.rc()); }
            @Override public Contig remove( final int index ) { return contigList.remove(index).rc(); }
        }
    }

    public static final class PathPart {
        private final Contig contig;
        private final int start;
        private int stop;

        public PathPart() { this(null, 0, 1); }
        public PathPart( final Contig contig, final int start ) { this(contig, start, start+1); }
        public PathPart( final Contig contig, final int start, final int stop ) {
            this.contig = contig;
            this.start = start;
            this.stop = stop;
        }

        public Contig getContig() { return contig; }
        public int getStart() { return start; }
        public int getStop() { return stop; }
        public void setStop( final int stop ) { this.stop = stop; }
        public boolean isGap() { return contig == null; }
        public int getLength() { return stop - start; }

        public void extendPath() { stop += 1; }
        public PathPart rc() {
            if ( contig == null ) return this;
            final int revBase = contig.size() - Kmer.KSIZE + 1;
            return new PathPart(contig.rc(), revBase - stop, revBase - start);
        }
    }

    public static final class Error {
        private final ContigImpl contig;
        private final int offset;
        private final byte call;
        private final byte quality;

        public Error( final Contig contig, final int offset, final byte call, final byte quality ) {
            this.contig = contig.canonical();
            this.offset = this.contig == contig ? offset : contig.size() - offset - 1;
            this.call = call;
            this.quality = quality;
        }

        public Contig getContig() { return contig; }
        public int getOffset() { return offset; }
        public byte getCall() { return call; }
        public byte getQuality() { return quality; }
    }

    public static final class GapFill {
        private final Contig start;
        private final Contig end;
        private final int distance;

        public GapFill( final Contig start, final Contig end, final int distance ) {
            this.start = start;
            this.end = end;
            this.distance = distance;
        }

        public Contig getStart() { return start; }
        public Contig getEnd() { return end; }
        public int getDistance() { return distance; }

        @Override public int hashCode() {
            return 47 * (47 * (47 * start.hashCode() + end.hashCode()) + distance);
        }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof GapFill && equals((GapFill)obj);
        }

        public boolean equals( final GapFill that ) {
            return this.start == that.start && this.end == that.end && this.distance == that.distance;
        }
    }

    public static final class Path {
        private final List<PathPart> parts;
        private final List<Error> errors;

        // RCing constructor
        private Path( final Path that ) {
            this.parts = new ArrayList<>();
            final List<PathPart> thoseParts = that.parts;
            for ( int idx = thoseParts.size() - 1; idx >= 0; --idx ) {
                parts.add(thoseParts.get(idx).rc());
            }
            this.errors = that.errors;
        }

        public Path( final byte[] readCalls,
                     final byte[] quals,
                     final KmerSet<KmerAdjacency> kmerAdjacencySet ) {
            parts = new ArrayList<>();
            List<Error> errs = null;
            byte[] calls = readCalls;
            long kVal = 0;
            int count = 0;
            PathPart currentPathPart = null;
            for ( int idx = 0; idx != calls.length; ++idx ) {
                final byte call = calls[idx];
                kVal <<= 2;
                switch ( call ) {
                    case 'C': case 'c': kVal += 1; break;
                    case 'G': case 'g': kVal += 2; break;
                    case 'T': case 't': kVal += 3; break;
                    case 'N': case 'n':
                        if ( readCalls == calls ) {
                            calls = Arrays.copyOf(readCalls, readCalls.length);
                        }
                        calls[idx] = 'A';
                        break;
                }
                if ( ++count >= Kmer.KSIZE ) {
                    final KmerAdjacency kmer = KmerAdjacencyImpl.find(kVal, kmerAdjacencySet);
                    Contig contig;
                    final int contigOffset;
                    // if we fail to look up the kmer (or if it's a suppressed kmer with no contig)
                    if ( kmer == null || (contig = kmer.getContig()) == null ) {
                        if ( currentPathPart == null ) {
                            // if there's no current path part, just create the 1st one as a NoKmer path part
                            // we'll try to backtrack if we run into a good kmer
                            currentPathPart = new PathPart();
                            parts.add(currentPathPart);
                        } else if ( currentPathPart.isGap() ) {
                            // if the current path part is NoKmer, just extend it
                            currentPathPart.extendPath();
                        } else if ( (contigOffset = currentPathPart.getStop() + Kmer.KSIZE -1) <
                                (contig = currentPathPart.getContig()).size() ) {
                            // if the current path part is on some contig, note the mismatch and extend it
                            if ( errs == null ) errs = new ArrayList<>();
                            errs.add(new Error(contig, contigOffset, call, quals[idx]));
                            currentPathPart.extendPath();
                            kVal &= ~3;
                            switch ( contig.getSequence().charAt(contigOffset) ) {
                                case 'C': case 'c': kVal += 1; break;
                                case 'G': case 'g': kVal += 2; break;
                                case 'T': case 't': kVal += 3; break;
                            }
                        } else if ( contig.getSuccessors().size() == 1 ) {
                            // at end of contig, but there's only one choice for successor contig
                            final Contig soleSuccessor = contig.getSuccessors().get(0);
                            if ( errs == null ) errs = new ArrayList<>();
                            errs.add(new Error(soleSuccessor, 0, call, quals[idx]));
                            currentPathPart = new PathPart(soleSuccessor, 0);
                            parts.add(currentPathPart);
                            kVal &= ~3;
                            switch ( soleSuccessor.getSequence().charAt(0) ) {
                                case 'C': case 'c': kVal += 1; break;
                                case 'G': case 'g': kVal += 2; break;
                                case 'T': case 't': kVal += 3; break;
                            }
                        } else {
                            // current path part is at the end of its contig -- create a new NoKmer path part
                            currentPathPart = new PathPart();
                            parts.add(currentPathPart);
                        }
                    } else {
                        // we've found our kmer
                        if ( currentPathPart == null ) {
                            // we've looked up a kmer, but don't have a current path part -- create one
                            currentPathPart = new PathPart(contig, kmer.getContigOffset());
                            parts.add(currentPathPart);
                        } else if ( contig == currentPathPart.getContig() ) {
                            // our lookup is on the current path part's contig -- extend it
                            if ( kmer.getContigOffset() == currentPathPart.getStop() ) {
                                currentPathPart.extendPath();
                            } else {
                                // weird:  kmer is non-contiguous.  start a new path part
                                currentPathPart = new PathPart(contig, kmer.getContigOffset());
                                parts.add(currentPathPart);
                            }
                        } else if ( !currentPathPart.isGap() ) {
                            // we're jumping to a new contig.  start a new path part
                            currentPathPart = new PathPart(contig, kmer.getContigOffset());
                            parts.add(currentPathPart);
                        } else if ( kmer.getContigOffset() == 0 && contig.getPredecessors().size() != 1 ) {
                            // we got our 1st good kmer lookup at the start of a contig after a chunk of NoKmers
                            // just add a new path part for it
                            currentPathPart = new PathPart(contig, 0);
                            parts.add(currentPathPart);
                        } else {
                            // we got our 1st good kmer lookup after a chunk of NoKmers, and we're not at the very start
                            // of the contig, so there's an upstream error to fix.
                            // we don't know how to fix errors in reverse, so rc the chunk in question,
                            // path it in the forward direction recursively, and rc that path.
                            parts.remove( parts.size() - 1);
                            final int end = idx + 1;
                            final int start = end - Kmer.KSIZE - currentPathPart.getStop();
                            final byte[] rcCalls = Arrays.copyOfRange(calls, start, end);
                            SequenceUtil.reverseComplement(rcCalls);
                            final byte[] rQuals = Arrays.copyOfRange(quals, start, end);
                            SequenceUtil.reverseQualities(rQuals);
                            final Path rcPath = new Path(rcCalls, rQuals, kmerAdjacencySet).rc();
                            parts.addAll(rcPath.getParts());
                            currentPathPart = parts.get(parts.size() - 1);
                        }
                    }
                }
            }
            this.errors = errs == null ? Collections.emptyList() : errs;
        }

        public List<PathPart> getParts() { return parts; }
        public List<Error> getErrors() { return errors; }
        public Path rc() { return new Path(this); }

        @Override public String toString() {
            if ( parts.size() == 0 ) return "";
            final StringBuilder sb = new StringBuilder();
            String prefix = "";
            final PathPart firstPart = parts.get(0);
            final PathPart lastPart = parts.get(parts.size() - 1);
            for ( final PathPart pp : parts ) {
                sb.append(prefix);
                prefix = ", ";
                if ( pp.isGap() ) {
                    sb.append("NoKmer(").append(pp.getLength()).append(")");
                } else {
                    final Contig contig = pp.getContig();
                    sb.append(contig);
                    final int maxStop = contig.size() - Kmer.KSIZE + 1;
                    if ( (pp != firstPart && pp.getStart() != 0) ||
                         (pp != lastPart && pp.getStop() != maxStop) ) {
                        sb.append('(').append(pp.getStart()).append('-').append(pp.getStop()).append('/');
                        sb.append(maxStop).append(')');
                    }
                }
            }
            return sb.toString();
        }
    }

    public static final class CutData {
        public static int nextNum;
        public int visitNum;
        public int minVisitNum;

        public CutData() {
            this.visitNum = ++nextNum;
            this.minVisitNum = visitNum;
        }
    }
}
