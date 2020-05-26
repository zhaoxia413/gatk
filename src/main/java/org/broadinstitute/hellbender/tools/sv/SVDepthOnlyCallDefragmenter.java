package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.util.stream.Collectors;

public class SVDepthOnlyCallDefragmenter extends LocatableClusterEngine<SVCallRecordWithEvidence> {

    private final double minSampleOverlap;
    private static final double PADDING_FRACTION = 0.2;
    private final List<GenomeLoc> coverageIntervals;
    private final TreeMap<GenomeLoc, Integer> genomicToBinMap;
    final GenomeLocParser parser;

    public SVDepthOnlyCallDefragmenter(final SAMSequenceDictionary dictionary) {
        this(dictionary, 0.9, null);
    }

    //for single-sample clustering case
    public SVDepthOnlyCallDefragmenter(final SAMSequenceDictionary dictionary, double minSampleOverlap, List<SimpleInterval> traversalIntervals) {
        super(dictionary, CLUSTERING_TYPE.SINGLE_LINKAGE);
        parser = new GenomeLocParser(this.dictionary);
        this.minSampleOverlap = minSampleOverlap;

        if (traversalIntervals != null) {
            this.coverageIntervals = IntervalUtils.genomeLocsFromLocatables(parser,traversalIntervals);
            genomicToBinMap = new TreeMap<>();
            for (int i = 0; i < traversalIntervals.size(); i++) {
                genomicToBinMap.put(parser.createGenomeLoc(traversalIntervals.get(i)),i);
            }
        } else {
            genomicToBinMap = null;
            coverageIntervals = null;
        }
    }

    /**
     * Find a single call representative of all the calls in the {@param cluster}
     * @param cluster   the events that are clustered together
     * @return  a call encompassing all the cluster's events and containing all the algorithms and genotypes
     */
    @Override
    protected SVCallRecordWithEvidence flattenCluster(final Collection<SVCallRecordWithEvidence> cluster) {
        final int newStart =  cluster.stream().mapToInt(SVCallRecordWithEvidence::getStart).min().getAsInt();
        final int newEnd = cluster.stream().mapToInt(SVCallRecordWithEvidence::getEnd).max().getAsInt();
        final SVCallRecordWithEvidence exampleCall = cluster.iterator().next();
        final int length = newEnd - newStart + 1;  //+1 because GATK intervals are inclusive
        final List<String> algorithms = cluster.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList()); //should be depth only
        final List<Genotype> clusterGenotypes = cluster.stream().flatMap(v -> v.getGenotypes().stream()).collect(Collectors.toList());
        return new SVCallRecordWithEvidence(exampleCall.getContig(), newStart, exampleCall.getStartStrand(),
                exampleCall.getEndContig(), newEnd, exampleCall.getEndStrand(), exampleCall.getType(), length, algorithms, clusterGenotypes,
                exampleCall.getStartSplitReadSites(), exampleCall.getEndSplitReadSites(), exampleCall.getDiscordantPairs());
    }

    /**
     * Determine if two calls should cluster based on their padded intervals and genotyped samples
     * @param a
     * @param b
     * @return true if the two calls should be in the same cluster
     */
    @Override
    protected boolean clusterTogether(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        if (!isDepthOnlyCall(a) || !isDepthOnlyCall(b)) return false;
        Utils.validate(a.getContig().equals(a.getEndContig()), "Call A is depth-only but interchromosomal");
        Utils.validate(b.getContig().equals(b.getEndContig()), "Call B is depth-only but interchromosomal");
        if (!a.getType().equals(b.getType())) return false;
        final Set<String> sharedSamples = new LinkedHashSet<>(a.getSamples());
        sharedSamples.retainAll(b.getSamples());
        final double sampleOverlap = Math.min(sharedSamples.size() / (double) a.getSamples().size(), sharedSamples.size() / (double) b.getSamples().size());
        if (sampleOverlap < minSampleOverlap) return false;
        //in the single-sample case, match copy number strictly
        boolean copyNumbersAgree = true;
        if (a.getGenotypes().size() == 1 && b.getGenotypes().size() == 1) {
            if (a.getGenotypes().get(0).hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) && b.getGenotypes().get(0).hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) &&
                a.getGenotypes().get(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) != b.getGenotypes().get(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                copyNumbersAgree = false;
            }
        }
        return getClusteringInterval(a, null)
                .overlaps(getClusteringInterval(b, null)) && copyNumbersAgree;
    }


    /**
     * Determine an overlap interval for clustering using {@value #PADDING_FRACTION} padding
     * Returned interval represents the interval in which the start position of a new event must fall in order to be added to the cluster (including {@param call})
     * @param call  new event to be clustered
     * @param currentClusterInterval    the cluster of interest, may be null
     * @return  an interval describing the cluster after {@param call} is added
     */
    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecordWithEvidence call, final SimpleInterval currentClusterInterval) {
        Utils.nonNull(call);
        final SimpleInterval callInterval = getCallInterval(call);
        final int paddedCallStart, paddedCallEnd;
        if (genomicToBinMap != null) {
            final GenomeLoc callStart = parser.createGenomeLoc(call.getContig(), call.getStart(), call.getStart());
            final GenomeLoc callEnd = parser.createGenomeLoc(call.getContig(), call.getEnd(), call.getEnd());
            final int callStartIndex = genomicToBinMap.ceilingEntry(callStart) == null ? 0 : genomicToBinMap.ceilingEntry(callStart) .getValue();
            final int callEndIndex = genomicToBinMap.floorEntry(callEnd) == null ? genomicToBinMap.size() - 1 : genomicToBinMap.floorEntry(callEnd).getValue();
            final int callBinLength = callEndIndex - callStartIndex + 1;
            final int paddedStartIndex = Math.max(callStartIndex - (int)Math.round(callBinLength * PADDING_FRACTION), 0);
            paddedCallStart = coverageIntervals.get(paddedStartIndex).getStart();
            final int paddedEndIndex = Math.min(callEndIndex + (int)Math.round(callBinLength * PADDING_FRACTION), genomicToBinMap.size() - 1);
            paddedCallEnd = coverageIntervals.get(paddedEndIndex).getEnd();
        } else {
            paddedCallStart = (int) (callInterval.getStart() - PADDING_FRACTION * callInterval.getLengthOnReference());
            paddedCallEnd = (int) (callInterval.getEnd() + PADDING_FRACTION * callInterval.getLengthOnReference());
        }
        final String currentContig = getCurrentContig();
        final int contigLength = dictionary.getSequence(currentContig).getSequenceLength();
        if (currentClusterInterval == null) {
            return IntervalUtils.trimIntervalToContig(currentContig, paddedCallStart, paddedCallEnd, contigLength);
        }
        //NOTE: this is an approximation -- padding should be based on the length of the call plus currentClusterIntervals
        final int newMinStart = Math.min(paddedCallStart, currentClusterInterval.getStart());
        final int newMaxEnd = Math.max(paddedCallEnd, currentClusterInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(currentContig, newMinStart, newMaxEnd, contigLength);
    }

    // Not used for single-linkage clustering
    @Override
    protected boolean itemsAreIdentical(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        return false;
    }

    // Not used for single-linkage clustering
    @Override
    protected SVCallRecordWithEvidence deduplicateIdenticalItems(final Collection<SVCallRecordWithEvidence> items) {
        return null;
    }

    private SimpleInterval getCallInterval(final SVCallRecordWithEvidence call) {
        return new SimpleInterval(call.getContig(), call.getStart(), call.getEnd());
    }

    public static boolean isDepthOnlyCall(final SVCallRecord call) {
        if (call.getAlgorithms().isEmpty()) return false;
        for (final String alg : call.getAlgorithms()) {
            if (!alg.equals(GATKSVVCFConstants.DEPTH_ALGORITHM)) return false;
        }
        return true;
    }

    @VisibleForTesting
    public static double getPaddingFraction() {
        return PADDING_FRACTION;
    }
}
