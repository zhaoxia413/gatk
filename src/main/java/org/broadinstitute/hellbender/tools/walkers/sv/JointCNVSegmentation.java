package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVClusterEngine;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordWithEvidence;
import org.broadinstitute.hellbender.tools.sv.SVDepthOnlyCallDefragmenter;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import shaded.cloud_nio.com.google.errorprone.annotations.Var;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@BetaFeature
@CommandLineProgramProperties(
        summary = "Gathers single-sample segmented gCNV VCFs, harmonizes breakpoints, and outputs a cohort VCF with genotypes.  Note that this tool assumes ploidy 2.",
        oneLineSummary = "Combined single-sample segmented gCNV VCFs.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class JointCNVSegmentation extends MultiVariantWalkerGroupedOnStart {

    private SortedSet<String> samples;
    private VariantContextWriter vcfWriter;
    private SAMSequenceDictionary dictionary;
    private SVDepthOnlyCallDefragmenter defragmenter;
    private SVClusterEngine clusterEngine;
    private List<GenomeLoc> callIntervals;

    private String currentContig;

    public static final String MIN_QUALITY_LONG_NAME = "minimum-qs-score";
    public static final String MODEL_CALL_INTERVALS = "model-call-intervals";

    @Argument(fullName = MIN_QUALITY_LONG_NAME, doc = "Minimum QS score to combine a variant segment")
    private int minQS = 20;

    @Argument(fullName = MODEL_CALL_INTERVALS, doc = "Intervals used for gCNV calls.  Should be preprocessed and filtered to line up with model calls. Required for exomes.")
    private File modelCallIntervalList;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined output file", optional=false)
    private File outputFile;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }

        final GenomeLocParser parser = new GenomeLocParser(this.dictionary);

        if (modelCallIntervalList == null) {
            callIntervals = null;
        } else {
        final List<GenomeLoc> inputCoverageIntervals = IntervalUtils.featureFileToIntervals(parser, modelCallIntervalList.getAbsolutePath());
        final List<GenomeLoc> inputTraversalIntervals = IntervalUtils.genomeLocsFromLocatables(parser,getTraversalIntervals());
            callIntervals = IntervalUtils.mergeListsBySetOperator(inputCoverageIntervals, inputTraversalIntervals, IntervalSetRule.INTERSECTION);
        }

        defragmenter = new SVDepthOnlyCallDefragmenter(dictionary, 0.8, callIntervals);
        clusterEngine = new SVClusterEngine(dictionary, true);

        vcfWriter = getVCFWriter();
    }

    private VariantContextWriter getVCFWriter() {
        samples = getSamplesForVariants();

        final VCFHeader inputVCFHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);

        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        headerLines.add(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        headerLines.add(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));

        VariantContextWriter writer = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = new IndexedSampleList(samples).asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    /**
     * @param variantContexts  VariantContexts from driving variants with matching start position
     *                         NOTE: This will never be empty
     * @param referenceContext ReferenceContext object covering the reference of the longest spanning VariantContext
     * @param readsContexts
     */
    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        if (currentContig == null) {
            currentContig = variantContexts.get(0).getContig(); //variantContexts should have identical start, so choose 0th arbitrarily
        } else if (!variantContexts.get(0).getContig().equals(currentContig)) {
            processClusters();
            currentContig = variantContexts.get(0).getContig();
        }
        for (final VariantContext vc : variantContexts) {
            final SVCallRecord record = SVCallRecord.createDepthOnlyFromGCNV(vc, minQS);
            if (record != null) {
                defragmenter.add(new SVCallRecordWithEvidence(record));
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        processClusters();
        return null;
    }

    private void processClusters() {
        final List<SVCallRecordWithEvidence> defragmentedCalls = defragmenter.getOutput();
        defragmentedCalls.stream().forEachOrdered(clusterEngine::add);
        //Jack and Isaac cluster first and then defragment
        final List<SVCallRecordWithEvidence> clusteredCalls = clusterEngine.getOutput();
        write(clusteredCalls);
    }

    private void write(final List<SVCallRecordWithEvidence> calls) {
        final List<VariantContext> sortedCalls = calls.stream()
                .sorted(Comparator.comparing(c -> c.getStartAsInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(this::buildVariantContext)
                .collect(Collectors.toList());
        Iterator<VariantContext> it = sortedCalls.iterator();
        ArrayList<VariantContext> overlappingVCs = new ArrayList<>();
        if (!it.hasNext()) {
            return;
        }
        VariantContext prev = it.next();
        overlappingVCs.add(prev);
        int clusterEnd = prev.getEnd();
        String clusterContig = prev.getContig();
        //gather groups of overlapping VCs and update the genotype copy numbers appropriately
        while (it.hasNext()) {
            final VariantContext curr = it.next();
            if (curr.getStart() < clusterEnd && curr.getContig().equals(clusterContig)) {
                overlappingVCs.add(curr);
                if (curr.getEnd() > clusterEnd) {
                    clusterEnd = curr.getEnd();
                }
            } else {
                final List<VariantContext> resolvedVCs = resolveVariantContexts(overlappingVCs);
                for (final VariantContext vc : resolvedVCs) { vcfWriter.add(vc); }
                overlappingVCs = new ArrayList<>();
                overlappingVCs.add(curr);
                clusterEnd = curr.getEnd();
                clusterContig = curr.getContig();
            }
        }
        //write out the last set of overlapping VCs
        final List<VariantContext> resolvedVCs = resolveVariantContexts(overlappingVCs);
        for (final VariantContext vc : resolvedVCs) { vcfWriter.add(vc); }
    }

    /**
     * Ensure genotype calls are consistent for overlapping variant contexts
     * Note that we assume that a sample will not occur twice with the same copy number because it should have been defragmented
     * @param overlappingVCs
     * @return
     */
    private List<VariantContext> resolveVariantContexts(List<VariantContext> overlappingVCs) {
        Utils.nonNull(overlappingVCs);
        final List<VariantContext> resolvedVCs = new ArrayList<>();
        final Iterator<VariantContext> it = overlappingVCs.iterator();

        final Map<String, MutablePair<Integer, Integer>> sampleCopyNumbers = new LinkedHashMap<>();  //sampleName, copyNumber, endPos -- it's safe to just use position because if the VCs overlap then they must be on the same contig
        while (it.hasNext()) {
            final VariantContext curr = it.next();
            //TODO: can we iterate over genotypes for consistency?
            for (final String s : curr.getSampleNames()) {
                //if this sample is in the table and we have a new variant for this sample, update the table
                //TODO: is it better to remove the sample?
                if (sampleCopyNumbers.containsKey(s) && sampleCopyNumbers.get(s).getLeft() != Integer.parseInt(curr.getGenotype(s).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString())) {
                    sampleCopyNumbers.get(s).setRight(Integer.parseInt(curr.getGenotype(s).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()));
                    sampleCopyNumbers.get(s).setLeft(curr.getEnd());
                }
            }
            resolvedVCs.add(updateGenotypes(curr, sampleCopyNumbers));
            //update copy number table for subsequent VCs using variants genotypes from input VCs
            for (final Genotype g : curr.getGenotypes()) {
                if (g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.CN)) {
                    sampleCopyNumbers.put(g.getSampleName(), new MutablePair<>(Integer.parseInt(g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN).toString()), curr.getAttributeAsInt(VCFConstants.END_KEY, curr.getStart())));
                }
            }
        }
        return resolvedVCs;
    }

    /**
     *
     * @param vc
     * @param sampleCopyNumbers may be modified to remove terminated variants
     * @return
     */
    private VariantContext updateGenotypes(final VariantContext vc, final Map<String, MutablePair<Integer, Integer>> sampleCopyNumbers) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final List<Genotype> newGenotypes = new ArrayList<>();
        for (final String sample : samples) {
            if (!sampleCopyNumbers.containsKey(sample) && !vc.hasGenotype(sample)) {
                final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
                genotypeBuilder.alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N));
                genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, HomoSapiensConstants.DEFAULT_PLOIDY);  //NOTE: ploidy 2 assumption
                newGenotypes.add(genotypeBuilder.make());
            } else {
                final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
                if (vc.hasGenotype(sample)) {
                    final int copyNumber = Integer.parseInt(vc.getGenotype(sample).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, HomoSapiensConstants.DEFAULT_PLOIDY).toString());
                    //only for CN=1 can we be sure there's no ref allele
                    if (copyNumber == 0) {
                        genotypeBuilder.alleles(Arrays.asList(GermlineCNVVariantComposer.DEL_ALLELE, GermlineCNVVariantComposer.DEL_ALLELE));
                    } else {
                        genotypeBuilder.alleles(vc.getGenotype(sample).getAlleles()); //will get an extra ref when converted to diploid
                    }
                    genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, copyNumber);
                } else {
                    genotypeBuilder.alleles(Lists.newArrayList(Allele.REF_N, Allele.REF_N));
                }
                if (sampleCopyNumbers.containsKey(sample)) {
                    if (sampleCopyNumbers.get(sample).getRight() > vc.getStart()) {
                        genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, sampleCopyNumbers.get(sample).getLeft());
                    } else {
                        if (!vc.hasGenotype(sample)) {
                            genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, HomoSapiensConstants.DEFAULT_PLOIDY);
                        }
                        sampleCopyNumbers.remove(sample);
                    }
                }
                newGenotypes.add(genotypeBuilder.make());
            }
        }
        return builder.genotypes(newGenotypes).make();
    }

    public VariantContext buildVariantContext(final SVCallRecordWithEvidence call) {
        Utils.nonNull(call);
        final Allele altAllele = Allele.create("<" + call.getType().name() + ">", false);
        final Allele refAllele = Allele.REF_N;
        final VariantContextBuilder builder = new VariantContextBuilder("", call.getContig(), call.getStart(), call.getEnd(),
                Lists.newArrayList(refAllele, altAllele));
        builder.attribute(VCFConstants.END_KEY, call.getEnd());
        builder.attribute(GATKSVVCFConstants.SVLEN, call.getLength());
        builder.attribute(VCFConstants.SVTYPE, call.getType());
        final List<Genotype> genotypes = new ArrayList<>();
        //TODO: I don't need this for loop
        for (final Genotype g : call.getGenotypes()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g);
            genotypeBuilder.alleles(Lists.newArrayList(refAllele, altAllele));
             if (g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.CN)) {
                genotypeBuilder.attribute(GermlineCNVSegmentVariantComposer.CN, g.getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN));
            }

            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);
        return builder.make();
    }

    @Override
    public void closeTool(){
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
