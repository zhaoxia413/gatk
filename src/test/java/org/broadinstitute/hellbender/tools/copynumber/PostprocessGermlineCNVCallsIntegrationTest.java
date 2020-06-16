package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration test for {@link PostprocessGermlineCNVCalls}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public final class PostprocessGermlineCNVCallsIntegrationTest extends CommandLineProgramTest {
    private static final BiConsumer<VariantContext, VariantContext> CHECK_VC_START = (actual, expected) ->
            Assert.assertEquals(actual.getStart(), expected.getStart());

    private static final BiConsumer<VariantContext, VariantContext> CHECK_VC_END = (actual, expected) -> {
            Assert.assertTrue(actual.hasAttribute(VCFConstants.END_KEY));
            Assert.assertTrue(expected.hasAttribute(VCFConstants.END_KEY));
            Assert.assertEquals(actual.getAttributeAsInt(VCFConstants.END_KEY, -1), expected.getAttributeAsInt(VCFConstants.END_KEY, -1));
        };

    private static final Predicate<VariantContext> HAS_EXPECTED_ATTRIBUTES = vc -> {
        Genotype g = vc.getGenotype(0);
        return g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.CN) &&
        g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.QS) &&
        g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.QA) &&
        g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.QSE) &&
        g.hasAnyAttribute(GermlineCNVSegmentVariantComposer.QSS);
    };

    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/gcnv-postprocess");

    private static final List<String> CALL_SHARDS = Arrays.asList(
            new File(TEST_SUB_DIR, "shard_0-calls").getAbsolutePath(),
            new File(TEST_SUB_DIR, "shard_1-calls").getAbsolutePath(),
            new File(TEST_SUB_DIR, "shard_2-calls").getAbsolutePath());

    private static final List<String> MODEL_SHARDS = Arrays.asList(
            new File(TEST_SUB_DIR, "shard_0-model").getAbsolutePath(),
            new File(TEST_SUB_DIR, "shard_1-model").getAbsolutePath(),
            new File(TEST_SUB_DIR, "shard_2-model").getAbsolutePath());

    private static final String PLOIDY_CALLS = new File(TEST_SUB_DIR, "ploidy-calls").getAbsolutePath();

    private static final int AUTOSOMAL_REF_COPY_NUMBER = 2;

    private final List<String> ALLOSOMAL_CONTIGS = Arrays.asList("X", "Y");

    private static final List<File> INTERVALS_VCF_CORRECT_OUTPUTS = Arrays.asList(
            new File(TEST_SUB_DIR, "intervals_output_SAMPLE_000.vcf.gz"),
            new File(TEST_SUB_DIR, "intervals_output_SAMPLE_001.vcf.gz"),
            new File(TEST_SUB_DIR, "intervals_output_SAMPLE_002.vcf.gz"));

    private static final List<File> SEGMENTS_VCF_CORRECT_OUTPUTS = Arrays.asList(
            new File(TEST_SUB_DIR, "segments_output_SAMPLE_000.vcf"),
            new File(TEST_SUB_DIR, "segments_output_SAMPLE_001.vcf"),
            new File(TEST_SUB_DIR, "segments_output_SAMPLE_002.vcf"));

    private static final List<File> DENOISED_COPY_RATIOS_OUTPUTS = Arrays.asList(
            new File(TEST_SUB_DIR, "denoised_copy_ratios_SAMPLE_000.tsv"),
            new File(TEST_SUB_DIR, "denoised_copy_ratios_SAMPLE_001.tsv"),
            new File(TEST_SUB_DIR, "denoised_copy_ratios_SAMPLE_002.tsv"));

    private static final int NUM_TEST_SAMPLES = 3;

    private static final File COMBINED_INTERVALS_VCF = new File(toolsTestDir + "copynumber/gcnv-postprocess/intervals.combined.vcf.gz");
    private static final File CLUSTERED_VCF = new File(toolsTestDir + "copynumber/clustering/threeSamples.vcf.gz");

    /**
     * Runs {@link PostprocessGermlineCNVCalls} for a single sample. If {@code segmentsOutputVCF} is null,
     * the tool will only generate intervals VCF output (which is the expected behavior).
     */
    private ArgumentsBuilder getArgsForSingleSample(final List<String> callShards,
                                        final List<String> modelShards,
                                        final int sampleIndex,
                                        final File intervalsOutputVCF,
                                        final File segmentsOutputVCF,
                                        final File denoisedCopyRatiosOutput,
                                        final List<String> allosomalContigs,
                                        final int refAutosomalCopyNumber) {
        final ArgumentsBuilder argumentsBuilder = new ArgumentsBuilder();
        argumentsBuilder.add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE,
                "false");

        argumentsBuilder.add(PostprocessGermlineCNVCalls.SAMPLE_INDEX_LONG_NAME,
                String.valueOf(sampleIndex));
        callShards.forEach(callDir -> argumentsBuilder.add(
                PostprocessGermlineCNVCalls.CALLS_SHARD_PATH_LONG_NAME, callDir));
        argumentsBuilder.add(PostprocessGermlineCNVCalls.OUTPUT_INTERVALS_VCF_LONG_NAME,
                intervalsOutputVCF.getAbsolutePath());
        allosomalContigs.forEach(contig -> argumentsBuilder.add(
                PostprocessGermlineCNVCalls.ALLOSOMAL_CONTIG_LONG_NAME, contig));
        argumentsBuilder.add(PostprocessGermlineCNVCalls.AUTOSOMAL_REF_COPY_NUMBER_LONG_NAME,
                String.valueOf(refAutosomalCopyNumber));

        /* add required arguments for segments VCF creation */
        argumentsBuilder.add(PostprocessGermlineCNVCalls.CONTIG_PLOIDY_CALLS_LONG_NAME,
                PLOIDY_CALLS);
        modelShards.forEach(modelDir -> argumentsBuilder.add(
                PostprocessGermlineCNVCalls.MODEL_SHARD_PATH_LONG_NAME, modelDir));
        argumentsBuilder.add(PostprocessGermlineCNVCalls.OUTPUT_SEGMENTS_VCF_LONG_NAME,
                segmentsOutputVCF.getAbsolutePath());

        /* add denoised copy ratio output file path */
        argumentsBuilder.add(PostprocessGermlineCNVCalls.OUTPUT_DENOISED_COPY_RATIOS_LONG_NAME,
                denoisedCopyRatiosOutput.getAbsolutePath());

        return argumentsBuilder;
    }

    private ArgumentsBuilder getArgsWithBreakpoints(final List<String> callShards,
                                        final List<String> modelShards,
                                        final int sampleIndex,
                                        final File intervalsOutputVCF,
                                        final File segmentsOutputVCF,
                                        final File denoisedCopyRatiosOutput,
                                        final List<String> allosomalContigs,
                                        final int refAutosomalCopyNumber,
                                        final File combinedIntervalsVCF,
                                        final File clusteredVCF) {
        ArgumentsBuilder args = getArgsForSingleSample(callShards, modelShards, sampleIndex, intervalsOutputVCF, segmentsOutputVCF, denoisedCopyRatiosOutput, allosomalContigs, refAutosomalCopyNumber);
        args.add(PostprocessGermlineCNVCalls.CLUSTERED_FILE_LONG_NAME, clusteredVCF);
        args.add(PostprocessGermlineCNVCalls.COMBINED_INTERVALS_LONG_NAME, combinedIntervalsVCF);
        return args;
    }

    @Test(dataProvider = "differentValidInput", groups = {"python"})
    public void testDifferentValidInput(final int sampleIndex,
                                        final List<String> callShards,
                                        final List<String> modelShards) throws IOException {
        final File actualIntervalsOutputVCF = createTempFile("intervals-output-vcf-" + sampleIndex, ".vcf");
        final File actualSegmentsOutputVCF = createTempFile("segments-output-vcf-" + sampleIndex, ".vcf");
        final File intervalsIndex = new File(actualIntervalsOutputVCF.getAbsolutePath() + FileExtensions.TRIBBLE_INDEX);
        Assert.assertFalse(intervalsIndex.exists());
        final File segmentsIndex = new File(actualIntervalsOutputVCF.getAbsolutePath() + FileExtensions.TRIBBLE_INDEX);
        Assert.assertFalse(segmentsIndex.exists());
        final File actualDenoisedCopyRatiosOutput = createTempFile("denoised-copy-ratios-output-" + sampleIndex, ".tsv");
        final File expectedIntervalsOutputVCF = INTERVALS_VCF_CORRECT_OUTPUTS.get(sampleIndex);
        final File expectedSegmentsOutputVCF = SEGMENTS_VCF_CORRECT_OUTPUTS.get(sampleIndex);
        final File expectedDenoisedCopyRatiosOutput = DENOISED_COPY_RATIOS_OUTPUTS.get(sampleIndex);
        final ArgumentsBuilder args = getArgsForSingleSample(callShards, modelShards, sampleIndex,
                actualIntervalsOutputVCF, actualSegmentsOutputVCF, actualDenoisedCopyRatiosOutput,
                ALLOSOMAL_CONTIGS, AUTOSOMAL_REF_COPY_NUMBER);
        runCommandLine(args);

        Assert.assertTrue(intervalsIndex.exists());
        Assert.assertTrue(segmentsIndex.exists());
        IntegrationTestSpec.assertEqualTextFiles(actualIntervalsOutputVCF, expectedIntervalsOutputVCF, "##");
        final VCFHeader intervalsHeader = VariantContextTestUtils.getVCFHeader(actualIntervalsOutputVCF.getAbsolutePath());
        Assert.assertTrue(intervalsHeader.getContigLines().size() > 0);
        IntegrationTestSpec.assertEqualTextFiles(actualSegmentsOutputVCF, expectedSegmentsOutputVCF, "##");
        final VCFHeader segmentsHeader = VariantContextTestUtils.getVCFHeader(actualIntervalsOutputVCF.getAbsolutePath());
        Assert.assertTrue(segmentsHeader.getContigLines().size() > 0);
        IntegrationTestSpec.assertEqualTextFiles(actualDenoisedCopyRatiosOutput, expectedDenoisedCopyRatiosOutput, "##");
    }

    @Test(dataProvider = "differentInvalidInput", expectedExceptions = IllegalArgumentException.class, groups = {"python"})
    public void testDifferentInvalidInput(final int sampleIndex,
                                          final List<String> callShards,
                                          final List<String> modelShards) throws IOException {
        final ArgumentsBuilder args = getArgsForSingleSample(callShards, modelShards, sampleIndex,
                createTempFile("intervals-output-vcf", ".vcf"),
                createTempFile("segments-output-vcf", ".vcf"),
                createTempFile("denoised-copy-ratios-output", ".tsv"),
                ALLOSOMAL_CONTIGS, AUTOSOMAL_REF_COPY_NUMBER);
        runCommandLine(args);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, groups = {"python"})
    public void testBadAutosomalContigs() {
        final ArgumentsBuilder args = getArgsForSingleSample(CALL_SHARDS, MODEL_SHARDS, 0,
                createTempFile("intervals-output-vcf", ".vcf"),
                createTempFile("segments-output-vcf", ".vcf"),
                createTempFile("denoised-copy-ratios-output", ".tsv"),
                Collections.singletonList("Z"), /* unknown contig */
                AUTOSOMAL_REF_COPY_NUMBER);
        runCommandLine(args);
    }

    @Test(groups = {"python"})
    public void testQualScoreCalculationWithBreakpoints() {
        //run one all-reference sample
        final File segmentsOutput = createTempFile("segments-output-vcf-0", ".vcf");
        final ArgumentsBuilder args = getArgsWithBreakpoints(CALL_SHARDS, MODEL_SHARDS, 0,
                createTempFile("intervals-output-vcf", ".vcf"),
                segmentsOutput,
                createTempFile("denoised-copy-ratios-output", ".tsv"),
                ALLOSOMAL_CONTIGS, 2, new File(TEST_SUB_DIR, "intervals_output_SAMPLE_000.vcf.gz"), CLUSTERED_VCF);
        runCommandLine(args);

        //...and one sample with variants
        final File segmentsOutput2 = createTempFile("segments-output-vcf", ".vcf");
        final ArgumentsBuilder args2 = getArgsWithBreakpoints(CALL_SHARDS, MODEL_SHARDS, 1,
                createTempFile("intervals-output-vcf", ".vcf"),
                segmentsOutput2,
                createTempFile("denoised-copy-ratios-output", ".tsv"),
                ALLOSOMAL_CONTIGS, 2, new File(TEST_SUB_DIR, "intervals_output_SAMPLE_001.vcf.gz"), CLUSTERED_VCF);
        runCommandLine(args2);

        final Pair<VCFHeader, List<VariantContext>> output = VariantContextTestUtils.readEntireVCFIntoMemory(segmentsOutput.getAbsolutePath());
        final Pair<VCFHeader, List<VariantContext>> output2 = VariantContextTestUtils.readEntireVCFIntoMemory(segmentsOutput2.getAbsolutePath());
        final Pair<VCFHeader, List<VariantContext>> sample001Segments = VariantContextTestUtils.readEntireVCFIntoMemory(SEGMENTS_VCF_CORRECT_OUTPUTS.get(1).getAbsolutePath());
        final Pair<VCFHeader, List<VariantContext>> clusteredBreakpoints = VariantContextTestUtils.readEntireVCFIntoMemory(CLUSTERED_VCF.getAbsolutePath());
        Assert.assertEquals(output.getRight().size(), clusteredBreakpoints.getRight().size());
        Assert.assertEquals(output2.getRight().size(), clusteredBreakpoints.getRight().size());
        Assert.assertTrue(output.getRight().stream().allMatch(HAS_EXPECTED_ATTRIBUTES));
        Assert.assertTrue(output2.getRight().stream().allMatch(HAS_EXPECTED_ATTRIBUTES));

        //clustered VCF starts with chromosome 2 because all samples are ref over chr 1
        Assert.assertEquals(clusteredBreakpoints.getRight().get(0).getContig(), output.getRight().get(0).getContig());

        //sample0 is all reference
        Assert.assertTrue(output.getRight().stream().allMatch(vc -> vc.getGenotypes().get(0).isHomRef()));

        //all VCs should have matching start and end
        GATKBaseTest.assertCondition(clusteredBreakpoints.getRight(), output.getRight(), CHECK_VC_START);
        GATKBaseTest.assertCondition(clusteredBreakpoints.getRight(), output.getRight(), CHECK_VC_END);
        GATKBaseTest.assertCondition(clusteredBreakpoints.getRight(), output2.getRight(), CHECK_VC_START);
        GATKBaseTest.assertCondition(clusteredBreakpoints.getRight(), output2.getRight(), CHECK_VC_END);

        //sample001 has a del and a dup on chromosome 2
        Assert.assertTrue(output2.getRight().stream().anyMatch(vc -> vc.getContig().equals("2") &&
                vc.getAlternateAllele(0).equals(GermlineCNVSegmentVariantComposer.DEL_ALLELE) && vc.getGenotypes().get(0).isHomVar() && //calls are "haploid"
                vc.getGenotype(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString().equals("0")));
        Assert.assertTrue(output2.getRight().stream().anyMatch(vc -> vc.getContig().equals("2") &&
                vc.getAlternateAllele(0).equals(GermlineCNVSegmentVariantComposer.DUP_ALLELE) && vc.getGenotypes().get(0).isHomVar() && //calls are "haploid"
                vc.getGenotype(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString().equals("4")));

        //compare sample001 with new QUALs against its single-sample segmented results
        // should have same Q score where breakpoints match, zero if breakpoint was moved in clustering
        for (final VariantContext segment : sample001Segments.getRight()) {
            if (segment.isVariant()) {
                final List<VariantContext> matches = output2.getRight().stream().filter(vc -> vc.getStart() == segment.getStart()).collect(Collectors.toList());
                if (matches.size() == 1) {
                    Assert.assertEquals(segment.getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN), matches.get(0).getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.CN));
                    Assert.assertEquals(segment.getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.QSS), matches.get(0).getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.QSS));
                    GATKBaseTest.assertEqualsIntSmart(Integer.parseInt(segment.getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.QS).toString()),
                            Integer.parseInt(matches.get(0).getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.QS).toString()), 2,  //I've seen some wiggle in a few variants
                            "QS score varies by more than 2.");
                    if (matches.get(0).getAttribute(VCFConstants.END_KEY).equals(segment.getAttribute(VCFConstants.END_KEY))) {
                        Assert.assertEquals(segment.getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.QSE), matches.get(0).getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.QSE));
                        Assert.assertEquals(segment.getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.NP), matches.get(0).getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.NP));
                    } else {
                        Assert.assertEquals(Integer.parseInt(matches.get(0).getGenotype(0).getExtendedAttribute(GermlineCNVSegmentVariantComposer.QSE).toString()), 0);
                    }
                }
            }
        }
    }

    @DataProvider(name = "differentValidInput")
    public Object[][] getDifferentValidInputTestData() {
        final List<Object[]> testCases = new ArrayList<>();
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* correct order */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* permutation */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(2), CALL_SHARDS.get(1), CALL_SHARDS.get(0)),
                        Arrays.asList(MODEL_SHARDS.get(2), MODEL_SHARDS.get(0), MODEL_SHARDS.get(1))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* permutation */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(1), CALL_SHARDS.get(2), CALL_SHARDS.get(0)),
                        Arrays.asList(MODEL_SHARDS.get(1), MODEL_SHARDS.get(0), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        return testCases.toArray(new Object[testCases.size()][]);
    }

    @DataProvider(name = "differentInvalidInput")
    public Object[][] getDifferentInvalidInput() {
        final List<Object[]> testCases = new ArrayList<>();
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* fewer model shards than call shards */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* fewer call shards than model shards */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[]{
                        /* repeated calls */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(2), CALL_SHARDS.get(2), CALL_SHARDS.get(1)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[]{
                        /* repeated models */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(1), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        return testCases.toArray(new Object[testCases.size()][]);
    }
}
