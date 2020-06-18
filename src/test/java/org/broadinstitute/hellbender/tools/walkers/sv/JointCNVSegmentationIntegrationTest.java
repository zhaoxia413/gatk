package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.GATKBaseTest.toolsTestDir;
import static org.broadinstitute.hellbender.testutils.BaseTest.createTempFile;
import static org.testng.Assert.*;

public class JointCNVSegmentationIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/gcnv-postprocess");

    private static final List<File> SEGMENTS_VCF_CORRECT_OUTPUTS = Arrays.asList(
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_000.vcf"),
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_001.vcf"),
    new File(TEST_SUB_DIR, "segments_output_SAMPLE_002.vcf"));

    @DataProvider
    public Object[][] postprocessOutputs() {
        return new Object[][] {
               new Object[]{SEGMENTS_VCF_CORRECT_OUTPUTS}
        };
    }

    @DataProvider
    public Object[][] overlappingSamples() {
        return new Object[][] {
            new Object[] {
                    Arrays.asList(new File(getToolTestDataDir() + "HG00365.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "HG01623.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "HG01789.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "HG02165.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "HG02221.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA07357.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA11829.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA12005.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA12046.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA12814.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA12873.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA18946.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA18997.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA19428.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA19456.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA20502.overlaps.vcf.gz"),
                    new File(getToolTestDataDir() + "NA21120.overlaps.vcf.gz"))
            }
        };
    }

    @Test(dataProvider = "postprocessOutputs")
    public void testThreeGCNVSamples(List<File> inputVcfs) {
        final File output = createTempFile("threeSamples", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(GATKBaseTest.b37Reference)
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "threeSamples.interval_list")
                .addIntervals(new File(getToolTestDataDir() + "threeSamples.interval_list"));

        inputVcfs.forEach(vcf -> args.addVCF(vcf));

        runCommandLine(args, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> withQStreshold = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(withQStreshold.getRight().size(), 7);

        final File output2 = createTempFile("threeSamples.noQSthreshold",".vcf");

        final ArgumentsBuilder args2 = new ArgumentsBuilder()
                .addOutput(output2)
                .addReference(GATKBaseTest.b37Reference)
                .add(JointCNVSegmentation.MIN_QUALITY_LONG_NAME, 0)
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "threeSamples.interval_list")
                .addIntervals(new File(getToolTestDataDir() + "threeSamples.interval_list"));
        inputVcfs.forEach(vcf -> args2.addVCF(vcf));

        runCommandLine(args2, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> withoutQStreshold = VariantContextTestUtils.readEntireVCFIntoMemory(output2.getAbsolutePath());
        Assert.assertEquals(withoutQStreshold.getRight().size(), 8);
        //extra variant at X:227988

        //another test to make sure adjacent events with different copy numbers don't get merged/defragmented?
    }

    @Test
    public void testDefragmentation() {
        final File output = createTempFile("defragmented",".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(GATKBaseTest.b37Reference)
                .addVCF(getToolTestDataDir() + "NA20533.fragmented.segments.vcf.gz")
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "intervals.chr13.interval_list")
                .addInterval("13:52951204-115064572");

        runCommandLine(args, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> defragmentedEvents = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(defragmentedEvents.getRight().size(), 1);
        Assert.assertEquals(defragmentedEvents.getRight().get(0).getAttributeAsInt(GATKSVVCFConstants.SVLEN,0), 62113369);
    }

    @Test(dataProvider = "overlappingSamples")
    public void testOverlappingEvents(final List<File> inputVcfs) {
        final File output = createTempFile("overlaps", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .addReference(GATKBaseTest.b37Reference)
                .add(JointCNVSegmentation.MODEL_CALL_INTERVALS, getToolTestDataDir() + "intervals.chr22.interval_list")
                .addInterval("22:22,538,114-23,538,437");

        inputVcfs.forEach(vcf -> args.addVCF(vcf));

        runCommandLine(args, JointCNVSegmentation.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> overlappingEvents = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(overlappingEvents.getRight().size(), 6);
        //do copy number checks on genotypes
        //at the start of the contig, all homRef genotypes should be CN2
        final VariantContext vc0 = overlappingEvents.getRight().get(0);
        for (final Genotype g : vc0.getGenotypes()) {
            if (g.isHomRef()) {
                Assert.assertEquals(Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 2);
            } else {
                Assert.assertNotEquals(Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()), 2);
            }
        }

        //TODO: more procedural checks
    }
}