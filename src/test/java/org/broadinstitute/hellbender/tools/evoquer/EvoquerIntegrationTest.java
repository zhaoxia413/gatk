package org.broadinstitute.hellbender.tools.evoquer;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.CombineGVCFs;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.vcf.SortVcf;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

public class EvoquerIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testChr20Dalio3ExomesWithGnarlyGenotyper() throws IOException {
        final String projectID = "broad-dsp-spec-ops";
        final String datasetMapString = "chr20   joint_genotyping_chr20_integration_test  pet_without_gq60_ir_c_sam_st vet";
        final String interval = "chr20";
        final File outputVCF = createTempFile("output", ".vcf");

        final File combineGVCFOutput = new File("src/test/resources/large/integration_test_3_sample.vcf");

        final File datasetMapFile = createTempFile("testChr203Exomes", ".dataset_map");
        try ( final PrintWriter writer = new PrintWriter(datasetMapFile) ) {
            writer.println(datasetMapString);
        }

        final String[] args = {
                "--project-id", projectID,
                "--dataset-map", datasetMapFile.getAbsolutePath(),
                "-R", hg38Reference,
                "-L", interval,
                "-O", outputVCF.getAbsolutePath(),
                "--run-query-only", "false",
                "--disable-gnarly-genotyper", "false"
        };

        runCommandLine(args);

        final List<VariantContext> variants = VariantContextTestUtils.streamVcf(outputVCF)
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        Assert.assertTrue(variants.size() > 0, "Output is empty.");

        final List<VariantContext> expectedVariants = VariantContextTestUtils.streamVcf(combineGVCFOutput)
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        //Assert.assertEquals(variants.size(), expectedVariants.size(), "Number of variants doesn't match expected.");
        final Iterator<VariantContext> expectedVariantsIterator = expectedVariants.iterator();

        for (VariantContext v : variants) {
            VariantContext expectedV = expectedVariantsIterator.next();
            Assert.assertEquals(v.getStart(), expectedV.getStart(), "Positions don't match.");
            //Assert.assertEquals(v.getAttributes(), expectedV.getAttributes());
            /*Assert.assertEquals(v.getAttribute("AC"), expectedV.getAttribute("AC"), "Position: " + v.getStart() + " does't have matching AC.");
            Assert.assertEquals(v.getAttribute("AF"), expectedV.getAttribute("AF"), "Position: " + v.getStart() + " does't have matching AF.");
            Assert.assertEquals(v.getAttribute("AN"), expectedV.getAttribute("AN"), "Position: " + v.getStart() + " does't have matching AN.");
            Assert.assertEquals(v.getAttribute("AS_AltDP"), expectedV.getAttribute("AS_AltDP"), "Position: " + v.getStart() + " does't have matching AS_AltDP.");
            Assert.assertEquals(v.getAttribute("AS_BaseQRankSum"), expectedV.getAttribute("AS_BaseQRankSum"), "Position: " + v.getStart() + " does't have matching AS_BaseQRankSum.");
            Assert.assertEquals(v.getAttribute("AS_FS"), expectedV.getAttribute("AS_FS"), "Position: " + v.getStart() + " does't have matching AS_FS.");
            Assert.assertEquals(v.getAttribute("AS_MQ"), expectedV.getAttribute("AS_MQ"), "Position: " + v.getStart() + " does't have matching AS_MQ.");
            Assert.assertEquals(v.getAttribute("AS_MQRankSum"), expectedV.getAttribute("AS_MQRankSum"), "Position: " + v.getStart() + " does't have matching AS_MQRankSum.");
            Assert.assertEquals(v.getAttribute("AS_QD"), expectedV.getAttribute("AS_QD"), "Position: " + v.getStart() + " does't have matching AS_QD.");
            Assert.assertEquals(v.getAttribute("AS_ReadPosRankSum"), expectedV.getAttribute("AS_ReadPosRankSum"), "Position: " + v.getStart() + " does't have matching AS_ReadPosRankSum.");
            Assert.assertEquals(v.getAttribute("AS_SOR"), expectedV.getAttribute("AS_SOR"), "Position: " + v.getStart() + " does't have matching AS_SOR.");
            Assert.assertEquals(v.getAttribute("BaseQRankSum"), expectedV.getAttribute("BaseQRankSum"), "Position: " + v.getStart() + " does't have matching BaseQRankSum.");
            Assert.assertEquals(v.getAttributeAsInt("DP", 0), expectedV.getAttributeAsInt("DP", 0), "Position: " + v.getStart() + " does't have matching DP.");
            Assert.assertEquals(v.getAttribute("ExcessHet"), expectedV.getAttribute("ExcessHet"), "Position: " + v.getStart() + " does't have matching ExcessHet.");
            Assert.assertEquals(v.getAttributeAsInt("FS", 0), expectedV.getAttributeAsInt("FS", 0), "Position: " + v.getStart() + " does't have matching FS.");
            Assert.assertEquals(v.getAttributeAsInt("MQ", 0), expectedV.getAttributeAsInt("MQ", 0), "Position: " + v.getStart() + " does't have matching MQ.");*/
            for (String s : v.getSampleNames()) {
                Genotype g = v.getGenotype(s);
                Genotype expectedG = expectedV.getGenotype(s);
                Assert.assertEquals(g.getGQ(), expectedG.getGQ(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching GQ.");
                if (!g.isHomRef()) {
                    Assert.assertEquals(g.getAnyAttribute("GT"), expectedG.getAnyAttribute("GT"), "Position: " + v.getStart() + " Sample: " + s + " does't have matching GT.");
                    Assert.assertEquals(g.getAD(), expectedG.getAD(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching AD.");
                    Assert.assertEquals(g.getDP(), expectedG.getDP(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching DP.");
                    Assert.assertEquals(g.getPL(), expectedG.getPL(), "Position: " + v.getStart() + " Sample: " + s + " does't have matching PL.");
                }
            }
        }
    }
}
