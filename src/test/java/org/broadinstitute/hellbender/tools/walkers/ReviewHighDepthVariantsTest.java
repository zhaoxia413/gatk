package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class ReviewHighDepthVariantsTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleLocusWalker() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                        " -R /Volumes/humgen_gsa-hpprojects/dev/mshand/SpecOps/Mitochondria/MitochondriaOnlyFastas/Homo_sapiens_assembly38.mt_only.fasta" +
                        " -I /Volumes/humgen_gsa-hpprojects/dev/mshand/SpecOps/Mitochondria/Filtering/IGV/mix_14/output.bam" +
                        " -V /Volumes/humgen_gsa-hpprojects/dev/mshand/SpecOps/Mitochondria/Filtering/IGV/mix_14/filtered_05.vcf" +
                        " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleLocusWalkerIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleLocusWalker", this);
    }

}
