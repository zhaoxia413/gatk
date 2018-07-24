package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.stream.StreamSupport;

public class SplitMultiAllelicSitesTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @Test
    public void testExampleVariantWalker() throws IOException {

        final File outputVcf = createTempFile("output", ".vcf");

        final String[] args = {
                "-V", "/Volumes/humgen_gsa-hpprojects/dev/mshand/SpecOps/Mitochondria/Filtering/IGV/09C96685_vs_163516/filtered_2.vcf",
                "-O", outputVcf.getAbsolutePath(),
        };

        runCommandLine(args);
        final long numVariants = StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false).count();
        Assert.assertTrue(numVariants < 4);

    }
}
