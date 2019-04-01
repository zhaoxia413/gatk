package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;

public class CollectPolymeraseSlippageDataIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test() throws Exception {
        final File output = createTempFile("output", ".txt");
        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37Reference,
                "-L", "20",
                "-O", output.getAbsolutePath()
        };

        runCommandLine(args);
        int g = 5;
    }
}