package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.*;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Snaps high AF calls to correct HC",
        oneLineSummary = "PostProcessingMitochondiraM2Calls",
        programGroup = VariantManipulationProgramGroup.class
)
public class PostProcessingM2Calls extends VariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    public File outFile = null;

    private VariantContextWriter vcfWriter = null;

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outFile);
        vcfWriter.writeHeader(getHeaderForVariants());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if(variant.isBiallelic()) {
            double AF = Double.parseDouble(variant.getGenotype(0).getAnyAttribute(VCFConstants.ALLELE_FREQUENCY_KEY).toString());
            if (AF < .5) {
                vcfWriter.add(variant);
            } else {
                final VariantContextBuilder builder = new VariantContextBuilder(variant);
                final GenotypeBuilder gb = new GenotypeBuilder(variant.getGenotype(0));
                int[] AD = variant.getGenotype(0).getAD();
                gb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, (double) AD[1]/(AD[0] + AD[1]));
                builder.genotypes(gb.make());
                vcfWriter.add(builder.make());
            }
        } else {
            vcfWriter.add(variant);
        }
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
