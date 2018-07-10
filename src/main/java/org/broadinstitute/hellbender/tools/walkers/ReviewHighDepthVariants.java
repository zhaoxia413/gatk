package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.util.List;


public class ReviewHighDepthVariants extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public String output;
    private SAMFileGATKReadWriter outputWriter;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files", optional = true)
    private List<FeatureInput<VariantContext>> variants;


    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(output), true);
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if(featureContext.hasBackingDataSource()) {
            List<VariantContext> vars = featureContext.getValues(variants);
            if(!vars.isEmpty()) {
                for (VariantContext variant : vars) {
                    if (variant.isBiallelic() && variant.isSNP()) {
                        ReadPileup pileup = alignmentContext.getBasePileup();
                        List<Integer> offsets = pileup.getOffsets();
                        List<GATKRead> reads = pileup.getReads();
                        for(int i=0; i<reads.size(); i++){
                            GATKRead read = reads.get(i);
                            byte base = read.getBase(offsets.get(i));
                            if( base != referenceContext.getBase() ) {
                                outputWriter.addRead(read);
                            }
                        }
                    }
                }
            }
        }
    }
}
