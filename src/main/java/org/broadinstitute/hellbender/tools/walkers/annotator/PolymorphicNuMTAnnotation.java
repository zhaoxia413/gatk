package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class PolymorphicNuMTAnnotation extends GenotypeAnnotation implements Annotation {

    public static final String POTENTIAL_POLYMORPHIC_NUMT = "POTENTIAL_POLYMORPHIC_NUMT";

    private double LOWER_BOUND_PROB = .05;

    private double LAMBDA;
    private PoissonDistribution AUTOSOMAL_COVERAGE = new PoissonDistribution(LAMBDA);
    private int MIN_AUTOSOMAL_HOM_ALT = AUTOSOMAL_COVERAGE.inverseCumulativeProbability(LOWER_BOUND_PROB);
    private int MAX_AUTOSOMAL_HOM_ALT = AUTOSOMAL_COVERAGE.inverseCumulativeProbability(1 - LOWER_BOUND_PROB);
    private int MIN_AUTOSOMAL_HET = MIN_AUTOSOMAL_HOM_ALT / 2;
    private int MAX_AUTOSOMAL_HET = MAX_AUTOSOMAL_HOM_ALT / 2;

    public PolymorphicNuMTAnnotation(final double lambda){
        this.LAMBDA = lambda;
    }

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        Utils.nonNull(likelihoods);

        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY, () -> null, -1);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
        final Allele altAlelle = vc.getAlternateAllele(indexOfMaxTumorLod);

        Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(g.getSampleName());
        final int numAltReads = (int) bestAlleles.stream().filter(ba -> ba.isInformative() && ba.allele.equals(altAlelle)).count();

        if ((numAltReads > MIN_AUTOSOMAL_HOM_ALT && numAltReads < MAX_AUTOSOMAL_HOM_ALT) || (numAltReads > MIN_AUTOSOMAL_HET && numAltReads < MAX_AUTOSOMAL_HET)) {
            gb.attribute(POTENTIAL_POLYMORPHIC_NUMT, true);
        }
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(POTENTIAL_POLYMORPHIC_NUMT, 1, VCFHeaderLineType.Flag, "Potentially a polymorphic NuMT false positive rather than a real mitochondrial variant."));
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(POTENTIAL_POLYMORPHIC_NUMT);
    }
}
