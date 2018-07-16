package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;


public class DiscordantPairArtifact extends GenotypeAnnotation implements StandardMutectAnnotation{

    public static final String NON_MT_OA = "NON_MT_OA";
    public static final String DISCORDANT_PAIRS = "DISC_PAIRS";

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        Utils.nonNull(likelihoods);

        // do not annotate the genotype of the normal sample
        if (g.isHomRef()){
            return;
        }

        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY, () -> null, -1);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
        final Allele altAlelle = vc.getAlternateAllele(indexOfMaxTumorLod);
        final Allele refAllele = vc.getReference();

        Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(g.getSampleName());

        int discordantAlt = (int) bestAlleles.stream().filter(ba -> ba.read.hasAttribute("OA") && ba.isInformative()
                && ba.allele.equals(altAlelle)
                && !ba.read.getAttributeAsString("OA").split(",")[0].equals(ba.read.getAttributeAsString("XO").split(",")[0])).count();

        int nonChrMAlt = (int) bestAlleles.stream().filter(ba -> ba.read.hasAttribute("OA") && ba.isInformative() &&
                ba.allele.equals(altAlelle) && !ba.read.getAttributeAsString("OA").split(",")[0].equals("chrM")).count();
        int nonChrMRef = (int) bestAlleles.stream().filter(ba -> ba.read.hasAttribute("OA") && ba.isInformative() &&
                ba.allele.equals(refAllele) && !ba.read.getAttributeAsString("OA").split(",")[0].equals("chrM")).count();

        final int[] nonChrMCounts = new int[2];
        nonChrMCounts[0] = nonChrMRef;
        nonChrMCounts[1] = nonChrMAlt;

        gb.attribute(NON_MT_OA, nonChrMCounts);
        gb.attribute(DISCORDANT_PAIRS, discordantAlt);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(NON_MT_OA, 2, VCFHeaderLineType.Integer, "number of non MT ref reads and number of non MT alt reads"),
            new VCFFormatHeaderLine(DISCORDANT_PAIRS, 1, VCFHeaderLineType.Integer, "number of discordant pairs in alt"));
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(NON_MT_OA);
    }
}
