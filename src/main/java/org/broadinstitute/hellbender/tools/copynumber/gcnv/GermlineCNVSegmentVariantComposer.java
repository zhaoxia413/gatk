package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.hadoop.yarn.webapp.hamlet.Hamlet;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Helper class for {@link PostprocessGermlineCNVCalls} for single-sample postprocessing of segmented
 * {@link GermlineCNVCaller} calls.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVSegmentVariantComposer extends GermlineCNVVariantComposer<IntegerCopyNumberSegment> {

    /* VCF FORMAT header keys */

    /**
     * Segment copy-number call
     */
    public static final String CN = "CN";

    /**
     * Number of points in the segment
     */
    public static final String NP = "NP";

    /**
     * Quality metric (some points called)
     */
    public static final String QS = "QS";

    /**
     * Quality metric (all points called)
     */
    public static final String QA = "QA";

    /**
     * Quality metric (segment start)
     */
    public static final String QSS = "QSS";

    /**
     * Quality metric (segment end)
     */
    public static final String QSE = "QSE";

    private final IntegerCopyNumberState refAutosomalCopyNumberState;
    private final Set<String> allosomalContigSet;
    private final ReferenceSequenceFile reference;

    /**
     * Constructor.
     *
     * @param outputWriter variant context writer
     * @param sampleName sample name
     * @param refAutosomalCopyNumberState ref copy-number state on autosomal contigs
     * @param allosomalContigSet set of allosomal contigs (ref copy-number allele be chosen according to
     *                           given contig baseline copy-number states)
     * @param reference may be null
     */
    public GermlineCNVSegmentVariantComposer(final VariantContextWriter outputWriter,
                                             final String sampleName,
                                             final IntegerCopyNumberState refAutosomalCopyNumberState,
                                             final Set<String> allosomalContigSet,
                                             final ReferenceSequenceFile reference) {
        super(outputWriter, sampleName);
        this.refAutosomalCopyNumberState = Utils.nonNull(refAutosomalCopyNumberState);
        this.allosomalContigSet = Utils.nonNull(allosomalContigSet);
        this.reference = reference;
    }

    @Override
    public void composeVariantContextHeader(final SAMSequenceDictionary sequenceDictionary,
                                            final Set<VCFHeaderLine> vcfDefaultToolHeaderLines) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Collections.singletonList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add default tool header lines */
        vcfDefaultToolHeaderLines.forEach(result::addMetaDataLine);

        result.setSequenceDictionary(sequenceDictionary);

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.Integer, "Segment genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN, 1,
                VCFHeaderLineType.Integer, "Segment most-likely copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(NP, 1,
                VCFHeaderLineType.Integer, "Number of points (i.e. targets or bins) in the segment"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QS, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that at least one point " +
                "(i.e. target or bin) in the segment agrees with the segment copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QA, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that all points " +
                "(i.e. targets or bins) in the segment agree with the segment copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QSS, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that the segment start " +
                "position is a genuine copy-number changepoint"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QSE, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that the segment end " +
                "position is a genuine copy-number changepoint"));

        /* INFO header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of the variant"));
        outputWriter.writeHeader(result);
    }

    /**
     * Compose a variant context from a given {@link IntegerCopyNumberSegment}
     *
     * @param segment an instance of {@link IntegerCopyNumberSegment}
     * @return composed variant context
     */
    @VisibleForTesting
    VariantContext composeVariantContext(final IntegerCopyNumberSegment segment) {
        final String contig = segment.getContig();
        final int start = segment.getStart();
        final int end = segment.getEnd();
        final int copyNumberCall = segment.getCallIntegerCopyNumberState().getCopyNumber();
        final Allele refAllele = reference == null ? REF_ALLELE : Allele.create(ReferenceUtils.getRefBaseAtPosition(reference, contig, start), true);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.chr(contig);
        variantContextBuilder.start(start);
        variantContextBuilder.stop(end);
        variantContextBuilder.id(String.format(VARIANT_PREFIX + "_%s_%d_%d", contig, start, end));

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final IntegerCopyNumberState refCopyNumber = allosomalContigSet.contains(contig)
                ? segment.getBaselineIntegerCopyNumberState()
                : refAutosomalCopyNumberState;
        genotypeBuilder.alleles(makeGenotypeAlleles(copyNumberCall, refCopyNumber.getCopyNumber(), refAllele));
        genotypeBuilder.attribute(CN, copyNumberCall);
        genotypeBuilder.attribute(NP, segment.getNumPoints());
        genotypeBuilder.attribute(QS, FastMath.round(segment.getQualitySomeCalled()));
        genotypeBuilder.attribute(QA, FastMath.round(segment.getQualityAllCalled()));
        genotypeBuilder.attribute(QSS, FastMath.round(segment.getQualityStart()));
        genotypeBuilder.attribute(QSE, FastMath.round(segment.getQualityEnd()));
        final Genotype genotype = genotypeBuilder.make();

        final Set<Allele> uniquifiedAlleles = new HashSet<>();
        uniquifiedAlleles.add(refAllele);
        if (copyNumberCall > refCopyNumber.getCopyNumber()) {
            uniquifiedAlleles.add(DUP_ALLELE);  //dupes need additional alts since their genotypes are no-call
        } else if (copyNumberCall < refCopyNumber.getCopyNumber()) {
            uniquifiedAlleles.add(DEL_ALLELE);  //dels may be no-call
        }
        variantContextBuilder.alleles(uniquifiedAlleles);
        variantContextBuilder.attribute(VCFConstants.END_KEY, end);
        variantContextBuilder.genotypes(genotype);
        variantContextBuilder.log10PError(segment.getQualitySomeCalled()/-10.0);
        return variantContextBuilder.make();
    }

    private List<Allele> makeGenotypeAlleles(final int copyNumberCall, final int refCopyNumber, final Allele refAllele) {
        final List<Allele> returnAlleles = new ArrayList<>();
        final Allele genotypeAllele = getAlleleForCopyNumber(copyNumberCall, refCopyNumber, refAllele);
        //some allosomes like Y can have ref copy number zero, in which case we just no-call
        if (refCopyNumber == 0) {
            return GATKVariantContextUtils.noCallAlleles(1);
        }
        //for only one haplotype we know which allele it has
        if (refCopyNumber == 1) {
           return Arrays.asList(getAlleleForCopyNumber(copyNumberCall, refCopyNumber, refAllele));
        //can't determine counts per haplotypes if there is a duplication
        } else if (genotypeAllele.equals(DUP_ALLELE)) {
            return GATKVariantContextUtils.noCallAlleles(refCopyNumber);
        //for homDels, hetDels or homRefs
        } else if (refCopyNumber == 2) {
            returnAlleles.add(genotypeAllele);
            if (copyNumberCall == 0) {
                returnAlleles.add(genotypeAllele);
            } else {
                returnAlleles.add(refAllele);
            }
            return returnAlleles;
        //multiploid dels
        } else {
            for (int i = 0; i < copyNumberCall; i++) {
                returnAlleles.add(refAllele);
            }
            for (int i = copyNumberCall; i < refCopyNumber; i++) {
                returnAlleles.add(DEL_ALLELE);
            }
            return returnAlleles;
        }
    }

    /**
     *
     * @param copyNumberCall
     * @param refCopyNumber
     * @param refAllele
     * @return variant allele if copyNumberCall != refCopyNumber, else refAllele
     */
    private Allele getAlleleForCopyNumber(final int copyNumberCall, final int refCopyNumber, final Allele refAllele) {
        final Allele allele;
        if (copyNumberCall > refCopyNumber) {
            allele = DUP_ALLELE;
        } else if (copyNumberCall < refCopyNumber) {
            allele = DEL_ALLELE;
        } else {
            allele = refAllele;
        }
        return allele;
    }
}
