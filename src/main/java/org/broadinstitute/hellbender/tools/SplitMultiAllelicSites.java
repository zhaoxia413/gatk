package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;

public class SplitMultiAllelicSites extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    public File outFile = null;

    private VariantContextWriter vcfWriter = null;
    public static final double SUM_GL_THRESH_NOCALL = -0.1; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.


    public enum GenotypeAssignmentMethod {
        /**
         * set all of the genotype GT values to NO_CALL
         */
        SET_TO_NO_CALL,

        /**
         * set all of the genotype GT values to NO_CALL and remove annotations
         */
        SET_TO_NO_CALL_NO_ANNOTATIONS,


        /**
         * Use the subsetted PLs to greedily assigned genotypes
         */
        USE_PLS_TO_ASSIGN,

        /**
         * Try to match the original GT calls, if at all possible
         *
         * Suppose I have 3 alleles: A/B/C and the following samples:
         *
         *       original_GT best_match to A/B best_match to A/C
         * S1 => A/A A/A A/A
         * S2 => A/B A/B A/A
         * S3 => B/B B/B A/A
         * S4 => B/C A/B A/C
         * S5 => C/C A/A C/C
         *
         * Basically, all alleles not in the subset map to ref.  It means that het-alt genotypes
         * when split into 2 bi-allelic variants will be het in each, which is good in some cases,
         * rather than the undetermined behavior when using the PLs to assign, which could result
         * in hom-var or hom-ref for each, depending on the exact PL values.
         */
        BEST_MATCH_TO_ORIGINAL,

        /**
         * do not even bother changing the GTs
         */
        DO_NOT_ASSIGN_GENOTYPES
    }

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outFile);
        vcfWriter.writeHeader(getHeaderForVariants());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if(variant.isBiallelic()) {
            vcfWriter.add(variant);
            return;
        }

        final List<VariantContext> biallelics = new LinkedList<>();

        for ( final Allele alt : variant.getAlternateAlleles() ) {
            final VariantContextBuilder builder = new VariantContextBuilder(variant);

            // make biallelic alleles
            final List<Allele> alleles = Arrays.asList(variant.getReference(), alt);
            builder.alleles(alleles);

            // since the VC has been subset, remove the invalid attributes
            for ( final String key : variant.getAttributes().keySet() ) {
                if ( !(key.equals(VCFConstants.ALLELE_COUNT_KEY) || key.equals(VCFConstants.ALLELE_FREQUENCY_KEY) || key.equals(VCFConstants.ALLELE_NUMBER_KEY))) {
                    builder.rmAttribute(key);
                }
            }

            // subset INFO field annotations if available if genotype is called
            addInfoFiledAnnotations(variant, builder, alt, true);

                    builder.genotypes(subsetAlleles(variant, alleles, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL));
                    final VariantContext trimmed = GATKVariantContextUtils.trimAlleles(builder.make(), false, true);
                    biallelics.add(trimmed);
        }

        for(VariantContext vc : biallelics) {
            vcfWriter.add(vc);
        }
    }

    /**
     * Subset the Variant Context to the specific set of alleles passed in (pruning the PLs appropriately)
     *
     * @param vc                 variant context with genotype likelihoods
     * @param allelesToUse       which alleles from the vc are okay to use; *** must be in the same relative order as those in the original VC ***
     * @param assignGenotypes    assignment strategy for the (subsetted) PLs
     * @return a new non-null GenotypesContext with subsetted alleles
     */
    public static GenotypesContext subsetAlleles(final VariantContext vc,
                                                 final List<Allele> allelesToUse,
                                                 final GenotypeAssignmentMethod assignGenotypes) {
        if ( vc == null ) throw new IllegalArgumentException("the VariantContext cannot be null");
        if ( allelesToUse == null ) throw new IllegalArgumentException("the alleles to use cannot be null");
        if ( allelesToUse.isEmpty() ) throw new IllegalArgumentException("must have alleles to use");
        if ( allelesToUse.get(0).isNonReference() ) throw new IllegalArgumentException("First allele must be the reference allele");
        if ( allelesToUse.size() == 1 ) throw new IllegalArgumentException("Cannot subset to only 1 alt allele");

        // optimization: if no input genotypes, just exit
        if (vc.getGenotypes().isEmpty()) return GenotypesContext.create();

        // find the likelihoods indexes to use from the used alternate alleles
        final List<List<Integer>> likelihoodIndexesToUse = determineLikelihoodIndexesToUse(vc, allelesToUse);

        // find the strand allele count indexes to use from the used alternate alleles
        final List<Integer> sacIndexesToUse = determineSACIndexesToUse(vc, allelesToUse);

        // create the new genotypes
        return createGenotypesWithSubsettedLikelihoods(vc.getGenotypes(), vc, allelesToUse, likelihoodIndexesToUse, sacIndexesToUse, assignGenotypes);
    }

    /**
     * Create the new GenotypesContext with the subsetted PLs, SACs and ADs
     *
     * @param originalGs               the original GenotypesContext
     * @param originalVC               the original VariantContext
     * @param allelesToUse             the actual alleles to use with the new Genotypes
     * @param likelihoodIndexesToUse   the indexes in the PL to use given the allelesToUse for each genotype (@see #determineLikelihoodIndexesToUse())
     * @param sacIndexesToUse          the indexes in the SAC to use given the allelesToUse (@see #determineSACIndexesToUse())
     * @param assignGenotypes          assignment strategy for the (subsetted) PLs
     * @return a new non-null GenotypesContext
     */
    private static GenotypesContext createGenotypesWithSubsettedLikelihoods(final GenotypesContext originalGs,
                                                                            final VariantContext originalVC,
                                                                            final List<Allele> allelesToUse,
                                                                            final List<List<Integer>> likelihoodIndexesToUse,
                                                                            final List<Integer> sacIndexesToUse,
                                                                            final GenotypeAssignmentMethod assignGenotypes) {

        if ( originalGs == null ) throw new IllegalArgumentException("the original GenotypesContext cannot be null");
        if ( originalVC == null ) throw new IllegalArgumentException("the original VariantContext cannot be null");
        if ( allelesToUse == null ) throw new IllegalArgumentException("the alleles to use cannot be null");

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());

        // the samples
        final List<String> sampleIndices = originalGs.getSampleNamesOrderedByName();

        // create the new genotypes
        for ( int k = 0; k < originalGs.size(); k++ ) {
            final Genotype g = originalGs.get(sampleIndices.get(k));
            final GenotypeBuilder gb = new GenotypeBuilder(g);

            // create the new likelihoods array from the used alleles
            double[] newLikelihoods;
            if ( !g.hasLikelihoods() ) {
                // we don't have any likelihoods, so we null out PLs and make G ./.
                newLikelihoods = null;
                gb.noPL();
            } else {
                // make sure we are seeing the expected number of likelihoods per sample
                final int expectedNumLikelihoods = GenotypeLikelihoods.numLikelihoods(originalVC.getNAlleles(), g.getPloidy());
                final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
                if ( likelihoodIndexesToUse == null ) {
                    newLikelihoods = originalLikelihoods;
                } else if ( originalLikelihoods.length != expectedNumLikelihoods ) {
                    newLikelihoods = null;
                } else {
                    newLikelihoods = new double[likelihoodIndexesToUse.get(k).size()];
                    int newIndex = 0;
                    for ( final int oldIndex : likelihoodIndexesToUse.get(k) )
                        newLikelihoods[newIndex++] = originalLikelihoods[oldIndex];

                    // might need to re-normalize
                    newLikelihoods = MathUtils.normalizeLog10(newLikelihoods, false, true);
                }

                if ( newLikelihoods == null || (originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) == 0 && likelihoodsAreUninformative(newLikelihoods) )) {
                    gb.noPL();
                    gb.noGQ();
                } else {
                    gb.PL(newLikelihoods);
                    final int PLindex = MathUtils.maxElementIndex(newLikelihoods);
                    gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
                }

            }

            // create the new strand allele counts array from the used alleles
            if ( g.hasExtendedAttribute(VCFConstants.STRAND_BIAS_KEY)){
                int[] newSACs = makeNewSACs(g, sacIndexesToUse);
                gb.attribute(VCFConstants.STRAND_BIAS_KEY, newSACs);
            }

            updateGenotypeAfterSubsetting(g.getAlleles(), g.getPloidy(), gb, assignGenotypes, newLikelihoods, allelesToUse);
            newGTs.add(gb.make());
        }

        return fixADFromSubsettedAlleles(newGTs, originalVC, allelesToUse);
    }

    /**
     * Fix the AD for the GenotypesContext of a VariantContext that has been subset
     *
     * @param originalGs       the original GenotypesContext
     * @param originalVC       the original VariantContext
     * @param allelesToUse     the new (sub)set of alleles to use
     * @return a new non-null GenotypesContext
     */
    public static GenotypesContext fixADFromSubsettedAlleles(final GenotypesContext originalGs, final VariantContext originalVC, final List<Allele> allelesToUse) {
        if (originalGs == null) throw new IllegalArgumentException("the original Gs cannot be null");
        if (originalVC == null) throw new IllegalArgumentException("the original VC cannot be null");
        if (allelesToUse == null) throw new IllegalArgumentException("the alleles to use list cannot be null");

        // the bitset representing the allele indexes we want to keep
        final BitSet alleleIndexesToUse = getAlleleIndexBitset(originalVC, allelesToUse);

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());

        // the samples
        final List<String> sampleIndices = originalGs.getSampleNamesOrderedByName();

        // create the new genotypes
        for ( int k = 0; k < originalGs.size(); k++ ) {
            final Genotype g = originalGs.get(sampleIndices.get(k));
            newGTs.add(fixAD(g, alleleIndexesToUse));
        }

        return newGTs;
    }

    /**
     * Fix the AD for the given Genotype
     *
     * @param genotype              the original Genotype
     * @param alleleIndexesToUse    a bitset describing whether or not to keep a given index
     * @return a non-null Genotype
     */
    private static Genotype fixAD(final Genotype genotype, final BitSet alleleIndexesToUse) {
        // if it ain't broke don't fix it
        if ( !genotype.hasAD() )
            return genotype;

        final GenotypeBuilder builder = new GenotypeBuilder(genotype);

        final int[] oldAD = genotype.getAD();
        final int[] newAD = new int[alleleIndexesToUse.cardinality()];

        int currentIndex = 0;
        for ( int i = alleleIndexesToUse.nextSetBit(0); i >= 0; i = alleleIndexesToUse.nextSetBit(i+1) ) {
            if ( i >= oldAD.length ) {
                throw new IllegalStateException("AD has " + oldAD.length + " items. It should have at least " + (i+1) + ".");
            }

            newAD[currentIndex++] = oldAD[i];
        }

        return builder.AD(newAD).make();
    }

    /**
     * Add the genotype call (GT) field to GenotypeBuilder using the requested algorithm assignmentMethod
     *
     * @param originalGT the original genotype calls, cannot be null
     * @param ploidy the number of sets of chromosomes
     * @param gb the builder where we should put our newly called alleles, cannot be null
     * @param assignmentMethod the method to use to do the assignment, cannot be null
     * @param newLikelihoods a vector of likelihoods to use if the method requires PLs, should be log10 likelihoods, cannot be null
     * @param allelesToUse the alleles we are using for our subsetting
     */
    public static void updateGenotypeAfterSubsetting(final List<Allele> originalGT,
                                                     final int ploidy,
                                                     final GenotypeBuilder gb,
                                                     final GenotypeAssignmentMethod assignmentMethod,
                                                     final double[] newLikelihoods,
                                                     final List<Allele> allelesToUse) {
        if ( originalGT == null ) throw new IllegalArgumentException("originalGT cannot be null");
        if ( gb == null ) throw new IllegalArgumentException("gb cannot be null");
        if ( allelesToUse.isEmpty() || allelesToUse == null ) throw new IllegalArgumentException("allelesToUse cannot be empty or null");

        switch ( assignmentMethod ) {
            case DO_NOT_ASSIGN_GENOTYPES:
                break;
            case SET_TO_NO_CALL:
                gb.alleles(noCallAlleles(ploidy));
                gb.noGQ();
                break;
            case SET_TO_NO_CALL_NO_ANNOTATIONS:
                gb.alleles(noCallAlleles(ploidy));
                gb.noGQ();
                gb.noAD();
                gb.noPL();
                gb.noAttributes();
                break;
            case USE_PLS_TO_ASSIGN:
                if ( newLikelihoods == null || likelihoodsAreUninformative(newLikelihoods) ) {
                    // if there is no mass on the (new) likelihoods, then just no-call the sample
                    gb.alleles(noCallAlleles(ploidy));
                    gb.noGQ();
                } else {
                    // find the genotype with maximum likelihoods
                    final int PLindex = MathUtils.maxElementIndex(newLikelihoods);
                    final List<Allele> alleles = new ArrayList<>();
                    for ( final Integer alleleIndex : GenotypeLikelihoods.getAlleles(PLindex, ploidy)) {
                        alleles.add(allelesToUse.get(alleleIndex) );
                    }
                    gb.alleles(alleles);
                    gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
                }
                break;
            case BEST_MATCH_TO_ORIGINAL:
                final List<Allele> best = new LinkedList<>();
                final Allele ref = allelesToUse.get(0);
                for ( final Allele originalAllele : originalGT ) {
                    best.add((allelesToUse.contains(originalAllele) || originalAllele.isNoCall()) ? originalAllele : ref);
                }
                gb.alleles(best);
                break;
        }
    }

    /**
     * Returns a {@link Allele#NO_CALL NO_CALL} allele list provided the ploidy.
     *
     * @param ploidy the required ploidy.
     *
     * @return never {@code null}, but an empty list if {@code ploidy} is equal or less than 0. The returned list
     *   might or might not be mutable.
     */
    public static List<Allele> noCallAlleles(final int ploidy) {
        if (NOCALL_LISTS.length <= ploidy)
            ensureNoCallListsCapacity(ploidy);
        return NOCALL_LISTS[ploidy];
    }

    /**
     * Cached NO_CALL immutable lists where the position ith contains the list with i elements.
     */
    private static List<Allele>[] NOCALL_LISTS = new List[] {
            Collections.emptyList(),
            Collections.singletonList(Allele.NO_CALL),
            Collections.nCopies(2,Allele.NO_CALL)
    };

    /**
     * Synchronized code to ensure that {@link #NOCALL_LISTS} has enough entries beyod the requested ploidy
     * @param capacity the requested ploidy.
     */
    private static synchronized void ensureNoCallListsCapacity(final int capacity) {
        final int currentCapacity = NOCALL_LISTS.length - 1;
        if (currentCapacity >= capacity)
            return;
        NOCALL_LISTS = Arrays.copyOf(NOCALL_LISTS,Math.max(capacity,currentCapacity << 1) + 1);
        for (int i = currentCapacity + 1; i < NOCALL_LISTS.length; i++)
            NOCALL_LISTS[i] = Collections.nCopies(i,Allele.NO_CALL);
    }



    /**
     * Make a new SAC array from the a subset of the genotype's original SAC
     *
     * @param g               the genotype
     * @param sacIndexesToUse the indexes in the SAC to use given the allelesToUse (@see #determineSACIndexesToUse())
     * @return subset of SACs from the original genotype, the original SACs if sacIndexesToUse is null
     */
    public static int[] makeNewSACs(final Genotype g, final List<Integer> sacIndexesToUse) {

        if (g == null) throw new IllegalArgumentException("the genotype cannot be null");

        final int[] oldSACs  = getSACs(g);

        if (sacIndexesToUse == null) {
            return oldSACs;
        } else {
            final int[] newSACs = new int[sacIndexesToUse.size()];
            int newIndex = 0;
            for (final int oldIndex : sacIndexesToUse) {
                newSACs[newIndex++] = oldSACs[oldIndex];
            }
            return newSACs;
        }
    }

    /**
     * Get the genotype SACs
     *
     * @param g the genotype
     * @return an arrays of SACs
     * @throws UserException if the type of the SACs is unexpected
     */
    private static int[] getSACs(final Genotype g) {

        if ( g == null ) throw new IllegalArgumentException("the Genotype cannot be null");
        if ( !g.hasExtendedAttribute(VCFConstants.STRAND_BIAS_KEY) )
            throw new IllegalArgumentException("Genotype must have SAC");

        if ( g.getExtendedAttributes().get(VCFConstants.STRAND_BIAS_KEY).getClass().equals(String.class) ) {
            final String SACsString = (String) g.getExtendedAttributes().get(VCFConstants.STRAND_BIAS_KEY);
            ArrayList<String> stringSACs = new ArrayList<>(Arrays.asList(SACsString.split(",")));
            final int[] intSACs = new int[stringSACs.size()];
            int i = 0;
            for (String sac : stringSACs)
                intSACs[i++] = Integer.parseInt(sac);

            return intSACs;
        }
        else if ( g.getExtendedAttributes().get(VCFConstants.STRAND_BIAS_KEY).getClass().equals(int[].class) )
            return (int[]) g.getExtendedAttributes().get(VCFConstants.STRAND_BIAS_KEY);
        else
            throw new UserException("Unexpected SAC type");
    }

    private static boolean likelihoodsAreUninformative(final double[] likelihoods) {
        return MathUtils.sum(likelihoods) > SUM_GL_THRESH_NOCALL;
    }

    /**
     * Find the strand allele count indexes to use for a selected set of alleles
     *
     * @param originalVC   the original VariantContext
     * @param allelesToUse the subset of alleles to use
     * @return a list of SAC indexes to use or null if none
     */
    public static List<Integer>  determineSACIndexesToUse(final VariantContext originalVC, final List<Allele> allelesToUse) {

        if ( originalVC == null ) throw new IllegalArgumentException("the original VC cannot be null");
        if ( allelesToUse == null ) throw new IllegalArgumentException("the alleles to use cannot be null");

        // the bitset representing the allele indexes we want to keep
        final BitSet alleleIndexesToUse = getAlleleIndexBitset(originalVC, allelesToUse);

        // an optimization: if we are supposed to use all (or none in the case of a ref call) of the alleles,
        // then we can keep the SACs as is; otherwise, we determine which ones to keep
        if (alleleIndexesToUse.cardinality() == alleleIndexesToUse.size())
            return null;

        return getSACIndexes(alleleIndexesToUse);
    }

    /**
     * Get the actual strand aleele counts indexes to use given the corresponding allele indexes
     *
     * @param alleleIndexesToUse    the bitset representing the alleles to use (@see #getAlleleIndexBitset)
     * @return a non-null List
     */
    private static List<Integer> getSACIndexes(final BitSet alleleIndexesToUse) {

        if (alleleIndexesToUse == null) throw new IllegalArgumentException("the alleles to use cannot be null");
        if (alleleIndexesToUse.isEmpty()) throw new IllegalArgumentException("cannot have no alleles to use");

        final List<Integer> result = new ArrayList<>(2 * alleleIndexesToUse.size());

        for (int SACindex = 0; SACindex < alleleIndexesToUse.size(); SACindex++) {
            if (alleleIndexesToUse.get(SACindex)) {
                result.add(2 * SACindex);
                result.add(2 * SACindex + 1);
            }
        }

        return result;
    }

    /**
     * Find the likelihood indexes to use for a selected set of alleles
     *
     * @param originalVC        the original VariantContext
     * @param allelesToUse      the subset of alleles to use
     * @return a list of PL indexes to use or null if none
     */
    private static List<List<Integer>> determineLikelihoodIndexesToUse(final VariantContext originalVC, final List<Allele> allelesToUse) {

        if ( originalVC == null) throw new IllegalArgumentException("the original VariantContext cannot be null");
        if ( allelesToUse == null ) throw new IllegalArgumentException("the alleles to use cannot be null");

        // the bitset representing the allele indexes we want to keep
        final BitSet alleleIndexesToUse = getAlleleIndexBitset(originalVC, allelesToUse);

        // an optimization: if we are supposed to use all (or none in the case of a ref call) of the alleles,
        // then we can keep the PLs as is; otherwise, we determine which ones to keep
        if ( alleleIndexesToUse.cardinality() == alleleIndexesToUse.size() )
            return null;

        return getLikelihoodIndexes(originalVC, alleleIndexesToUse);
    }

    /**
     * Get the actual likelihoods indexes to use given the corresponding diploid allele indexes
     *
     * @param originalVC           the original VariantContext
     * @param alleleIndexesToUse   the bitset representing the alleles to use (@see #getAlleleIndexBitset)
     * @return likelihoods indexes for each genotype
     */
    private static List<List<Integer>> getLikelihoodIndexes(final VariantContext originalVC, BitSet alleleIndexesToUse) {

        final List<List<Integer>> likelihoodIndexesPerGenotype = new ArrayList<List<Integer>>(10);

        for (final Genotype g : originalVC.getGenotypes()) {
            final int numLikelihoods = GenotypeLikelihoods.numLikelihoods(originalVC.getNAlleles(), g.getPloidy());
            final List<Integer> likelihoodIndexes = new ArrayList<>(30);
            for ( int PLindex = 0; PLindex < numLikelihoods; PLindex++ ) {
                GenotypeLikelihoods.initializeAnyploidPLIndexToAlleleIndices(originalVC.getNAlleles() - 1, g.getPloidy());
                // consider this entry only if all the alleles are good
                if ( GenotypeLikelihoods.getAlleles(PLindex, g.getPloidy()).stream().allMatch(i -> alleleIndexesToUse.get(i)) )
                    likelihoodIndexes.add(PLindex);
            }
            likelihoodIndexesPerGenotype.add(likelihoodIndexes);
        }

        return likelihoodIndexesPerGenotype;
    }

    /**
     * Add the VCF INFO field annotations for the used alleles
     *
     * @param vc                original variant context
     * @param builder           variant context builder with subset of original variant context's alleles
     * @param altAllele         alternate allele
     * @param keepOriginalChrCounts keep the orignal chromosome counts before subsetting
     * @return variant context builder with updated INFO field attribute values
     */
    private static void addInfoFiledAnnotations(final VariantContext vc, final VariantContextBuilder builder,  final Allele altAllele,
                                                final boolean keepOriginalChrCounts){

        if (vc == null) throw new IllegalArgumentException("the variant context cannot be null");
        if (builder == null) throw new IllegalArgumentException("the variant context builder cannot be null");
        if (builder.getAlleles() == null) throw new IllegalArgumentException("the variant context builder alleles cannot be null");

        final List<Allele> alleles = builder.getAlleles();
        if (alleles.size() < 2) throw new IllegalArgumentException("the variant context builder must contain at least 2 alleles");

        // don't have to subset, the original vc has the same number and hence, the same alleles
        boolean keepOriginal = ( vc.getAlleles().size() == builder.getAlleles().size() );

        if ( keepOriginalChrCounts ) {
            if (vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY))
                builder.attribute(VCFConstants.ALLELE_COUNT_KEY, keepOriginal ?
                        vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY) : getAltAlleleInfoFieldValue(VCFConstants.ALLELE_COUNT_KEY, vc, altAllele));
            if (vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY))
                builder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, keepOriginal ?
                        vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) :  getAltAlleleInfoFieldValue(VCFConstants.ALLELE_FREQUENCY_KEY, vc, altAllele));
            if (vc.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY)) {
                builder.attribute(VCFConstants.ALLELE_NUMBER_KEY, vc.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
            }
        }

        VariantContextUtils.calculateChromosomeCounts(builder, true);
    }

    /**
     * Get the alternate allele INFO field value
     *
     * @param infoFieldName     VCF INFO field name
     * @param vc                variant context
     * @param altAllele         the alternate allele
     * @return alternate allele INFO field value
     */
    private static Object getAltAlleleInfoFieldValue(final String infoFieldName, final VariantContext vc, final Allele altAllele) {

        //final String[] splitOriginalField = vc.getAttribute(infoFieldName).toString().split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
        final Object[] splitOriginalField = getVAttributeValues(vc.getAttribute(infoFieldName));

        // subset the field
        final BitSet alleleIndexesToUse = getAlleleIndexBitset(vc, Arrays.asList(altAllele));

        // skip the first allele, which is the reference
        for (int i = 1; i < alleleIndexesToUse.size(); i++) {
            if (alleleIndexesToUse.get(i))
                return splitOriginalField[i-1];
        }

        throw new UserException("Alternate allele " + altAllele.toString() + " not in Variant Context " + vc.toString());
    }

    /**
     * Pulls out the appropriate values for the INFO field attribute
     *
     * @param attribute    INFO field attribute
     * @return tokenized attribute values
     */
    private static Object[] getVAttributeValues(final Object attribute) {

        if (attribute == null) throw new IllegalArgumentException("the attribute cannot be null");

        // break the original attributes into separate tokens
        final Object[] tokens;
        if ( attribute.getClass().isArray() )
            tokens = (Object[])attribute;
        else if ( List.class.isAssignableFrom(attribute.getClass()) )
            tokens = ((List)attribute).toArray();
        else
            tokens = attribute.toString().split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);

        return tokens;
    }

    /**
     * Given an original VariantContext and a list of alleles from that VC to keep,
     * returns a bitset representing which allele indexes should be kept
     *
     * @param originalVC   the original VC
     * @param allelesToUse the list of alleles to keep
     * @return non-null bitset
     */
    private static BitSet getAlleleIndexBitset(final VariantContext originalVC, final List<Allele> allelesToUse) {

        if (originalVC == null) throw new IllegalArgumentException("the original VC cannot be null");
        if (allelesToUse == null) throw new IllegalArgumentException("the alleles to use cannot be null");

        final int numOriginalAltAlleles = originalVC.getNAlleles() - 1;
        final BitSet alleleIndexesToKeep = new BitSet(numOriginalAltAlleles + 1);

        // the reference Allele is always used
        alleleIndexesToKeep.set(0);
        for (int i = 0; i < numOriginalAltAlleles; i++) {
            if (allelesToUse.contains(originalVC.getAlternateAllele(i)))
                alleleIndexesToKeep.set(i + 1);
        }

        return alleleIndexesToKeep;
    }



    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }


}
