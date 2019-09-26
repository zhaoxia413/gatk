package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.util.QualityUtil;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

public final class FragmentUtils {
    private FragmentUtils() {}

    public final static double DEFAULT_PCR_SNV_ERROR_RATE = 1e-4;
    public final static int DEFAULT_PCR_SNV_ERROR_QUAL = QualityUtil.getPhredScoreFromErrorProbability(DEFAULT_PCR_SNV_ERROR_RATE);
    public final static int HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL = DEFAULT_PCR_SNV_ERROR_QUAL / 2;

    /**
     * Fix two overlapping reads from the same fragment by adjusting base qualities, if possible
     *
     *  Looks at the bases and alignment, and tries its best to create adjusted base qualities so that the observations
     * are not treated independently.  Sets the qualities of firstRead and secondRead to mimic a merged read or
     * nothing if the algorithm cannot create a meaningful one
     * @param pair two overlapping paired reads
     * @param pcrSnvQual phred-scaled quality of substitution errors from PCR. May not be negative.
     * @param pcrIndelQual phred-scaled quality of indel errors from PCR. May not be negative.
     */
    public static void adjustQualsOfOverlappingPairedFragments(final Pair<GATKRead, GATKRead> pair, final int pcrSnvQual, final int pcrIndelQual) {
        final boolean inOrder = pair.getLeft().getSoftStart() < pair.getRight().getSoftStart();
        final GATKRead firstRead = inOrder ? pair.getLeft() : pair.getRight();
        final GATKRead secondRead = inOrder ? pair.getRight() : pair.getLeft();

        Utils.nonNull(firstRead);
        Utils.nonNull(secondRead);
        Utils.validateArg(firstRead.getName().equals(secondRead.getName()), () ->
                "attempting to merge two reads with different names " + firstRead + " and " + secondRead);

        // don't adjust fragments that do not overlap
        if (firstRead.getEnd() < secondRead.getStart() || !firstRead.getContig().equals(secondRead.getContig())) {
            return;
        }

        final Pair<Integer, Boolean> offset = ReadUtils.getReadCoordinateForReferenceCoordinate(firstRead, secondRead.getStart());
        final int firstReadStop = (offset.getRight() ? offset.getLeft() + 1 : offset.getLeft());
        final int numOverlappingBases = Math.min(firstRead.getLength() - firstReadStop, secondRead.getLength());

        final byte[] firstReadBases = firstRead.getBases();
        final byte[] firstReadQuals = firstRead.getBaseQualities();
        final byte[] secondReadBases = secondRead.getBases();
        final byte[] secondReadQuals = secondRead.getBaseQualities();

        final int halfOfPcrErrorQual = pcrSnvQual / 2;
        final int halfOfPcrIndelQual = pcrIndelQual / 2;

        for (int i = 0; i < numOverlappingBases; i++) {

            final int firstReadIndex = firstReadStop + i;
            final byte firstReadBase = firstReadBases[firstReadIndex];
            final byte secondReadBase = secondReadBases[i];

            if (firstReadBase == secondReadBase) {
                firstReadQuals[firstReadIndex] = (byte) Math.min(firstReadQuals[firstReadIndex], halfOfPcrErrorQual);
                secondReadQuals[i] = (byte) Math.min(secondReadQuals[i], halfOfPcrErrorQual);
            }
        }
        firstRead.setBaseQualities(firstReadQuals);
        secondRead.setBaseQualities(secondReadQuals);


        final byte[] firstReadInsertionQuals = ReadUtils.getBaseInsertionQualities(firstRead);
        final byte[] firstReadDeletionQuals = ReadUtils.getBaseInsertionQualities(firstRead);
        final byte[] secondReadInsertionQuals = ReadUtils.getBaseInsertionQualities(secondRead);
        final byte[] secondReadDeletionQuals = ReadUtils.getBaseInsertionQualities(secondRead);

        for (int i = 0; i < numOverlappingBases; i++) {
            final int firstReadIndex = firstReadStop + i;
            firstReadDeletionQuals[firstReadIndex] = (byte) Math.min(firstReadDeletionQuals[firstReadIndex], halfOfPcrIndelQual);
            firstReadInsertionQuals[firstReadIndex] = (byte) Math.min(firstReadInsertionQuals[firstReadIndex], halfOfPcrIndelQual);
            secondReadDeletionQuals[i] = (byte) Math.min(secondReadDeletionQuals[i], halfOfPcrIndelQual);
            secondReadInsertionQuals[i] = (byte) Math.min(secondReadInsertionQuals[i], halfOfPcrIndelQual);
        }

        ReadUtils.setDeletionBaseQualities(firstRead, firstReadDeletionQuals);
        ReadUtils.setInsertionBaseQualities(firstRead, firstReadInsertionQuals);
        ReadUtils.setDeletionBaseQualities(secondRead, secondReadDeletionQuals);
        ReadUtils.setInsertionBaseQualities(secondRead, secondReadInsertionQuals);
    }
}
