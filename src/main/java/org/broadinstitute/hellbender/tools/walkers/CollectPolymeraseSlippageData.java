package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.PolymeraseSlippageRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;

@CommandLineProgramProperties(
        summary = "Collect data on polymerase slippage in STRs",
        oneLineSummary = "Collect data on polymerase slippage in STRs",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class CollectPolymeraseSlippageData extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file")
    private File outputFile;

    private static final int MIN_STR_LENGTH = 5;
    private static final int MAX_REPEAT_UNIT_SIZE = 6;
    private static final int MIN_NUM_REPEATS = 3;
    private static final int MAX_NUM_REPEATS = 100;
    private static final int MAX_STR_SIZE = 30;
    private static final int REF_CACHE_SIZE = 1_000_000;


    private String lastContig = null;
    private byte[] refContig = null;
    private int endOfLastSTR = 0;

    final List<PolymeraseSlippageRecord> output = new ArrayList<>();

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final String contig = alignmentContext.getContig();
        final int position = (int) alignmentContext.getPosition();

        // TODO: check if contig == null?
        if (!contig.equals(lastContig)) {
            final SAMSequenceDictionary dict = getBestAvailableSequenceDictionary();
            final SimpleInterval wholeContig = new SimpleInterval(contig, 1, dict.getSequence(contig).getSequenceLength());
            refContig = referenceContext.getBases(wholeContig);
            endOfLastSTR = 0;
            lastContig = contig;
        } else {
            if (position <= endOfLastSTR) {
                return; // skip positions within last detected STR
            }
        }

        //TODO: this doesn't quite handle STRs with multipe possible repeat units, such as
        // TODO: TATTATTATT, which has three TATs starting from index 0 and three ATTs starting from 1

        int maxSTRLengthFound = 0;
        int bestUnitLength = 1;
        int bestNumRepeats = 1;
        for (int unitLength = 1; unitLength <= MAX_REPEAT_UNIT_SIZE; unitLength++) {
            final int numRepeats = countRepeats(refContig, position, unitLength);
            final int strLength = numRepeats * unitLength;
            if (strLength > maxSTRLengthFound && numRepeats >= MIN_NUM_REPEATS) {
                maxSTRLengthFound = strLength;
                bestUnitLength = unitLength;
                bestNumRepeats = numRepeats;
            }
        }

        if (maxSTRLengthFound < MIN_STR_LENGTH) {
            return;
        }

        final byte[] repeatUnit = Arrays.copyOfRange(refContig, position, position + bestUnitLength);

        final PolymeraseSlippageRecord record = new PolymeraseSlippageRecord(new String(repeatUnit), bestNumRepeats);

        for (final PileupElement pe : alignmentContext.getBasePileup()) {
            final GATKRead read = pe.getRead();
            final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), position, ReadUtils.ClippingTail.RIGHT_TAIL, true);
            if (offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED || offset < 0 || offset > read.getLength()) {
                continue;
            }

            final byte[] readBases = read.getBasesNoCopy();
            final byte[] readBasesFromPosition = Arrays.copyOfRange(readBases, offset + 1, Math.min(readBases.length, offset + MAX_STR_SIZE));

            final int numRepeats = GATKVariantContextUtils.findNumberOfRepetitions(repeatUnit, readBasesFromPosition, true);
            if (offset + repeatUnit.length * numRepeats + 1>= readBases.length) {
                continue; // uninformative -- read ends in STR
            }

            record.incrementCount(numRepeats);
        }
        output.add(record);
        endOfLastSTR = position + maxSTRLengthFound;
    }

    @Override
    public Object onTraversalSuccess() {
        PolymeraseSlippageRecord.writeToFile(output, outputFile);
        return "Success";
    }

    private int countRepeats(final byte[] ref, final int position, final int unitLength) {
        final int maxNumRepeats = Math.min((ref.length - position) / unitLength, MAX_NUM_REPEATS);

        int repeats;
        for (repeats = 1; repeats < maxNumRepeats; repeats++) {
            final int offset = position + repeats * unitLength;
            for (int n = 0; n < unitLength; n++) {
                if (ref[position + n] != ref[offset + n]) {
                    return repeats;
                }
            }
        }

        return repeats;
    }

}
