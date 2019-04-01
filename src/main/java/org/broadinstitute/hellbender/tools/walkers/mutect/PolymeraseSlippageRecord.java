package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

/**
 * Created by David Benjamin on 2/13/17.
 */
public class PolymeraseSlippageRecord {
    private final String repeatUnit;
    private final int numRefRepeats;

    private int refCount = 0;

    private int oneDeletionCount = 0;
    private int twoDeletionCount = 0;
    
    private int oneInsertionCount = 0;
    private int twoInsertionCount = 0;

    public PolymeraseSlippageRecord(final String repeatUnit, final int numRefRepeats) {
        this.repeatUnit = repeatUnit;
        this.numRefRepeats = numRefRepeats;
    }

    public String getRepeatUnit() { return repeatUnit; }

    public int getNumRefRepeats() { return numRefRepeats; }

    public int getTwoDeletionCount() { return twoDeletionCount; }

    public int getOneDeletionCount() { return oneDeletionCount; }

    public int getRefCount() { return refCount; }

    public int getOneInsertionCount() { return oneInsertionCount; }

    public int getTwoInsertionCount() { return twoInsertionCount; }

    public void incrementCount(final int numRepeats) {
           if (numRepeats == numRefRepeats) {
               refCount++;
           } else if (numRepeats == numRefRepeats - 1) {
               oneDeletionCount++;
           } else if (numRepeats == getNumRefRepeats() + 1) {
               oneInsertionCount++;
           } else if (numRepeats == numRefRepeats - 2) {
               twoDeletionCount++;
           } else if (numRepeats == numRefRepeats + 2) {
               twoInsertionCount++;
           }
    }


    //----- The following two public static methods read and write contamination files
    public static void writeToFile(final List<PolymeraseSlippageRecord> records, final File outputTable) {
        try ( PolymeraseSlippageTableWriter writer = new PolymeraseSlippageTableWriter(IOUtils.fileToPath(outputTable)) ) {
            writer.writeAllRecords(records);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable));
        }
    }

    public static List<PolymeraseSlippageRecord> readFromFile(final File tableFile) {
        try( PolymeraseSlippageTableReader reader = new PolymeraseSlippageTableReader(IOUtils.fileToPath(tableFile)) ) {
            return reader.toList();
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", tableFile));
        }
    }

    //-------- The following methods are boilerplate for reading and writing contamination tables
    private static class PolymeraseSlippageTableWriter extends TableWriter<PolymeraseSlippageRecord> {
        private PolymeraseSlippageTableWriter(final Path output) throws IOException {
            super(output, PolymeraseSlippageTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final PolymeraseSlippageRecord record, final DataLine dataLine) {
            dataLine.set(PolymeraseSlippageTableColumn.REPEAT_UNIT.toString(), record.getRepeatUnit())
                    .set(PolymeraseSlippageTableColumn.NUM_REPEATS.toString(), record.getNumRefRepeats())
                    .set(PolymeraseSlippageTableColumn.REF.toString(), record.getRefCount())
                    .set(PolymeraseSlippageTableColumn.ONE_DELETION.toString(), record.getOneDeletionCount())
                    .set(PolymeraseSlippageTableColumn.TWO_DELETION.toString(), record.getTwoDeletionCount())
                    .set(PolymeraseSlippageTableColumn.ONE_INSERTION.toString(), record.getOneInsertionCount())
                    .set(PolymeraseSlippageTableColumn.TWO_INSERTION.toString(), record.getTwoInsertionCount());
        }
    }

    private static class PolymeraseSlippageTableReader extends TableReader<PolymeraseSlippageRecord> {
        public PolymeraseSlippageTableReader(final Path path) throws IOException {
            super(path);
        }

        @Override
        protected PolymeraseSlippageRecord createRecord(final DataLine dataLine) {
            final String repeatUnit = dataLine.get(PolymeraseSlippageTableColumn.REPEAT_UNIT);
            final int numRepeats = dataLine.getInt(PolymeraseSlippageTableColumn.NUM_REPEATS);

            final int ref = dataLine.getInt(PolymeraseSlippageTableColumn.REF);
            final int oneDeletion = dataLine.getInt(PolymeraseSlippageTableColumn.ONE_DELETION);
            final int twoDeletion = dataLine.getInt(PolymeraseSlippageTableColumn.TWO_DELETION);
            final int oneInsertion = dataLine.getInt(PolymeraseSlippageTableColumn.ONE_DELETION);
            final int twoInsertion = dataLine.getInt(PolymeraseSlippageTableColumn.TWO_DELETION);

            final PolymeraseSlippageRecord result = new PolymeraseSlippageRecord(repeatUnit, numRepeats);
            result.refCount = ref;
            result.oneDeletionCount = oneDeletion;
            result.twoDeletionCount = twoDeletion;
            result.oneInsertionCount = oneInsertion;
            result.twoInsertionCount = twoInsertion;
            return result;
        }
    }

    private enum PolymeraseSlippageTableColumn {
        REPEAT_UNIT("unit"),
        NUM_REPEATS("repeats"),
        REF("ref"),
        ONE_DELETION("1del"),
        TWO_DELETION("2del"),
        ONE_INSERTION("1ins"),
        TWO_INSERTION("2ins");

        private final String columnName;

        PolymeraseSlippageTableColumn(final String columnName) {
            this.columnName = Utils.nonNull(columnName);
        }

        @Override
        public String toString() {
            return columnName;
        }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection(REPEAT_UNIT, NUM_REPEATS, REF, ONE_DELETION, TWO_DELETION, ONE_INSERTION, TWO_INSERTION);
    }
}
