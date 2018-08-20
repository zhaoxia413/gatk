package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Random;

public class XGBoostUtilsUnitTest  extends GATKBaseTest {
    private static final String SV_UTILS_TEST_DIR = toolsTestDir + "spark/sv/utils/";
    private static final String TEST_MATRIX_DATA_FILE = SV_UTILS_TEST_DIR + "agaricus-integers.csv.gz";
    private static final RealMatrix TEST_MATRIX = MachineLearningUtils.loadCsvFile(TEST_MATRIX_DATA_FILE);
    private static final XGBoostUtils.GATKXGBooster classifier = new XGBoostUtils.GATKXGBooster();
    private static final int NUM_TRAINING_ROUNDS = XGBoostUtils.DEFAULT_NUM_TRAINING_ROUNDS;
    private static final int EARLY_STOPPING_ROUNDS = XGBoostUtils.DEFAULT_EARLY_STOPPING_ROUNDS;
    private static final int NUM_CROSSVALIDATION_FOLDS = MachineLearningUtils.DEFAULT_NUM_CROSSVALIDATION_FOLDS;
    private static final int NUM_TUNING_ROUNDS = XGBoostUtils.DEFAULT_NUM_TUNING_ROUNDS; // keep tests quick
    private static final double TUNING_FRACTION = 1.0 / (1.0 + NUM_TUNING_ROUNDS);
    private static final boolean MAXIMIZE_EVAL_METRIC = false;
    private static final Map<String, Object> CLASSIFIER_PARAMS = XGBoostUtils.DEFAULT_CLASSIFIER_PARAMETERS;
    static {
        CLASSIFIER_PARAMS.put(XGBoostUtils.EVAL_METRIC_KEY, "logloss");
    }
    private static final Random random = new Random();

    private static final double MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY = 0.99;

    private static void assertLabelsEqual(final int[] actuals, final int[] expecteds, final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], "at index=" + index + ": " + message);
        }
    }

    @Test(groups = "sv")
    protected void testBasicTrain() throws IOException {
        // check that a classifier can be trained from data
        classifier.train(CLASSIFIER_PARAMS, TEST_MATRIX);

        // check that the classifier can make predictions on the same data;
        final int[] predictedLabels = classifier.predictClassLabels(TEST_MATRIX);

        // check that the predictions match the data (this data is easy to predict, they should)
        assertLabelsEqual(predictedLabels, MachineLearningUtils.getClassLabels(TEST_MATRIX),
                "Predicted labels not identical to actual labels");

        // predict probabilities of whole matrix
        final float[][] probabilities = classifier.predictProbability(TEST_MATRIX);
        // check that you get indentical results when predicting row-by-row
        for(int row = 0; row < TEST_MATRIX.getRowDimension(); ++row) {
            final float[] rowProbability = classifier.predictProbability(TEST_MATRIX.getRow(row));
            assertArrayEquals(rowProbability, probabilities[row], 0,
                    "Row " + row + ": different probabilities predicted for matrix and row-by-row");
        }

        // save classifier to temporary file
        File tempFile = File.createTempFile("gatk-xgboost-classifier", "kryo");
        classifier.save(tempFile.getAbsolutePath());
        // load classifier from temporary file

        final MachineLearningUtils.GATKClassifier loadedClassifier = MachineLearningUtils.GATKClassifier.load(
                tempFile.getAbsolutePath()
        );
        // check that you get identical results to first probability predictions
        final float[][] loadedProbabilities = loadedClassifier.predictProbability(TEST_MATRIX);
        assertMatrixEquals(loadedProbabilities, probabilities, 0.0,
                "Probabilities predicted by loaded classifier not equal to original");
    }

    @Test(groups = "sv")
    protected void testGetTrainingTrace() {
        final int[] stratify = MachineLearningUtils.getClassLabels(TEST_MATRIX);

        final MachineLearningUtils.TrainTestSplit hyperSplit = MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                TUNING_FRACTION, TEST_MATRIX.getRowDimension(), random, stratify
        );
        final RealMatrix trainMatrix = MachineLearningUtils.sliceRows(TEST_MATRIX, hyperSplit.trainRows);
        final RealMatrix validateMatrix = MachineLearningUtils.sliceRows(TEST_MATRIX, hyperSplit.testRows);

        final float[] trainingTrace = classifier.trainAndReturnQualityTrace(CLASSIFIER_PARAMS, trainMatrix, validateMatrix,
                NUM_TRAINING_ROUNDS, EARLY_STOPPING_ROUNDS, MAXIMIZE_EVAL_METRIC);
        Assert.assertEquals(
                trainingTrace.length, NUM_TRAINING_ROUNDS,
                "Training trace did not have requested number of rounds (" + trainingTrace.length
                        + " instead of " + NUM_TRAINING_ROUNDS
        );
    }

    @Test(groups = "sv")
    protected void testCrossvalidatedTuneAndTrain() {
        final int[] stratify = MachineLearningUtils.getClassLabels(TEST_MATRIX);

        final MachineLearningUtils.TrainTestSplit hyperSplit = MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                TUNING_FRACTION, TEST_MATRIX.getRowDimension(), random, stratify
        );
        final RealMatrix tuneMatrix = MachineLearningUtils.sliceRows(TEST_MATRIX, hyperSplit.trainRows);
        final int[] tuneStratify = MachineLearningUtils.slice(stratify, hyperSplit.trainRows);
        final RealMatrix validateMatrix = MachineLearningUtils.sliceRows(TEST_MATRIX, hyperSplit.testRows);
        final int[] validateStratify = MachineLearningUtils.slice(stratify, hyperSplit.testRows);

        // set the number of threads
        classifier.chooseNumThreads(CLASSIFIER_PARAMS, XGBoostUtils.NUM_THREADS_KEY, tuneMatrix);

        final Map<String, Object> bestClassifierParameters = classifier.tuneClassifierParameters(
                CLASSIFIER_PARAMS, XGBoostUtils.DEFAULT_TUNING_PARAMETERS,
                MachineLearningUtils.ClassifierTuningStrategy.RANDOM,
                tuneMatrix, random, tuneStratify, NUM_CROSSVALIDATION_FOLDS, NUM_TRAINING_ROUNDS,
                EARLY_STOPPING_ROUNDS, NUM_TUNING_ROUNDS, MAXIMIZE_EVAL_METRIC
        );

        final int[] predictedTestLabels = classifier.crossvalidatePredict(
                validateMatrix, bestClassifierParameters, random, validateStratify, NUM_CROSSVALIDATION_FOLDS
        );

        final double accuracy = MachineLearningUtils.getPredictionAccuracy(
                predictedTestLabels, MachineLearningUtils.getClassLabels(validateMatrix)
        );

        Assert.assertTrue(accuracy >= MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY,
                "Crossvalidated prediction accuracy (" + accuracy + ") less than passing ("
                        + MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY + ")");
    }

    private static void assertArrayEquals(final float[] actuals, final float[] expecteds, final double tol,
                                          final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], tol, "at index=" + index + ": " + message);
        }
    }

    private static void assertMatrixEquals(final float[][] actuals, final float[][] expecteds, final double tol,
                                           final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Number of rows not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            assertArrayEquals(actuals[index], expecteds[index], tol, "at row=" + index + ": " + message);
        }
    }
}
