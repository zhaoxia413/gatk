package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.KryoSerializable;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class MachineLearningUtils {
    public static final String NUM_TRAINING_ROUNDS_KEY = "num_training_rounds";
    public static final int DEFAULT_NUM_CROSSVALIDATION_FOLDS = 5;
    public static final int DEFAULT_NUM_TUNING_ROUNDS = 100;

    private static final int NUM_CALIBRATION_TRAINING_ROUNDS = 500;
    private static final int NUM_CALIBRATION_CLASS_ROWS = 500;
    public static final int CLASS_LABEL_COLUMN = 0;
    private static final Logger localLogger = LogManager.getLogger(MachineLearningUtils.class);

    public static Array2DRowRealMatrix loadCsvFile(final String filename) {
        return loadCsvFile(filename, ",", "#", CLASS_LABEL_COLUMN);
    }

    /**
     * Load matrix from (possibly gzipped) csv file. Note: the first column of said
     * @param filename
     * @return
     */
    public static Array2DRowRealMatrix loadCsvFile(final String filename, final String delimiter,
                                                   final String commentCharacter, final int classLabelsColumn) {
        final int numColumns;
        final List<double[]> rowsList;
        try (final BufferedReader reader = IOUtil.openFileForBufferedReading(new File(filename))) {
            rowsList = new ArrayList<>();
            if(!reader.ready()) {
                throw new GATKException("Unable to read matrix from " + filename);
            }
            rowsList.add(getNextCsvFeaturesRow(reader, delimiter, commentCharacter, classLabelsColumn));
            numColumns = rowsList.get(0).length;
            while(reader.ready()) {
                final double[] features = getNextCsvFeaturesRow(reader, delimiter, commentCharacter, classLabelsColumn);
                if(features.length != numColumns) {
                    throw new GATKException("filename does not encode a matrix, rows do not all have the same length");
                }
                rowsList.add(features);
            }
        } catch(IOException err) {
            throw new GATKException(err.getMessage());
        }
        return new Array2DRowRealMatrix(rowsList.toArray(new double[0][]), false);
    }

    private static double[] getNextCsvFeaturesRow(final BufferedReader reader, final String delimiter,
                                                  final String commentCharacter, final int classLabelsColumn) throws IOException {
        String line = reader.readLine();
        while(line.startsWith(commentCharacter) || line.isEmpty()) {
            line = reader.readLine();
        }
        final String[] words = line.split(delimiter, -1);
        final double[] features = Arrays.stream(words).mapToDouble(Double::valueOf).toArray();
        if(classLabelsColumn != CLASS_LABEL_COLUMN) {
            final double temp = features[classLabelsColumn];
            features[classLabelsColumn] = features[CLASS_LABEL_COLUMN];
            features[CLASS_LABEL_COLUMN] = temp;
        }
        return features;
    }

    public static RealMatrix sliceRows(final RealMatrix matrix, final int[] sliceRows) {
        final int[] allColumns = getRange(matrix.getColumnDimension());
        return matrix.getSubMatrix(sliceRows, allColumns);
    }


    public static int[] getClassLabels(final RealMatrix matrix) {
        return getClassLabels(matrix, CLASS_LABEL_COLUMN);
    }

    public static int[] getClassLabels(final RealMatrix matrix, final int classLabelColumn) {
        final int numRows = matrix.getRowDimension();
        final int[] classLabels = new int[numRows];
        for(int row = 0; row < numRows; ++row) {
            classLabels[row] = (int)Math.round(matrix.getEntry(row, classLabelColumn));
        }
        return classLabels;
    }

    public static abstract class GATKClassifier implements Serializable, KryoSerializable {
        private final double[][] singlePredictWrapper = new double [1][];

        // train classifier. Note, function should return "this"
        public abstract GATKClassifier train(final Map<String, Object> classifierParameters,
                                             final RealMatrix trainingMatrix);

        public abstract float[] trainAndReturnQualityTrace(
                final Map<String, Object> classifierParameters, final RealMatrix trainingMatrix,
                final RealMatrix evaluationMatrix, final int maxTrainingRounds, final int earlyStoppingRounds,
                final boolean maximizeEvalMetric);

        public abstract float[][] predictProbability(final RealMatrix matrix);

        public float[] predictProbability(final double[] featureVector) {
            singlePredictWrapper[0] = featureVector;
            final RealMatrix matrix = new Array2DRowRealMatrix(singlePredictWrapper, false);
            return (predictProbability(matrix)[0]);
        }

        public int[] predictClassLabels(final RealMatrix matrix) {
            final float [][] predictedProbabilities = predictProbability(matrix);
            final int [] predictedLabels = new int[predictedProbabilities.length];
            if(predictedProbabilities.length == 0) {
                return predictedLabels;
            }
            final int numColumns = predictedProbabilities[0].length;
            if(numColumns == 1) {
                // binary classifier, reporting only probability of class == 1 (or "true")
                for(int row = 0; row < predictedProbabilities.length; ++row) {
                    predictedLabels[row] = predictedProbabilities[row][0] >= 0.5 ? 1 : 0;
                }

            } else {
                // multiclass classifier (or at binary independently reporting probability of class 0 or 1)
                for (int row = 0; row < predictedProbabilities.length; ++row) {
                    predictedLabels[row] = argmax(predictedProbabilities[row]);
                }
            }
            return predictedLabels;
        }

        public void save(final String saveFilePath) throws IOException {
            try(FileOutputStream fileOutputStream = new FileOutputStream(saveFilePath)) {
                save(fileOutputStream);
            }
        }

        public void save(final FileOutputStream fileOutputStream) {
            final Kryo kryo = new Kryo();
            final Output output = new Output(fileOutputStream);
            kryo.writeClassAndObject(output, this);
            output.close();
        }

        public static GATKClassifier load(final String saveFilePath) throws IOException {
            try(FileInputStream fileInputStream = new FileInputStream(saveFilePath)) {
                return load(fileInputStream);
            }
        }

        public static GATKClassifier load(final FileInputStream fileInputStream) {
            final Kryo kryo = new Kryo();
            final Input input = new Input(fileInputStream);
            final GATKClassifier classifier = (GATKClassifier)kryo.readClassAndObject(input);
            input.close();
            return classifier;
        }

        public void chooseNumThreads(final Map<String, Object> classifierParameters, final String numThreadsKey,
                                     final RealMatrix trainingMatrix) {
            final int numCalibrationRows = NUM_CALIBRATION_CLASS_ROWS * 2;
            localLogger.info("numCalibrationRows = " + numCalibrationRows);
            localLogger.info("numTrainingRows = " + trainingMatrix.getRowDimension());
            final RealMatrix calibrationMatrix;
            if(trainingMatrix.getRowDimension() <= numCalibrationRows) {
                calibrationMatrix = trainingMatrix;
            } else {
                final int[] stratify = getClassLabels(trainingMatrix);
                localLogger.info("trainingFraction = " + numCalibrationRows / (float)trainingMatrix.getRowDimension());
                final TrainTestSplit trainTestSplit = TrainTestSplit.getTrainTestSplit(trainingMatrix,
                        numCalibrationRows / (float)trainingMatrix.getRowDimension(), new Random(), stratify);
                localLogger.info("numTrain = " + trainTestSplit.trainRows.length);
                localLogger.info("numTest = " + trainTestSplit.testRows.length);
                calibrationMatrix = sliceRows(trainingMatrix, trainTestSplit.trainRows);
            }
            final Map<String, Object> calibrationParams = new HashMap<>(classifierParameters);
            calibrationParams.put(NUM_TRAINING_ROUNDS_KEY, NUM_CALIBRATION_TRAINING_ROUNDS);
            final int maxNumThreads =
                    classifierParameters.containsKey(numThreadsKey) && (int)classifierParameters.get(numThreadsKey) > 0 ?
                            (int)classifierParameters.get(numThreadsKey) : Runtime.getRuntime().availableProcessors();
            if(maxNumThreads == 1) {
                classifierParameters.put(numThreadsKey, 1);
                return;
            }
            long bestElapsedTime = Long.MAX_VALUE;
            for(int numThreads = 1; numThreads < maxNumThreads; ++numThreads) {
                calibrationParams.put(numThreadsKey, numThreads);
                final long elapsedTime = getTrainingTime(calibrationParams, calibrationMatrix);
                if(elapsedTime < bestElapsedTime) {
                    bestElapsedTime = elapsedTime;
                    classifierParameters.put(numThreadsKey, numThreads);
                }
            }

        }

        private long getTrainingTime(final Map<String, Object> classifierParameters, final RealMatrix trainingMatrix) {
            final long startTime = System.nanoTime();
            train(classifierParameters, trainingMatrix);
            return System.nanoTime() - startTime;
        }

        private float[][] getCrossvalidatedTrainingTraces(final Map<String, Object> classifierParameters,
                                                          final RealMatrix trainMatrix, final List<TrainTestSplit> splits,
                                                          final int maxTrainingRounds, final int earlyStoppingRounds,
                                                          final boolean maximizeEvalMetric) {
            final int numCrossvalidationFolds = splits.size();
            float[][] trainingTraces = new float[numCrossvalidationFolds][];
            for(int fold = 0; fold < numCrossvalidationFolds; ++fold) {
                final TrainTestSplit split = splits.get(fold);
                trainingTraces[fold] = trainAndReturnQualityTrace(
                        classifierParameters, sliceRows(trainMatrix, split.trainRows),
                        sliceRows(trainMatrix, split.testRows), maxTrainingRounds, earlyStoppingRounds,
                        maximizeEvalMetric
                );
            }
            return trainingTraces;
        }

        private double getTrainingScore(final Map<String, Object> classifierParameters,
                                        final float[][] trainingTraces, final boolean maximizeEvalMetric) {
            // find the index with the best total (i.e. mean) score across rounds. This yields the optimal number of
            // rounds of training
            float bestTotalScore = maximizeEvalMetric ? Float.MIN_VALUE : Float.MAX_VALUE;
            int bestRoundIndex = -1;
            final int maxTrainingRounds = trainingTraces[0].length;
            for(int roundIndex = 0; roundIndex < maxTrainingRounds; ++roundIndex) {
                float roundScore = trainingTraces[0][roundIndex];
                for(int fold = 1; fold < trainingTraces.length; ++fold) {
                    roundScore += trainingTraces[fold][roundIndex];
                }
                if(maximizeEvalMetric ? roundScore > bestTotalScore : roundScore < bestTotalScore) {
                    // the score at this round of training is the best so far
                    bestTotalScore = roundScore;
                    bestRoundIndex = roundIndex;
                }
            }
            final int numTrainingRounds = bestRoundIndex + 1;

            // report the overall score for this set of parameters as the score of the *worst* trace at the optimal training
            // index. Selecting the worst trace demands high reliability from the classifier across similar data sets.
            float trainingScore = trainingTraces[0][bestRoundIndex];
            if(maximizeEvalMetric) {
                for (int fold = 1; fold < trainingTraces.length; ++fold) {
                    trainingScore = Math.min(trainingScore, trainingTraces[fold][bestRoundIndex]);
                }
            } else {
                for (int fold = 1; fold < trainingTraces.length; ++fold) {
                    trainingScore = Math.max(trainingScore, trainingTraces[fold][bestRoundIndex]);
                }
            }

            classifierParameters.put(NUM_TRAINING_ROUNDS_KEY, numTrainingRounds);
            return trainingScore;
        }

        public Map<String, Object> tuneClassifierParameters(final Map<String, Object> classifierParameters,
                                                            final Map<String, ClassifierParamRange<?>> tuneClassifierParameters,
                                                            final ClassifierTuningStrategy classifierTuningStrategy,
                                                            final RealMatrix trainMatrix, final Random random,
                                                            final int[] stratify, final int numCrossvalidationFolds,
                                                            final int maxTrainingRounds, final int earlyStoppingRounds,
                                                            final int numTuningRounds, final boolean maximizeEvalMetric) {
            final List<TrainTestSplit> splits = new ArrayList<>(numCrossvalidationFolds);
            TrainTestSplit.getCrossvalidationSplits(trainMatrix, numCrossvalidationFolds, random, stratify)
                    .forEachRemaining(splits::add);

            final ClassifierTuner classifierTuner;
            switch(classifierTuningStrategy) {
                case RANDOM:
                    classifierTuner = new RandomClassifierTuner(
                            this, tuneClassifierParameters, numTuningRounds, maximizeEvalMetric, random
                    );
                    break;
                default:
                    throw new IllegalStateException("Invalid ClassifierTuningStrategy: " + classifierTuningStrategy);
            }

            return classifierTuner.getBestParameters(classifierParameters, trainMatrix, splits, maxTrainingRounds, earlyStoppingRounds);

        }

        public int[] crossvalidatePredict(final RealMatrix dataMatrix, final Map<String, Object> classifierParameters,
                                          final Random random, final int[] stratify, final int numCrossvalidationFolds) {
            final int[] predictedLabels = new int[dataMatrix.getRowDimension()];
            final Iterator<TrainTestSplit> splitIterator = TrainTestSplit.getCrossvalidationSplits(
                    dataMatrix, numCrossvalidationFolds, random, stratify
            );
            while(splitIterator.hasNext()) {
                final TrainTestSplit split = splitIterator.next();
                // train on training data from this crossvalidation split
                train(classifierParameters, sliceRows(dataMatrix, split.trainRows));
                // predict on testing data from this split
                final int[] predictedTestLabels = predictClassLabels(sliceRows(dataMatrix, split.testRows));
                // and assign those values into the final predictions
                sliceAssign(predictedLabels, split.testRows, predictedTestLabels);
            }
            return predictedLabels;
        }
    }

    public static class TrainTestSplit {
        final int[] trainRows;
        final int[] testRows;

        TrainTestSplit(final int[] trainRows, final int[] testRows) {
            this.trainRows = trainRows;
            this.testRows = testRows;
        }

        public static TrainTestSplit getTrainTestSplit(final RealMatrix trainMatrix, final double trainingFraction,
                                                       final Random random, final int[] stratify) {
            final int numRows = trainMatrix.getRowDimension();
            if(stratify != null && stratify.length != numRows) {
                throw new IllegalArgumentException("stratify.length (" + stratify.length + ") not equal to num rows ("
                        + numRows + ")");
            }
            final int numTrain = (int)Math.round(numRows * trainingFraction);
            final int numTest = numRows - numTrain;
            if(numTrain <= 0) {
                if(trainingFraction < 0) {
                    throw new IllegalArgumentException(
                            "trainingFraction (" + trainingFraction + ") must be in range [0, 1]"
                    );
                }
                return new TrainTestSplit(new int[0], getRange(numRows));
            } else if(numTest <= 0) {
                if(trainingFraction > 1) {
                    throw new IllegalArgumentException(
                            "trainingFraction (" + trainingFraction + ") must be in range[0, 1]"
                    );
                }
                return new TrainTestSplit(getRange(numRows), new int[0]);
            }
            final int[] split_index_ordering = stratify == null ?
                    getRandomPermutation(random, numRows)
                    : TrainTestSplit.getStratfiedIndexOrdering(random, stratify);
            final int[] trainRows = new int[numTrain];
            final int[] testRows = new int[numTest];
            int nextTrainInd = 1;
            int nextTestInd = 1;
            for(int i = 0; i < numRows; ++i) {
                if(numTrain * nextTestInd >= numTest * nextTrainInd) {
                    // training set gets next index
                    trainRows[nextTrainInd - 1] = split_index_ordering[i];
                    ++nextTrainInd;
                } else {
                    testRows[nextTestInd - 1] = split_index_ordering[i];
                    ++nextTestInd;
                }
            }
            Arrays.sort(trainRows);
            Arrays.sort(testRows);
            return new TrainTestSplit(trainRows, testRows);
        }

        public static Iterator<TrainTestSplit> getCrossvalidationSplits(final RealMatrix trainMatrix,
                                                                 final int numCrossvalidationFolds,
                                                                 final Random random, final int[] stratify) {
            if(numCrossvalidationFolds < 2) {
                throw new IllegalArgumentException("numCrossvalidationFolds (" + numCrossvalidationFolds + ") must be >= 2");
            }
            final int numRows = stratify == null ? trainMatrix.getRowDimension() : stratify.length;
            final int[] split_index_ordering = stratify == null ?
                    getRandomPermutation(random, numRows)
                    : TrainTestSplit.getStratfiedIndexOrdering(random, stratify);

            return new FoldSplitIterator(split_index_ordering, numCrossvalidationFolds);
        }

        private static int[] getStratfiedIndexOrdering(final Random random, final int[] stratify) {
            /*
            logical (but memory inefficient) process
            1. make a random permutation, and use it to permute stratify
            final int[] permutation = getRange(stratify.length);
            final int[] permuted_stratify = slice(stratify, permutation);
            2. find the indices that would sort the permuted stratify array. The permutation ensures that entries with
               the same value in stratify will be in random order. However, they point to the wrong stratify values.
            final int[] permuted_sort_indices = argsort(permuted_stratify);
            3. unpermute the sort_indices to permute to the correct stratify values. Because the sort visited the values
               in random order, the ordering of stratify indices will be permuted (between indices pointing to equal
               stratify values).
            final int[] stratify_inds = slice(permutation, permuted_sort_indices);
            */
            final int[] permutation = getRandomPermutation(random, stratify.length);
            return slice(permutation, argsort(slice(stratify, permutation)));
        }

        private static class FoldSplitIterator implements Iterator<TrainTestSplit> {
            private final int[] split_index_ordering;
            private final int numFolds;
            private int fold;

            FoldSplitIterator(final int[] split_index_ordering, final int numFolds) {
                this.split_index_ordering = split_index_ordering;
                this.numFolds = numFolds;
                this.fold = 0;
            }

            @Override
            public boolean hasNext() {
                return fold < numFolds;
            }

            @Override
            public TrainTestSplit next() {
                final int numTest = (split_index_ordering.length - 1 - fold) / numFolds;
                final int numTrain = split_index_ordering.length - numTest;
                int[] testRows = new int[numTest];
                int[] trainRows = new int[numTrain];
                int trainIndex;
                int orderingIndex;
                if(fold > 0) {
                    for(trainIndex = 0; trainIndex < fold; ++trainIndex) {
                        trainRows[trainIndex] = split_index_ordering[trainIndex];
                    }
                    orderingIndex = trainIndex;
                } else {
                    orderingIndex = 0;
                    trainIndex = 0;
                }
                for(int testIndex = 0; testIndex < testRows.length; ++testIndex) {
                    testRows[testIndex] = split_index_ordering[orderingIndex];
                    final int orderingStop = Math.min(orderingIndex + numFolds, split_index_ordering.length);
                    for(++orderingIndex; orderingIndex < orderingStop; ++orderingIndex, ++trainIndex) {
                        trainRows[trainIndex] = split_index_ordering[orderingIndex];
                    }
                }

                ++fold;
                Arrays.sort(trainRows);
                Arrays.sort(testRows);
                return new TrainTestSplit(trainRows, testRows);
            }
        }
    }

    // To-do: write ClassifierTuner with strategy = BAYES
    public enum ClassifierTuningStrategy { RANDOM }

    static abstract class ClassifierTuner {
        private final GATKClassifier classifier;
        protected final Map<String, ClassifierParamRange<?>> tuneParameters;
        protected final int numTuningRounds;
        protected final boolean maximizeEvalMetric;
        protected final List<Map<String, Object>> hyperparameterSets;
        protected final List<Double> hyperparameterScores;

        ClassifierTuner(final GATKClassifier classifier, final Map<String, ClassifierParamRange<?>> tuneParameters,
                        final int numTuningRounds, final boolean maximizeEvalMetric) {
            if(numTuningRounds < 1) {
                throw new IllegalArgumentException("numTuningRounds (" + numTuningRounds + ") must be >= 1");
            }
            this.classifier = classifier;
            this.tuneParameters = tuneParameters;
            this.numTuningRounds = numTuningRounds;
            this.hyperparameterSets = new ArrayList<> ();
            this.hyperparameterScores = new ArrayList<>();
            this.maximizeEvalMetric = maximizeEvalMetric;
        }

        abstract protected Map<String, Object> chooseNextHyperparameters();

        Map<String, Object> getBestParameters(final Map<String, Object> classifierParameters,
                                              final RealMatrix trainMatrix, final List<TrainTestSplit> splits,
                                              final int maxTrainingRounds, final int earlyStoppingRounds) {
            Map<String, Object> bestParameters = null;
            double bestScore = maximizeEvalMetric ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
            for(int i = 0; i < numTuningRounds; ++i) {
                final Map<String, Object> hyperparameters = chooseNextHyperparameters();
                hyperparameterSets.add(hyperparameters);
                final Map<String, Object> testParameters = new HashMap<>(classifierParameters);
                testParameters.putAll(hyperparameters);
                final float[][] trainingTraces = classifier.getCrossvalidatedTrainingTraces(
                        testParameters, trainMatrix, splits, maxTrainingRounds, earlyStoppingRounds,
                        maximizeEvalMetric
                );
                final double score = classifier.getTrainingScore(testParameters, trainingTraces, maximizeEvalMetric);
                hyperparameterScores.add(score);
                // This is the new best score if a) it is better than the previous best score so far -OR-
                //                               b) it exactly ties the best score, but uses fewer rounds of training
                if((maximizeEvalMetric ? score > bestScore : score < bestScore)
                        || (score == bestScore &&
                            (int)testParameters.get(NUM_TRAINING_ROUNDS_KEY) < (int)bestParameters.get(NUM_TRAINING_ROUNDS_KEY))) {
                    bestScore = score;
                    bestParameters = testParameters;
                }
            }
            return bestParameters;
        }
    }

    static class RandomClassifierTuner extends ClassifierTuner {
        private final Map<String, Object[]> randomParameters;
        RandomClassifierTuner(final GATKClassifier classifier, final Map<String, ClassifierParamRange<?>> tuneParameters,
                              final int numTuningRounds, final boolean maximizeEvalMetric, final Random random) {
            super(classifier, tuneParameters, numTuningRounds, maximizeEvalMetric);
            randomParameters = tuneParameters.entrySet().stream().collect(
                    Collectors.toMap(
                            Map.Entry::getKey,
                            x -> x.getValue().getRandomSamples(random, numTuningRounds)
                    )
            );
        }

        @Override
        protected Map<String, Object> chooseNextHyperparameters() {
            final int index = hyperparameterScores.size();
            return randomParameters.entrySet().stream().collect(
                    Collectors.toMap(
                            Map.Entry::getKey,
                            x -> x.getValue()[index]
                    )
            );
        }

    }


    public interface ClassifierParamRange<T> {
        T[] getRandomSamples(final Random random, final int numSamples);
    }

    public static class ClassifierLinearParamRange implements ClassifierParamRange<Double> {
        private final double low;
        private final double high;

        public ClassifierLinearParamRange(double low, double high) {
            this.low = low;
            this.high = high;
        }

        public Double[] getRandomSamples(final Random random, final int numSamples) {
            final Double[] samples = new Double[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (high + low) / 2.0;
                }
                return samples;
            }
            final double delta = (high - low) / (numSamples - 1);
            double val = low;
            final int[] permutation = getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val += delta;
                samples[permutation[i]] = val;
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    public static class ClassifierLogParamRange implements ClassifierParamRange<Double> {
        private final double low;
        private final double high;

        public ClassifierLogParamRange(double low, double high) {
            if(low * high <= 0) {
                throw new IllegalArgumentException("low (" + low + ") and high (" + high + ") must be the same sign, and non-zero");
            }
            this.low = low;
            this.high = high;
        }

        public Double[] getRandomSamples(final Random random, final int numSamples) {
            final Double[] samples = new Double[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = Math.sqrt(high * low);
                }
                return samples;
            }
            final double delta = Math.pow(high / low, 1.0 / (numSamples - 1));
            double val = low;
            final int[] permutation = getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val *= delta;
                samples[permutation[i]] = val;
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    public static class ClassifierIntegerLinearParamRange implements ClassifierParamRange<Integer> {
        private final int low;
        private final int high;

        public ClassifierIntegerLinearParamRange(int low, int high) {
            this.low = low;
            this.high = high;
        }

        public Integer[] getRandomSamples(final Random random, final int numSamples) {
            final Integer[] samples = new Integer[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (int)Math.round((high + low) / 2.0);
                }
                return samples;
            }
            final double delta = (high - low) / (double)(numSamples - 1);
            double val = low;
            final int[] permutation = getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val += delta;
                samples[permutation[i]] = (int)Math.round(val);
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    public static class ClassifierIntegerLogParamRange implements ClassifierParamRange<Integer> {
        private final int low;
        private final int high;

        public ClassifierIntegerLogParamRange(int low, int high) {
            if(low * high <= 0) {
                throw new IllegalArgumentException("low (" + low + ") and high (" + high + ") must be the same sign, and non-zero");
            }
            this.low = low;
            this.high = high;
        }

        public Integer[] getRandomSamples(final Random random, final int numSamples) {
            final Integer[] samples = new Integer[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (int)Math.round(Math.sqrt((double)high * low));
                }
                return samples;
            }
            final double delta = Math.pow(high / (double)low, 1.0 / (numSamples - 1));
            double val = low;
            final int[] permutation = getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val *= delta;
                samples[permutation[i]] = (int)Math.round(val);
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    public static Integer[] getRange(final Integer numElements) {
        if(numElements < 0) {
            throw new IllegalArgumentException("numElements must be >= 0");
        }
        final Integer[] range = new Integer[numElements];
        for(Integer i = 0; i < numElements; ++i) {
            range[i] = i;
        }
        return range;
    }

    public static int[] getRange(final int numElements) {
        if(numElements < 0) {
            throw new IllegalArgumentException("numElements must be >= 0");
        }
        return IntStream.range(0, numElements).toArray();
    }

    public static int[] slice(final int[] arr, final int[] indices) {
        final int[] sliced_arr = new int[indices.length];
        for(int i = 0; i < indices.length; ++i) {
            sliced_arr[i] = arr[indices[i]];
        }
        return sliced_arr;
    }

    public static void sliceAssign(final int[] arr, final int[] indices, final int[] newValues) {
        if(indices.length != newValues.length) {
            throw new IllegalArgumentException("length of indices does not match length of newValues");
        }
        for(int i = 0; i < indices.length; ++i) {
            arr[indices[i]] = newValues[i];
        }
    }

    private static int argmax(final float[] arr) {
        int bestIndex = -1;
        float bestValue = Float.MIN_VALUE;
        for(int index = 0; index < arr.length; ++index) {
            if(arr[index] > bestValue) {
                bestIndex = index;
                bestValue = arr[index];
            }
        }
        return bestIndex;
    }

    public static int[] argsort(final int[] arr) {
        final Integer[] sortIndices = getRange((Integer)arr.length);
        Arrays.sort(sortIndices, Comparator.comparingInt(ind -> arr[ind]));
        return ArrayUtils.toPrimitive(sortIndices);
    }

    private static <T> int[] uniqueLabels(final T[] arr) {
        int[] labels = new int[arr.length];
        Map<T, Integer> values = new HashMap<>();
        for(int i = 0; i < arr.length; ++i) {
            final T value = arr[i];
            if(values.containsKey(value)) {
                labels[i] = values.get(value);
            } else {
                labels[i] = values.size();
                values.put(value, values.size());
            }
        }
        return labels;
    }

    private static int[] uniqueLabels(final float[] arr) {
        int[] labels = new int[arr.length];
        Map<Float, Integer> values = new HashMap<>();
        for(int i = 0; i < arr.length; ++i) {
            final float value = arr[i];
            if(values.containsKey(value)) {
                labels[i] = values.get(value);
            } else {
                labels[i] = values.size();
                values.put(value, values.size());
            }
        }
        return labels;
    }

    public static int[] getRandomPermutation(final Random random, final int numElements) {
        // Knuth shuffle
        final int[] permutation = getRange(numElements);
        for(int i = numElements - 1; i > 0; --i) {
            final int swap_ind = random.nextInt(i);
            final int swap_val = permutation[swap_ind];
            permutation[swap_ind] = permutation[i];
            permutation[i] = swap_val;
        }
        return permutation;
    }

    public static double getPredictionAccuracy(final int[] predictedLabels, final int[] correctLabels) {
        int numCorrect = 0;
        for(int row = 0; row < correctLabels.length; ++row) {
            if(predictedLabels[row] == correctLabels[row]) {
                ++numCorrect;
            }
        }
        return numCorrect / (double)correctLabels.length;
    }
}

