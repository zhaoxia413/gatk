package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.KryoSerializable;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import scala.Int;

import javax.crypto.Mac;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
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
        private static final long serialVersionUID = 1L;

        private final double[][] singlePredictWrapper = new double [1][];

        // train classifier. Note, function should return "this"
        public abstract GATKClassifier train(final Map<String, Object> classifierParameters,
                                             final RealMatrix trainingMatrix);

        public abstract double[] trainAndReturnQualityTrace(
                final Map<String, Object> classifierParameters, final RealMatrix trainingMatrix,
                final RealMatrix evaluationMatrix, final int maxTrainingRounds, final int earlyStoppingRounds);

        public abstract double[][] predictProbability(final RealMatrix matrix);

        public abstract boolean getMaximizeEvalMetric(final Map<String, Object> classifierParameters);

        /*
        public defined methods
         */
        public double[] predictProbability(final double[] featureVector) {
            singlePredictWrapper[0] = featureVector;
            final RealMatrix matrix = new Array2DRowRealMatrix(singlePredictWrapper, false);
            return (predictProbability(matrix)[0]);
        }

        public int[] predictClassLabels(final RealMatrix matrix) {
            final double [][] predictedProbabilities = predictProbability(matrix);
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

        public Map<String, Object> tuneClassifierParameters(final Map<String, Object> classifierParameters,
                                                            final Map<String, ClassifierParamRange<?>> tuneClassifierParameters,
                                                            final ClassifierTuningStrategy classifierTuningStrategy,
                                                            final RealMatrix trainMatrix, final Random random,
                                                            final int[] stratify, final int numCrossvalidationFolds,
                                                            final int maxTrainingRounds, final int earlyStoppingRounds,
                                                            final int numTuningRounds) {
            final List<TrainTestSplit> splits = new ArrayList<>(numCrossvalidationFolds);
            TrainTestSplit.getCrossvalidationSplits(
                    numCrossvalidationFolds, trainMatrix.getRowDimension(), random, stratify
            ).forEachRemaining(splits::add);

            final ClassifierTuner classifierTuner;
            switch(classifierTuningStrategy) {
                case RANDOM:
                    classifierTuner = new RandomClassifierTuner(
                            this, classifierParameters, tuneClassifierParameters, numTuningRounds, random
                    );
                    break;
                default:
                    throw new IllegalStateException("Invalid ClassifierTuningStrategy: " + classifierTuningStrategy);
            }

            return classifierTuner.getBestParameters(trainMatrix, splits, maxTrainingRounds, earlyStoppingRounds);

        }

        public int[] crossvalidatePredict(final RealMatrix dataMatrix, final Map<String, Object> classifierParameters,
                                          final Random random, final int[] stratify, final int numCrossvalidationFolds) {
            final int[] predictedLabels = new int[dataMatrix.getRowDimension()];
            final Iterator<TrainTestSplit> splitIterator = TrainTestSplit.getCrossvalidationSplits(
                    numCrossvalidationFolds, dataMatrix.getRowDimension(), random, stratify
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
                localLogger.info("trainingFraction = " + numCalibrationRows / (double)trainingMatrix.getRowDimension());
                final TrainTestSplit trainTestSplit = TrainTestSplit.getTrainTestSplit(
                        numCalibrationRows / (double)trainingMatrix.getRowDimension(),
                        trainingMatrix.getRowDimension(), new Random(), stratify
                );
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

        /*
        private methods
         */
        private long getTrainingTime(final Map<String, Object> classifierParameters, final RealMatrix trainingMatrix) {
            final long startTime = System.nanoTime();
            train(classifierParameters, trainingMatrix);
            return System.nanoTime() - startTime;
        }

        private double[][] getCrossvalidatedTrainingTraces(final Map<String, Object> classifierParameters,
                                                          final RealMatrix trainMatrix, final List<TrainTestSplit> splits,
                                                          final int maxTrainingRounds, final int earlyStoppingRounds) {
            final int numCrossvalidationFolds = splits.size();
            double[][] trainingTraces = new double[numCrossvalidationFolds][];
            for(int fold = 0; fold < numCrossvalidationFolds; ++fold) {
                final TrainTestSplit split = splits.get(fold);
                trainingTraces[fold] = trainAndReturnQualityTrace(
                        classifierParameters, sliceRows(trainMatrix, split.trainRows),
                        sliceRows(trainMatrix, split.testRows), maxTrainingRounds, earlyStoppingRounds
                );
            }
            return trainingTraces;
        }

        private double getTrainingScore(final Map<String, Object> classifierParameters,
                                        final double[][] trainingTraces) {
            // find the index with the best total (i.e. mean) score across rounds. This yields the optimal number of
            // rounds of training
            final boolean maximizeEvalMetric = getMaximizeEvalMetric(classifierParameters);
            double bestTotalScore = maximizeEvalMetric ? Double.MIN_VALUE : Double.MAX_VALUE;
            int bestRoundIndex = -1;
            final int maxTrainingRounds = trainingTraces[0].length;
            for(int roundIndex = 0; roundIndex < maxTrainingRounds; ++roundIndex) {
                double roundScore = trainingTraces[0][roundIndex];
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
            double trainingScore = trainingTraces[0][bestRoundIndex];
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
    }

    public static class TrainTestSplit {
        public final int[] trainRows;
        public final int[] testRows;

        TrainTestSplit(final int[] trainRows, final int[] testRows) {
            this.trainRows = trainRows;
            this.testRows = testRows;
        }


        public static TrainTestSplit getTrainTestSplit(final double trainingFraction, final int numRows,
                                                       final Random random, final int[] stratify) {
            if(stratify != null && stratify.length != numRows) {
                throw new IllegalArgumentException(
                        "stratify.length (" + stratify.length + ") != numRows (" + numRows + ")"
                );
            }

            final long numTrain = Math.round(numRows * trainingFraction);
            final long numTest = numRows - numTrain;
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
            final int[] trainRows = new int[(int)numTrain];
            final int[] testRows = new int[(int)numTest];
            int nextTrainInd = 1;
            int nextTestInd = 1;
            for(final int split_index : split_index_ordering) {
                if(numTrain * nextTestInd >= numTest * nextTrainInd) {
                    // training set gets next index
                    trainRows[nextTrainInd - 1] = split_index;
                    ++nextTrainInd;
                } else {
                    testRows[nextTestInd - 1] = split_index;
                    ++nextTestInd;
                }
            }
            Arrays.sort(trainRows);
            Arrays.sort(testRows);
            return new TrainTestSplit(trainRows, testRows);
        }

        public static Iterator<TrainTestSplit> getCrossvalidationSplits(final int numCrossvalidationFolds, final int numRows,
                                                                        final Random random, final int[] stratify) {
            if(numCrossvalidationFolds < 2) {
                throw new IllegalArgumentException("numCrossvalidationFolds (" + numCrossvalidationFolds + ") must be >= 2");
            }
            if(stratify != null && stratify.length != numRows) {
                throw new IllegalArgumentException(
                        "stratify.length (" + stratify.length + ") != numRows (" + numRows + ")"
                );
            }
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
                final int numTest = 1 + (split_index_ordering.length - 1 - fold) / numFolds;
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
        protected final Map<String, Object> classifierParameters;
        protected final Map<String, ClassifierParamRange<?>> tuneParameters;
        protected final int numTuningRounds;
        protected final boolean maximizeEvalMetric;
        protected final List<Map<String, Object>> hyperparameterSets;
        protected final List<Double> hyperparameterScores;

        ClassifierTuner(final GATKClassifier classifier, final Map<String, Object> classifierParameters,
                        final Map<String, ClassifierParamRange<?>> tuneParameters,
                        final int numTuningRounds) {
            if(numTuningRounds < 1) {
                throw new IllegalArgumentException("numTuningRounds (" + numTuningRounds + ") must be >= 1");
            }
            this.classifier = classifier;
            this.tuneParameters = tuneParameters;
            this.numTuningRounds = numTuningRounds;
            this.hyperparameterSets = new ArrayList<> ();
            this.hyperparameterScores = new ArrayList<>();
            this.classifierParameters = classifierParameters;
            this.maximizeEvalMetric = classifier.getMaximizeEvalMetric(classifierParameters);
        }

        abstract protected Map<String, Object> chooseNextHyperparameters();

        Map<String, Object> getBestParameters(final RealMatrix trainMatrix, final List<TrainTestSplit> splits,
                                              final int maxTrainingRounds, final int earlyStoppingRounds) {
            Map<String, Object> bestParameters = null;
            double bestScore = maximizeEvalMetric ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;

            MachineLearningUtils.localLogger.info("Getting best parameters");
            ConsoleProgressBar progress = new ConsoleProgressBar(numTuningRounds);
            for(int i = 0; i < numTuningRounds; ++i) {
                final Map<String, Object> hyperparameters = chooseNextHyperparameters();
                hyperparameterSets.add(hyperparameters);
                final Map<String, Object> testParameters = new HashMap<>(classifierParameters);
                testParameters.putAll(hyperparameters);
                final double[][] trainingTraces = classifier.getCrossvalidatedTrainingTraces(
                        testParameters, trainMatrix, splits, maxTrainingRounds, earlyStoppingRounds
                );
                final double score = classifier.getTrainingScore(testParameters, trainingTraces);
                hyperparameterScores.add(score);
                // This is the new best score if a) it is better than the previous best score so far -OR-
                //                               b) it exactly ties the best score, but uses fewer rounds of training
                if((maximizeEvalMetric ? score > bestScore : score < bestScore)
                        || (score == bestScore &&
                            (int)testParameters.get(NUM_TRAINING_ROUNDS_KEY) < (int)bestParameters.get(NUM_TRAINING_ROUNDS_KEY))) {
                    bestScore = score;
                    bestParameters = testParameters;
                }
                progress.update(1);
            }
            return bestParameters;
        }
    }

    static class RandomClassifierTuner extends ClassifierTuner {
        private final Map<String, Object[]> randomParameters;
        RandomClassifierTuner(final GATKClassifier classifier, final Map<String, Object> classifierParameters,
                              final Map<String, ClassifierParamRange<?>> tuneParameters, final int numTuningRounds,
                              final Random random) {
            super(classifier, classifierParameters, tuneParameters, numTuningRounds);
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
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
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
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
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
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
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
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
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

    public static class ConsoleProgressBar {
        private static final long MIN_UPDATE_INTERVAL_NS = 500000000; // = 0.5 sec
        private static final int NUM_BAR_CHARACTERS = 10;

        private static final double SECONDS_IN_MINUTE = 60.0;
        private static final long MINUTES_IN_HOUR = 60;
        private static final long HOURS_IN_DAY = 24;
        private static final String CARRIAGE_RETURN = "\r";

        private final long workToDo;
        private final long bornTime;
        private long nextUpdateTime;
        private long workDone;
        private long workRemaining;
        private int maxOutputLength;
        private final String updateInfoFormat;

        ConsoleProgressBar(final long workToDo) {
            if(workToDo <= 0) {
                throw new IllegalArgumentException("workToDo must be > 0");
            }
            this.workToDo = workToDo;
            workDone = 0;
            workRemaining = workToDo;
            bornTime = System.nanoTime();
            nextUpdateTime = bornTime + MIN_UPDATE_INTERVAL_NS;
            maxOutputLength = 0;
            final int workToDoLength = String.format("%d", workToDo).length();
            final String workDoneFormat = String.format("%%%dd/%%%dd", workToDoLength, workToDoLength);
            updateInfoFormat = workDoneFormat + " %.1f%% elapsed %s, remaining %s";
            drawBar(0.0, Double.NaN);
        }

        public void update(final long workJustCompleted) {
            if(workJustCompleted <= 0) {
                throw new IllegalArgumentException("workJustCompleted must be > 0");
            }
            workRemaining -= workJustCompleted;
            workDone += workJustCompleted;
            final long now = System.nanoTime();
            if(now < nextUpdateTime && workRemaining > 0) {
                return;  // avoid thrashing to the screen
            } else {
                nextUpdateTime = now + MIN_UPDATE_INTERVAL_NS;
            }
            final double elapsedTimeSec = 1.0e-9 * (now - bornTime);
            final double workPerSec = workDone / elapsedTimeSec;
            final double remainingTimeSec = workRemaining / workPerSec;
            System.out.flush();
            drawBar(elapsedTimeSec, remainingTimeSec);
        }

        private void drawBar(final double elapsedTimeSec, final double remainingTimeSec) {
            // NOTE: carriage return means each output will start from the beginning of the line
            final StringBuilder stringBuilder = new StringBuilder(CARRIAGE_RETURN);
            // draw actual progress bar
            final int numBarFilled = (int)(NUM_BAR_CHARACTERS * workDone / workToDo);
            final int numBarUnfilled = NUM_BAR_CHARACTERS - numBarFilled;
            stringBuilder.append("|");
            stringBuilder.append(StringUtils.repeat('#', numBarFilled));
            stringBuilder.append(StringUtils.repeat(' ', numBarUnfilled));
            stringBuilder.append("| ");
            // write summary statistics on completion amount, times
            stringBuilder.append(
                    workDone > 0 ?
                    String.format(
                        updateInfoFormat, workDone, workToDo, workDone * 100.0 / workToDo,
                        secondsToTimeString(elapsedTimeSec), secondsToTimeString(remainingTimeSec)
                    )
                    : String.format(updateInfoFormat, 0, workToDo, 0.0, secondsToTimeString(0.0), "???")
            );
            // do any necessary padding
            if(stringBuilder.length() > maxOutputLength) {
                maxOutputLength = stringBuilder.length();
            } else if(stringBuilder.length() < maxOutputLength) {
                // pad with spaces to obliterate previous message
                stringBuilder.append(StringUtils.repeat(' ', maxOutputLength - stringBuilder.length()));
            }
            // if we're done, add newline
            if(workRemaining <= 0) {
                stringBuilder.append("\n");
            }
            // write out and flush
            System.out.print(stringBuilder.toString());
            System.out.flush();
        }

        private static String secondsToTimeString(double seconds) {
            if(seconds < SECONDS_IN_MINUTE) {
                return String.format("%.1fs", seconds);
            }
            long minutes = (int)Math.floor(seconds / SECONDS_IN_MINUTE);
            seconds = seconds % SECONDS_IN_MINUTE;
            long hours = minutes / MINUTES_IN_HOUR;
            if(hours <= 0) {
                return String.format("%dm %.1fs", minutes, seconds);
            }
            minutes = minutes % MINUTES_IN_HOUR;
            long days = hours / HOURS_IN_DAY;
            if(days <= 0) {
                return String.format("%dh %dm %.1fs", hours, minutes, seconds);
            } else {
                hours = hours % HOURS_IN_DAY;
                return String.format("%dd %dh %dm %.1fs", days, hours, minutes, seconds);
            }
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

    public static int argmax(final double[] arr) {
        int bestIndex = -1;
        double bestValue = Double.MIN_VALUE;
        for(int index = 0; index < arr.length; ++index) {
            if(arr[index] > bestValue) {
                bestIndex = index;
                bestValue = arr[index];
            }
        }
        return bestIndex;
    }

    public static RealMatrix concatenateColumns(final RealMatrix matrixA, final RealMatrix matrixB) {
        if(matrixA.getRowDimension() != matrixB.getRowDimension()) {
            throw new IllegalArgumentException("matrixA and matrixB do not have same number of rows.");
        }
        final RealMatrix matrixC = matrixA.createMatrix(
                matrixA.getRowDimension(),
                matrixA.getColumnDimension() + matrixB.getColumnDimension()
        );
        matrixC.setSubMatrix(matrixA.getData(), 0, 0);
        matrixC.setSubMatrix(matrixB.getData(), 0, matrixA.getColumnDimension());
        return matrixC;
    }

    public static int[] argsort(final int[] arr) {
        final Integer[] sortIndices = getRange((Integer)arr.length);
        Arrays.sort(sortIndices, Comparator.comparingInt(ind -> arr[ind]));
        return ArrayUtils.toPrimitive(sortIndices);
    }

    public static <T extends Comparable<? super T>> int[] argsort(final T[] arr) {
        final Integer[] sortIndices = getRange((Integer)arr.length);
        Arrays.sort(sortIndices, Comparator.comparing(ind -> arr[ind]));
        return ArrayUtils.toPrimitive(sortIndices);
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


    public static int[] stratifyMatrixToStratifyArray(final RealMatrix stratifyMatrix, final int numBins,
                                                      final int minCountsPerStratifyValue) {
        return stratifyMatrixToStratifyArray(stratifyMatrix, numBins, minCountsPerStratifyValue, new ArrayList<>());
    }

    public static int[] stratifyMatrixToStratifyArray(final RealMatrix stratifyMatrix, final int numBins,
                                                      final int minCountsPerStratifyValue,
                                                      final Collection<Integer> categoricalColumns) {
        // first form binned version of each column
        final int[][] binnedStratifyMatrix = new int[stratifyMatrix.getRowDimension()][stratifyMatrix.getColumnDimension()];
        for(int columnIndex = 0; columnIndex < stratifyMatrix.getColumnDimension(); ++columnIndex) {
            final int[] columnBins = categoricalColumns.contains(columnIndex) ?
                    getCategoryCodes(stratifyMatrix.getColumn(columnIndex))
                    : getBinnedColumn(stratifyMatrix.getColumn(columnIndex), numBins);
            for(int row = 0; row < stratifyMatrix.getRowDimension(); ++row) {
                binnedStratifyMatrix[row][columnIndex] = columnBins[row];
            }
        }

        // Get a set of unique rows
        Set<List<Integer>> resultsToCheck = getUniqueRows(
                binnedStratifyMatrix,
                IntStream.range(0, stratifyMatrix.getColumnDimension()).boxed().collect(Collectors.toList()),
                stratifyMatrix.getColumnDimension()
        );
        // While there are rows that have too few instances, decrease the number of columns under consideration (only
        // for those undersized rows)
        final Set<List<Integer>> uniqueResults = new HashSet<>();
        for(int useNumColumns = stratifyMatrix.getColumnDimension() - 1; useNumColumns > 0; --useNumColumns) {
            // Add undersized rows to a set to be reprocessed. Add properly sized rows to final uniqueResults set.
            final Set<List<Integer>> undersizedStratify = new HashSet<>();
            for(final List<Integer> uniqueResult : resultsToCheck) {
                if(uniqueResult.size() < minCountsPerStratifyValue) {
                    undersizedStratify.add(uniqueResult);
                } else {
                    uniqueResults.add(uniqueResult);
                }
            }
            if(undersizedStratify.isEmpty()) {
                break; // no problematic values, we're done
            } else if(useNumColumns == 1) {
                // can't consider fewer columns, lump remaining rows into the stratify value with the fewest members, to
                // make an "odd-ball" value.
                if(uniqueResults.isEmpty()) {
                    // handle edge case where every unique row has too few counts
                    uniqueResults.add(new ArrayList<>());
                }
                final List<Integer> smallest = uniqueResults.stream().min(Comparator.comparing(List::size))
                        .orElseThrow(NoSuchElementException::new);
                for(final List<Integer> uniqueResult : undersizedStratify) {
                    smallest.addAll(uniqueResult);
                }
                Collections.sort(smallest);
            } else {
                // find unique rows from subset that are too small, looking at one fewer column
                final List<Integer> tooSmall = undersizedStratify.stream().flatMap(Collection::stream).sorted()
                        .collect(Collectors.toList());
                resultsToCheck = getUniqueRows(binnedStratifyMatrix, tooSmall, useNumColumns);
            }
        }


        final int[] stratifyArray = new int[binnedStratifyMatrix.length];
        int stratifyValue = 0;
        for(final List<Integer> uniqueResult : uniqueResults) {
            for(final int index : uniqueResult) {
                stratifyArray[index] = stratifyValue;
            }
            ++stratifyValue;
        }
        return stratifyArray;
    }

    private static int[] getCategoryCodes(final double[] column) {
        final int[] categoryCodes = new int[column.length];
        final Map<Double, Integer> codesMap = new HashMap<>();
        for(int i = 0; i < column.length; ++i) {
            final double value = column[i];
            final Integer code = codesMap.getOrDefault(value, null);
            if(code == null) {
                categoryCodes[i] = codesMap.size();
                codesMap.put(value, codesMap.size());
            } else {
                categoryCodes[i] = code;
            }
        }
        return categoryCodes;
    }

    private static int[] getBinnedColumn(final double[] column, final int numBins) {
        final int numUniqueColumnValues = (int)DoubleStream.of(column).distinct().count();
        if(numUniqueColumnValues <= numBins) {
            // Too much repetition to bin the data into the requested number of bins, just return unique values.
            return getCategoryCodes(column);
        }

        // Attempt to get requested number of percentiles, insisting that all percentiles are unique. If some values are
        // repeated, "percentiles" will not be evenly spaced, and it may not return exactly the requested number.
        // NOTE: the last percentile will not be used for binning (percentiles are "posts", bins are "fence") so request
        // one more percentile than bins
        final double[] percentiles = getUniquePercentiles(column, numBins + 1);

        // bin column to specified percentiles, dumping NaN into the last bin if it is present. Note that the percentiles
        // will be properly sized to account for the presence of NaN values.
        // NOTE: binarySearch is limited to ignore last percentile, because we don't want the max value in the array
        // being mapped past the maximum requested number of bins.
        return DoubleStream.of(column).mapToInt(
                v -> Double.isNaN(v) ? numBins - 1 : Arrays.binarySearch(percentiles, 0, percentiles.length, v)
        ).toArray();
    }

    private static double[] getUniquePercentiles(final double[] column, final int numPercentiles) {
        final int numNaN = (int)Arrays.stream(column).filter(Double::isNaN).count();
        final int numSortablePercentiles = numNaN == 0 ? numPercentiles : numPercentiles - 1;
        final Percentile percentileEvaluator = new Percentile().withEstimationType(Percentile.EstimationType.R_1);
        percentileEvaluator.setData(column);
        double[] percentiles = percentileSpace(percentileEvaluator, numSortablePercentiles);
        if(percentiles.length == numSortablePercentiles) {
            return percentiles; // should be the case for typical non-repeating data
        }

        final int numSortable = column.length - numNaN;
        int numRequestHigh = numSortablePercentiles;
        while(percentiles.length < numSortablePercentiles && numRequestHigh < numSortable) {
            numRequestHigh = Math.min(2 * numRequestHigh, numSortable);
            percentiles = percentileSpace(percentileEvaluator, numRequestHigh);
        }
        if(percentiles.length == numSortablePercentiles) {
            return percentiles;
        }
        int numRequestLow = numSortablePercentiles;
        int range = numRequestHigh - numRequestLow;
        while(range > 1) {
            final int numRequest = numRequestLow + range / 2;
            percentiles = percentileSpace(percentileEvaluator, numRequest);
            if(percentiles.length < numSortablePercentiles) {
                numRequestLow = numRequest;
            } else if(percentiles.length > numSortablePercentiles){
                numRequestHigh = numRequest;
            } else {
                return percentiles;
            }
            range = numRequestHigh - numRequestLow;
        }

        // I'm not convinced that it's possible to get down to this point. It implies that you checked for one more
        // percentile but got 2 or more new unique values back. Just in case, insist on having *fewer* than requested,
        // since other approximate cases always yield that outcome.
        percentiles = percentileSpace(percentileEvaluator, numRequestLow);

        return percentiles;
    }

    private static final double[] percentileSpace(final Percentile percentileEvaluator, final int numPercentiles) {
        final double low = 50.0 / percentileEvaluator.getData().length;
        final double high = 100.0 - low;
        final double coef = (high - low) / (numPercentiles - 1);
        return IntStream.range(0, numPercentiles).mapToDouble(i -> percentileEvaluator.evaluate(low + i * coef))
                .distinct().toArray();
    }


    private static Set<List<Integer>> getUniqueRows(final int[][] binnedStratifyMatrix,
                                                    final List<Integer> checkRows,
                                                    final int useNumColumns) {
        final Map<String, List<Integer>> uniqueResultMap = new HashMap<>();

        for(final int rowIndex : checkRows) {
            final int[] row = binnedStratifyMatrix[rowIndex];
            final String rowString = Arrays.toString(Arrays.copyOfRange(row, 0, useNumColumns));
            final List<Integer> uniqueResult = uniqueResultMap.getOrDefault(rowString, null);
            if(uniqueResult == null) {
                final List<Integer> newResult = new ArrayList<>();
                newResult.add(rowIndex);
                uniqueResultMap.put(rowString, newResult);
            } else {
                uniqueResult.add(rowIndex);
            }
        }

        return new HashSet<>(uniqueResultMap.values());
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

