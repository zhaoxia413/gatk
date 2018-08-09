package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.ml.classification.Classifier;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

/**
 * (Internal) Trains classifier for filtering BreakpointEvidence, and produces supporting files for testing"
 *
 * <p>This tool trains a classifier used to filter BreakpointEvidence, discarding evidence that is unlikely to overlap
 * the break point of a structural variant, and passing the rest (which will be combined via overlap to form
 * assembly intervals).</p>
 * After training the classifier, the model file is saved, along with supporting files that are used in unit" +
 * or integration tests.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A SAM/BAM/CRAM file of paired-end, aligned and coordinate-sorted reads.</li>
 *     <li>A BWA index image for the reference.
 *         You can use BwaMemIndexImageCreator to create the index image file.</li>
 *     <li>A list of ubiquitous kmers to ignore.
 *         You can use FindBadGenomicGenomicKmersSpark to create the list of kmers to ignore.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A file of aligned contigs.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk FindBreakpointEvidenceSpark \
 *     -I input_reads.bam \
 *     --aligner-index-image reference.img \
 *     --kmers-to-ignore ignored_kmers.txt \
 *     -O assemblies.sam
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 *
 * <h3>Caveats</h3>
 * <p>Expected input is a paired-end, coordinate-sorted BAM with around 30x coverage.
 * Coverage much lower than that probably won't work well.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Trains classifier for filtering BreakpointEvidence, and produces supporting files for testing",
        summary =
                "This tool trains a classifier used to filter BreakpointEvidence, discarding evidence that is unlikely to overlap" +
                        " the break point of a structural variant, and passing the rest (which will be combined via overlap to form" +
                        " assembly intervals)." +
                        " After training the classifier, the model file is saved, along with supporting files that are used in unit" +
                        " or integration tests.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class TrainEvidenceFilterClassifier extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(StructuralVariationDiscoveryPipelineSpark.class);

    /*
    @ArgumentCollection
    private final StructuralVariationDiscoveryArgumentCollection.TrainEvidenceFilterArgumentCollection params =
            new StructuralVariationDiscoveryArgumentCollection.TrainEvidenceFilterArgumentCollection();
    */
    /**
     * Demo stuff STARTS here....
     */
    @Argument(doc = "path to SVM-light sparse data files, used to load demo classifier data for testing and training."+
                " Folder should contain one train.txt file and one test.txt file.",
            fullName = "demo-data-dir", optional = true)
    private String demoDataDir = "~/Documents/data/demo_classifier_data/agaricus";

    @Argument(doc="Stop classifier training if score does not improve for this many consecutive rounds.",
            fullName = "early-stopping-rounds", optional = true)
    private final int earlyStoppingRounds = 10;

    @Argument(doc="Train classifier for at most this many rounds.",
            fullName = "max-training-rounds", optional = true)
    private final int maxTrainingRounds = 1000;

    @Argument(doc="Use this metric to evaluate performance of the classifier.",
            fullName = "eval-metric", optional = true)
    private final String evalMetric = "logloss";

    @Argument(doc="Classifier tries to maximize eval-metric if true, minimize it if false.",
            fullName = "maximize-eval-metric", optional = true)
    private final boolean maximizeEvalMetric = false;

    @Argument(doc="Seed for random numbers. If null, initialize randomly",
            fullName = "random-seed", optional = true)
    private final Long seed = 0L;


    enum ClassifierTuningStrategy { RANDOM }

    abstract class ClassifierTuner {
        protected final Map<String, ClassifierParamRange<?>> tuneParameters;
        protected final int numSamples;
        protected final boolean maximizeEvalMetric;
        protected final List<Map<String, Object>> hyperparameterSets;
        protected final List<Double> hyperparameterScores;

        ClassifierTuner(final Map<String, ClassifierParamRange<?>> tuneParameters, final int numSamples,
                        final boolean maximizeEvalMetric) {
            if(numSamples < 1) {
                throw new IllegalArgumentException("numSamples (" + numSamples + ") must be >= 1");
            }
            this.tuneParameters = tuneParameters;
            this.numSamples = numSamples;
            this.hyperparameterSets = new ArrayList<> ();
            this.hyperparameterScores = new ArrayList<>();
            this.maximizeEvalMetric = maximizeEvalMetric;
        }

        abstract protected Map<String, Object> chooseNextHyperparameters();

        Map<String, Object> getBestParameters(final Map<String, Object> classifierParameters, final DMatrix trainMatrix,
                                              final TrainTestSplit[] splits, final int maxTrainingRounds,
                                              final int earlyStoppingRounds) {
            Map<String, Object> bestParameters = null;
            double bestScore = maximizeEvalMetric ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
            for(int i = 0; i < numSamples; ++i) {
                final Map<String, Object> hyperparameters = chooseNextHyperparameters();
                // need to store classifierParameters, update with hyperparameters
                classifierParameters.putAll(hyperparameters);
                final double score = getTrainingScore(classifierParameters, trainMatrix, splits, maxTrainingRounds,
                        earlyStoppingRounds, maximizeEvalMetric);
                if(maximizeEvalMetric ? score > bestScore : score < bestScore) {
                    bestScore = score;
                    bestParameters = classifierParameters;
                }
            }
            return bestParameters;
        }
    }

    class RandomClassifierTuner extends ClassifierTuner {

    }

    interface ClassifierParamRange<T> {
        T[] getRandomSamples(final Random random, final int numSamples);
    }

    class ClassifierLinearParamRange implements ClassifierParamRange<Double> {
        private final double low;
        private final double high;

        ClassifierLinearParamRange(double low, double high) {
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

    class ClassifierLogParamRange implements ClassifierParamRange<Double> {
        private final double low;
        private final double high;

        ClassifierLogParamRange(double low, double high) {
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

    class ClassifierIntegerLinearParamRange implements ClassifierParamRange<Integer> {
        private final int low;
        private final int high;

        ClassifierIntegerLinearParamRange(int low, int high) {
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
            final double delta = (high - low) / (numSamples - 1);
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

    class ClassifierIntegerLogParamRange implements ClassifierParamRange<Integer> {
        private final int low;
        private final int high;

        ClassifierIntegerLogParamRange(int low, int high) {
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
                    samples[0] = (int)Math.round(Math.sqrt(high * low));
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
                samples[permutation[i]] = (int)Math.round(val);
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }


    @Argument(doc="Tuning strategy for choosing classifier hyperparameters",
            fullName = "classifier-tuning-strategy", optional = true)
    final ClassifierTuningStrategy classifierTuningStrategy = ClassifierTuningStrategy.RANDOM;

    private final Random random = (seed == null ? new Random() : new Random(seed));
    /**
     * Demo stuff ENDS here....
     */

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        localLogger.info("Loading demo data");
        final DMatrix trainMatrix = loadDMatrix(demoDataDir + "/train.svm.txt");
        final DMatrix testMatrix = loadDMatrix(demoDataDir + "/test.svm.txt");
        // To-do: set scalePosWeight

        localLogger.info("Setting fit parameters");
        //specify parameters
        @SuppressWarnings("serial")
        Map<String, Object> classifierParameters = new HashMap<String, Object>() {
            {
                put("eta", 1.0);
                put("max_depth", 2);
                put("silent", 1);
                put("objective", "binary:logistic");
                put("nthread", 4);
                put("eval_metric", evalMetric);
                put("seed", (seed == null ? random.nextLong() : seed));
            }
        };

        localLogger.info("Setting watches");
        // specify data sets to evaluate
        @SuppressWarnings("serial")
        Map<String, DMatrix> watches = new HashMap<String, DMatrix> () {
            {
                put("train", trainMatrix);
                put("test", testMatrix);
            }
        };

        localLogger.info("Fitting training data");
        final Booster booster;
        try {
            booster = XGBoost.train(trainMatrix, classifierParameters, maxTrainingRounds, watches, null, null);
        } catch (XGBoostError err) {
            throw new GATKException(err.getMessage());
        }

        localLogger.info("Predicting class labels on test data");
        final float[] predictedTestLabels = predictLabel(booster, testMatrix);

        localLogger.info("Calculating prediction accuracy");
        final double accuracy;
        try {
             accuracy = predictionAccuracy(predictedTestLabels, testMatrix.getLabel());
        } catch (XGBoostError err) {
            throw new GATKException(err.getMessage());
        }
        localLogger.info("Accuracy predicted on test set = " + accuracy);

        final float[] trainingTrace = getTrainingTrace(classifierParameters, trainMatrix, testMatrix, maxTrainingRounds, earlyStoppingRounds, maximizeEvalMetric);
        localLogger.info("trainingTrace:");
        for(int index = 0; index < trainingTrace.length; ++index) {
            localLogger.info(index + ": " + trainingTrace[index]);
        }
    }

    private DMatrix loadDMatrix(final String dataPath) {
        try {
            return new DMatrix(dataPath);
        } catch(XGBoostError err) {
            throw new GATKException(err.getMessage());
        }
    }

    private float[] predictLabel(final Booster booster, final DMatrix dmatrix) {
        final float [][] predictions;
        try {
            predictions = booster.predict(dmatrix);
        } catch (XGBoostError err) {
            throw new GATKException(err.getMessage());
        }
        final float [] predictedLabels = new float[predictions.length];
        if(predictions.length == 0) {
            return predictedLabels;
        }
        final int numColumns = predictions[0].length;
        if(numColumns == 1) {
            // binary classifier, reporting only probability of class == 1 (or "true")
            for(int row = 0; row < predictions.length; ++row) {
                predictedLabels[row] = (float)(predictions[row][0] >= 0.5 ? 1.0 : 0.0);
            }

        } else {
            // multiclass classifier (or at binary independently reporting probability of class 0 or 1)
            for (int row = 0; row < predictions.length; ++row) {
                predictedLabels[row] = (float)argmax(predictions[row]);
            }
        }
        return predictedLabels;
    }

    private static Integer[] getRange(final Integer numElements) {
        final Integer[] range = new Integer[numElements];
        for(Integer i = 0; i < numElements; ++i) {
            range[i] = i;
        }
        return range;
    }

    private static int[] getRange(final int numElements) {
        final int[] range = new int[numElements];
        for(int i = 0; i < numElements; ++i) {
            range[i] = i;
        }
        return range;
    }

    private static int[] slice(final int[] arr, final int[] indices) {
        final int[] sliced_arr = new int[indices.length];
        for(int i = 0; i < indices.length; ++i) {
            sliced_arr[i] = arr[indices[i]];
        }
        return sliced_arr;
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

    private static int[] argsort(final int[] arr) {
        final Integer[] sortIndices = getRange((Integer)arr.length);
        Arrays.sort(sortIndices, Comparator.comparingInt(ind -> arr[ind]));
        return ArrayUtils.toPrimitive(sortIndices);
    }

    private static int[] getRandomPermutation(final Random random, final int numElements) {
        // Knuth shuffle
        final int[] permutation = getRange(numElements);
        for(int i = numElements; i > 0; --i) {
            final int swap_ind = random.nextInt(i);
            final int swap_val = permutation[swap_ind];
            permutation[swap_ind] = permutation[i];
            permutation[i] = swap_val;
        }
        return permutation;
    }

    private static double predictionAccuracy(final float[] predictedLabels, final float[] correctLabels) {
        int numCorrect = 0;
        for(int row = 0; row < correctLabels.length; ++row) {
            if(predictedLabels[row] == correctLabels[row]) {
                ++numCorrect;
            }
        }
        return numCorrect / (double)correctLabels.length;
    }
    
    private Map<String, Object> tuneClassifierParameters(final Map<String, Object> classifierParameters,
                                                         final Map<String, ClassifierParamRange<?>> tuneClassifierParameters,
                                                         final ClassifierTuningStrategy classifierTuningStrategy,
                                                         final DMatrix trainMatrix, final int[] stratify,
                                                         final int numCrossvalidationFolds, final int maxTrainingRounds,
                                                         final int earlyStoppingRounds, final boolean maximizeEvalMetric,
                                                         final int numTuningSamples) {
        final List<TrainTestSplit> splits = new ArrayList<>();
        TrainTestSplit.getCrossvalidationSplits(trainMatrix, numCrossvalidationFolds, random, stratify)
                .forEachRemaining(splits::add);

        final ClassifierTuner classifierTuner;
        switch(classifierTuningStrategy) {
            case RANDOM:
                classifierTuner = new RandomClassifierTuner(tuneClassifierParameters, numTuningSamples, maximizeEvalMetric);

            default:
                throw new IllegalStateException("Invalid ClassifierTuningStrategy: " + classifierTuningStrategy);
        }

        return classifierTuner.getBestParameters(classifierParameters, trainMatrix, splits, maxTrainingRounds, earlyStoppingRounds);

    }


    private double getTrainingScore(final Map<String, Object> classifierParameters,
                                  final DMatrix trainMatrix, final TrainTestSplit[] splits,
                                  final int maxTrainingRounds, final int earlyStoppingRounds,
                                  final boolean maximizeEvalMetric) {
        final float[][] trainingTraces = getCrossvalidatedTrainingTraces(
                classifierParameters, trainMatrix, splits, maxTrainingRounds, earlyStoppingRounds, maximizeEvalMetric
        );

        // find the index with the best total (i.e. mean) score across rounds. This yields the optimal number of
        // rounds of training
        float bestTotalScore = maximizeEvalMetric ? Float.MIN_VALUE : Float.MAX_VALUE;
        int bestRoundIndex = -1;
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

        classifierParameters.put("training_rounds", numTrainingRounds);
        return trainingScore;
    }


    private float[][] getCrossvalidatedTrainingTraces(final Map<String, Object> classifierParameters,
                                                      final DMatrix trainMatrix, final TrainTestSplit[] splits,
                                                      final int maxTrainingRounds, final int earlyStoppingRounds,
                                                      final boolean maximizeEvalMetric) {
        final int numCrossvalidationFolds = splits.length;
        float[][] trainingTraces = new float[numCrossvalidationFolds][];
        for(int fold = 0; fold < numCrossvalidationFolds; ++fold) {
            final TrainTestSplit split = splits[fold];
            try {
                trainingTraces[fold] = getTrainingTrace(
                        classifierParameters, trainMatrix.slice(split.trainRows), trainMatrix.slice(split.testRows),
                        maxTrainingRounds, earlyStoppingRounds, maximizeEvalMetric
                );
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
            }

        }
        return trainingTraces;
    }


    private float[] getTrainingTrace(final Map<String, Object> classifierParameters,
                                     final DMatrix trainMatrix, final DMatrix testMatrix,
                                     final int maxTrainingRounds, final int earlyStoppingRounds,
                                     final boolean maximizeEvalMetric) {
        final DMatrix[] evalMatrices = {testMatrix};
        final String[] evalNames = {"test"};
        final float[] metricsOut = new float[1];
        final float[] trainingTrace = new float [maxTrainingRounds];
        float bestTraceValue;
        int stopRound = earlyStoppingRounds;
        try {
            Map<String, DMatrix> watches = new HashMap<> ();
            final Booster booster = XGBoost.train(trainMatrix, classifierParameters, 1, watches, null, null);
            booster.evalSet(evalMatrices, evalNames, 0, metricsOut);
            trainingTrace[0] = metricsOut[0];
            bestTraceValue = trainingTrace[0];
            for(int trainingRound = 1; trainingRound < maxTrainingRounds; ++trainingRound) {
                booster.update(trainMatrix, trainingRound);
                booster.evalSet(evalMatrices, evalNames, trainingRound, metricsOut);
                trainingTrace[trainingRound] = metricsOut[0];
                if(maximizeEvalMetric ? metricsOut[0] > bestTraceValue : metricsOut[0] < bestTraceValue) {
                    // got new bestVal
                    bestTraceValue = metricsOut[0];
                    stopRound = trainingRound + earlyStoppingRounds;
                } else if(trainingRound >= stopRound) {
                    // condition for early stopping has been met. Fill out remaining trace with the most recent value
                    for(int setIndex = trainingRound + 1; setIndex < maxTrainingRounds; ++setIndex) {
                        trainingTrace[setIndex] = trainingTrace[trainingRound];
                    }
                    break;
                }
            }
        } catch(XGBoostError err) {
            throw new GATKException(err.getMessage());
        }
        return trainingTrace;
    }


    static class TrainTestSplit {
        final int[] trainRows;
        final int[] testRows;

        TrainTestSplit(final int[] trainRows, final int[] testRows) {
            this.trainRows = trainRows;
            this.testRows = testRows;
        }

        static TrainTestSplit getTrainTestSplit(final DMatrix trainMatrix, final double trainingFraction,
                                                final Random random, final int[] stratify) {
            final int numRows;
            try {
                numRows = stratify == null ? (int)trainMatrix.rowNum() : stratify.length;
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
            }
            if(trainingFraction < 0.5 * numRows) {
                if(trainingFraction < 0) {
                    throw new IllegalArgumentException("trainingFraction (" + trainingFraction + ") must be in range [0, 1]");
                }
                return new TrainTestSplit(new int[0], getRange(numRows));
            } else if(trainingFraction > 1.0 - 0.5 * numRows) {
                if(trainingFraction > 1) {
                    throw new IllegalArgumentException("trainingFraction (" + trainingFraction + ") must be in range[0, 1]");
                }
                return new TrainTestSplit(getRange(numRows), new int[0]);
            }
            final int[] split_index_ordering = stratify == null ?
                    getRandomPermutation(random, numRows)
                    : TrainTestSplit.getStratfiedIndexOrdering(random, stratify);
            final int numTrain = (int)Math.round(numRows * trainingFraction);
            final int numTest = numRows - numTrain;
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
            return new TrainTestSplit(trainRows, testRows);
        }

        static Iterator<TrainTestSplit> getCrossvalidationSplits(final DMatrix trainMatrix, final int numCrossvalidationFolds,
                                                                 final Random random, final int[] stratify) {
            if(numCrossvalidationFolds < 2) {
                throw new IllegalArgumentException("numCrossvalidationFolds (" + numCrossvalidationFolds + ") must be >= 2");
            }
            final int numRows;
            try {
                numRows = stratify == null ? (int)trainMatrix.rowNum() : stratify.length;
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
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
                final int numTrain = (split_index_ordering.length - 1 - fold) / numFolds;
                final int numTest = split_index_ordering.length - numTrain;
                int[] trainRows = new int[numTrain];
                int[] testRows = new int[numTest];
                int testIndex;
                int orderingIndex;
                if(fold > 0) {
                    for(testIndex = 0; testIndex < fold; ++testIndex) {
                        testRows[testIndex] = split_index_ordering[testIndex];
                    }
                    orderingIndex = testIndex;
                } else {
                    orderingIndex = 0;
                    testIndex = 0;
                }
                for(int trainIndex = 0; trainIndex < trainRows.length; ++trainIndex) {
                    trainRows[trainIndex] = split_index_ordering[orderingIndex];
                    final int orderingStop = Math.min(orderingIndex + numFolds, split_index_ordering.length);
                    for(++orderingIndex; orderingIndex < orderingStop; ++orderingIndex, ++testIndex) {
                        testRows[testIndex] = split_index_ordering[orderingIndex];
                    }
                }

                ++fold;
                return new TrainTestSplit(trainRows, testRows);
            }
        }
    }
}
