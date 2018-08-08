package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
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
            fullName = "maximize-eval-metric", optional = false)
    private final boolean maximizeEvalMetric = false;

    @Argument(doc="Seed for random numbers. If null, initialize randomly",
            fullName = "random-seed", optional = true)
    private final Long seed = 0L;


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
        Map<String, Object> fitParams = new HashMap<String, Object>() {
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
            booster = XGBoost.train(trainMatrix, fitParams, maxTrainingRounds, watches, null, null);
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

        final float[] trainingTrace = getTrainingTrace(fitParams, trainMatrix, testMatrix, maxTrainingRounds, earlyStoppingRounds, maximizeEvalMetric);
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

    private float[] getTrainingTrace(final Map<String, Object> fitParams,
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
            final Booster booster = XGBoost.train(trainMatrix, fitParams, 1, watches, null, null);
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

    private float[][] getCrossvalidatedTrainingTraces(final Map<String, Object> fitParams,
                                                      final DMatrix trainMatrix, final TrainTestSplit[] splits,
                                                      final int maxTrainingRounds, final int earlyStoppingRounds,
                                                      final boolean maximizeEvalMetric) {
        final int numCrossvalidationFolds = splits.length;
        float[][] trainingTraces = new float[numCrossvalidationFolds][];
        for(int fold = 0; fold < numCrossvalidationFolds; ++fold) {
            final TrainTestSplit split = splits[fold];
            try {
                trainingTraces[fold] = getTrainingTrace(
                        fitParams, trainMatrix.slice(split.trainRows), trainMatrix.slice(split.testRows),
                        maxTrainingRounds, earlyStoppingRounds, maximizeEvalMetric
                );
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
            }

        }
        return trainingTraces;
    }

    private void setTrainingScore(final Map<String, Object> fitParams,
                                   final DMatrix trainMatrix, final TrainTestSplit[] splits,
                                   final int maxTrainingRounds, final int earlyStoppingRounds,
                                   final boolean maximizeEvalMetric) {
        final float[][] trainingTraces = getCrossvalidatedTrainingTraces(
                fitParams, trainMatrix, splits, maxTrainingRounds, earlyStoppingRounds, maximizeEvalMetric
        );
        float bestScore = maximizeEvalMetric ? Float.MIN_VALUE : Float.MAX_VALUE;
        int bestRoundIndex = -1;
        for(int roundIndex = 0; roundIndex < maxTrainingRounds; ++roundIndex) {
            float roundScore = trainingTraces[0][roundIndex];
            for(int fold = 1; fold < trainingTraces.length; ++fold) {
                roundScore += trainingTraces[fold][roundIndex];
            }
            if(maximizeEvalMetric ? roundScore > bestScore : roundScore < bestScore) {
                // the score at this round of training is the best so far
                bestScore = roundScore;
                bestRoundIndex = roundIndex;
            }
        }
        fitParams.put("training_score", bestScore);
        fitParams.put("training_rounds", bestRoundIndex + 1);
    }

    private TrainTestSplit[] getCrossvalidationSplits(final DMatrix trainMatrix, final int numCrossvalidationFolds, final int[] stratify) {
        return null;
    }

    private TrainTestSplit getTrainTestSplit(final DMatrix trainMatrix, final double trainingFraction, final int[] stratify) {
        final int numRows;
        try {
            numRows = stratify == null ? (int)trainMatrix.rowNum() : stratify.length;
        } catch(XGBoostError err) {
            throw new GATKException(err.getMessage());
        }
        final int[] stratify_inds = stratify == null ?
                getRandomPermutation(random, numRows)
                : TrainTestSplit.getStratfiedIndexOrdering(random, stratify);
        final double testingFraction = 1.0 - trainingFraction;
        int nTrain = (int)Math.floor(numRows * trainingFraction);
        int nTest = (int)Math.floor(numRows * testingFraction);
        if(nTrain + nTest < numRows) {
            if(trainingFraction * (nTest + 1) >= testingFraction * (nTrain + 1)) {
                ++nTrain;
            } else {
                ++nTest;
            }
        }
        final List<Integer> trainList = new ArrayList<>(nTrain);
        final List<Integer> testList = new ArrayList<>(nTest);
        int nextNumTrain = 1;
        int nextNumTest = 1;
        for(int i = 0; i < numRows; ++i) {
            if(trainingFraction * nextNumTest >= testingFraction * nextNumTrain) {
                // training set gets next index
                trainList.add(stratify_inds[i]);
                ++nextNumTrain;
            } else {
                testList.add(stratify_inds[i]);
                ++nextNumTest;
            }
        }
        return new TrainTestSplit(
                ArrayUtils.toPrimitive(trainList.toArray(new Integer[0])),
                ArrayUtils.toPrimitive(testList.toArray(new Integer[0]))
        );
    }

    static class TrainTestSplit {
        final int[] trainRows;
        final int[] testRows;

        TrainTestSplit(final int[] trainRows, final int[] testRows) {
            this.trainRows = trainRows;
            this.testRows = testRows;
        }

        static int[] getStratfiedIndexOrdering(final Random random, final int[] stratify) {
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
    }
}
