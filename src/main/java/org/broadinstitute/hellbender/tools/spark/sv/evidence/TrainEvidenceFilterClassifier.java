package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import ml.dmlc.xgboost4j.java.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.MachineLearningUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.XGBoostUtils;

import java.util.*;


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

    @Argument(doc="When performing cross-validation, use this many folds.",
            fullName = "num-crossvalidation-folds", optional = true)
    private final int numCrossvalidationFolds = 5;

    @Argument(doc="When optimizing hyperparameters, search for this many rounds.",
            fullName = "num-hyperparameter-optimization-rounds", optional = true)
    private final int numTuningRounds = 100;

    @Argument(doc="Use this metric to evaluate performance of the classifier.",
            fullName = "eval-metric", optional = true)
    private final String evalMetric = "logloss";

    @Argument(doc="Classifier tries to maximize eval-metric if true, minimize it if false.",
            fullName = "maximize-eval-metric", optional = true)
    private final boolean maximizeEvalMetric = false;

    @Argument(doc="Seed for random numbers. If null, initialize randomly",
            fullName = "random-seed", optional = true)
    private final Long seed = 0L;

    @Argument(doc="Number of threads to use for training classifier. It is optimal to have one thread per available" +
                  " physical processor (not hyperthread)",
            fullName = "nthread", optional = true)
    private final int nthread = Runtime.getRuntime().availableProcessors();




    @Argument(doc="Tuning strategy for choosing classifier hyperparameters",
            fullName = "classifier-tuning-strategy", optional = true)
    final MachineLearningUtils.ClassifierTuningStrategy classifierTuningStrategy
            = MachineLearningUtils.ClassifierTuningStrategy.RANDOM;

    private final Random random = (seed == null ? new Random() : new Random(seed));
    /**
     * Demo stuff ENDS here....
     */

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        localLogger.info("Loading demo data");
        final XGBoostUtils.GATKDMatrix trainMatrix = XGBoostUtils.loadSvmFile(demoDataDir + "/train.svm.txt");
        final XGBoostUtils.GATKDMatrix testMatrix = XGBoostUtils.loadSvmFile(demoDataDir + "/test.svm.txt");
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
                put("nthread", nthread);
                put("eval_metric", evalMetric);
                put("seed", (seed == null ? random.nextLong() : seed));
            }
        };

        // specify tunable parameters
        @SuppressWarnings("serial")
        Map<String, MachineLearningUtils.ClassifierParamRange<?>> tuneClassifierParameters
                = new HashMap<String, MachineLearningUtils.ClassifierParamRange<?>>() {
            {
                put("eta", new MachineLearningUtils.ClassifierLogParamRange(0.01, 1.0));
                put("max_depth", new MachineLearningUtils.ClassifierIntegerLinearParamRange(2, 10));
            }
        };

        localLogger.info("Setting watches");
        // specify data sets to evaluate
        @SuppressWarnings("serial")
        Map<String, DMatrix> watches = new HashMap<String, DMatrix> () {
            {
                put("test", testMatrix.dMatrix);
            }
        };

        localLogger.info("Fitting training data");
        final Booster booster;
        float[][] metrics = new float[1][maxTrainingRounds];
        try {
            booster = XGBoost.train(trainMatrix.dMatrix, classifierParameters, maxTrainingRounds, watches, metrics, null, null, earlyStoppingRounds);
        } catch (XGBoostError err) {
            throw new GATKException(err.getMessage());
        }

        localLogger.info("Fitting training data with new interface");
        MachineLearningUtils.GATKClassifier classifier = new XGBoostUtils.GATKXGBooster();
        classifier.chooseNumThreads(classifierParameters, "nthread", trainMatrix);
        localLogger.info("chose nthread = " + (int)classifierParameters.get("nthread"));
        classifierParameters.put(MachineLearningUtils.NUM_TRAINING_ROUNDS_KEY, maxTrainingRounds);
        classifier.train(classifierParameters, trainMatrix);
        classifierParameters.remove(MachineLearningUtils.NUM_TRAINING_ROUNDS_KEY);


        localLogger.info("Predicting class labels on test data");
        //final float[] predictedTestLabels = predictLabel(booster, testMatrix);
        final int[] predictedTestLabels = classifier.predictClassLabels(testMatrix);

        localLogger.info("Calculating prediction accuracy");
        final double accuracy = predictionAccuracy(predictedTestLabels, testMatrix.getClassLabels());
        localLogger.info("Accuracy predicted on test set = " + accuracy);

        final float[] trainingTrace = classifier.trainAndReturnQualityTrace(classifierParameters, trainMatrix, testMatrix,
                maxTrainingRounds, earlyStoppingRounds, maximizeEvalMetric);
        localLogger.info("trainingTrace:");
        for(int index = 0; index < trainingTrace.length; ++index) {
            localLogger.info(index + ": " + trainingTrace[index]);
        }

        final int[] stratify = trainMatrix.getClassLabels();
        final Map<String, Object> bestClassifierParameters = classifier.tuneClassifierParameters(
                classifierParameters, tuneClassifierParameters, MachineLearningUtils.ClassifierTuningStrategy.RANDOM,
                trainMatrix, random, stratify, numCrossvalidationFolds, maxTrainingRounds, earlyStoppingRounds,
                numTuningRounds, maximizeEvalMetric
        );
        localLogger.info("bestClassifierParameters: " + bestClassifierParameters.toString());
    }

    private static double predictionAccuracy(final int[] predictedLabels, final int[] correctLabels) {
        int numCorrect = 0;
        for(int row = 0; row < correctLabels.length; ++row) {
            if(predictedLabels[row] == correctLabels[row]) {
                ++numCorrect;
            }
        }
        return numCorrect / (double)correctLabels.length;
    }
}
