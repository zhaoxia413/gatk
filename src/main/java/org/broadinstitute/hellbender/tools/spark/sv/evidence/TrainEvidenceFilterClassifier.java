package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.MachineLearningUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.XGBoostUtils;

import java.io.File;
import java.io.IOException;
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
public class TrainEvidenceFilterClassifier extends CommandLineProgram {
    private static final long serialVersionUID = 1L;
    private static final Logger localLogger = LogManager.getLogger(TrainEvidenceFilterClassifier.class);

    private static final String CURRENT_DIRECTORY = System.getProperty("user.dir");
    private static final String gatkDirectory = System.getProperty("gatkdir", CURRENT_DIRECTORY) + "/";
    private static final String publicTestDirRelative = "src/test/resources/";
    public static final String publicTestDir = new File(gatkDirectory, publicTestDirRelative).getAbsolutePath() + "/";


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
            fullName = "demo-data-file", optional = true)
    private String demoDataFile = publicTestDir + "org/broadinstitute/hellbender/tools/spark/sv/utils/agaricus-integers.csv.gz";

    @Argument(doc = "path to SVM-light sparse data files, used to load demo classifier data for testing and training."+
            " Folder should contain one train.txt file and one test.txt file.",
            fullName = "classifier-model-file", optional = false)
    private String classifierModelFile;

    @Argument(doc="Stop classifier training if score does not improve for this many consecutive rounds.",
            fullName = "early-stopping-rounds", optional = true)
    private final int earlyStoppingRounds = XGBoostUtils.DEFAULT_EARLY_STOPPING_ROUNDS;

    @Argument(doc="Train classifier for at most this many rounds.",
            fullName = "max-training-rounds", optional = true)
    private final int maxTrainingRounds = XGBoostUtils.DEFAULT_NUM_TRAINING_ROUNDS;

    @Argument(doc="When performing cross-validation, use this many folds.",
            fullName = "num-crossvalidation-folds", optional = true)
    private final int numCrossvalidationFolds = MachineLearningUtils.DEFAULT_NUM_CROSSVALIDATION_FOLDS;

    @Argument(doc="When optimizing hyperparameters, search for this many rounds.",
            fullName = "num-hyperparameter-optimization-rounds", optional = true)
    private final int numTuningRounds = MachineLearningUtils.DEFAULT_NUM_TUNING_ROUNDS;

    @Argument(doc="When optimizing hyperparameters, reserve this proportion of data for tuning hyperparameters.",
            fullName = "hyperparameter-tuning-proportion", optional = true)
    private Double hyperparameterTuningProportion = null;

    @Argument(doc="Use this metric to evaluate performance of the classifier.",
            fullName = "eval-metric", optional = true)
    private final String evalMetric = XGBoostUtils.DEFAULT_EVAL_METRIC;

    @Argument(doc="Seed for random numbers. If null, initialize randomly",
            fullName = "random-seed", optional = true)
    private final Long seed = XGBoostUtils.DEFAULT_SEED;

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
    //protected void runTool(final JavaSparkContext ctx) {
    protected Object doWork() {
        localLogger.info("Loading demo data");
        final RealMatrix dataMatrix = MachineLearningUtils.loadCsvFile(demoDataFile);

        if(hyperparameterTuningProportion == null) {
            hyperparameterTuningProportion = 1.0 / (1.0 + numTuningRounds);
        }
        //final int[] stratify = MachineLearningUtils.getClassLabels(dataMatrix);
        final int[] stratify = MachineLearningUtils.stratifyMatrixToStratifyArray(dataMatrix, 5,
                (int)Math.ceil(numCrossvalidationFolds / hyperparameterTuningProportion));

        localLogger.info("Splitting data");
        final MachineLearningUtils.TrainTestSplit hyperSplit = MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                hyperparameterTuningProportion, dataMatrix.getRowDimension(), random, stratify
        );
        final RealMatrix tuneMatrix = MachineLearningUtils.sliceRows(dataMatrix, hyperSplit.trainRows);
        final int[] tuneStratify = MachineLearningUtils.slice(stratify, hyperSplit.trainRows);
        final RealMatrix validateMatrix = MachineLearningUtils.sliceRows(dataMatrix, hyperSplit.testRows);
        final int[] validateStratify = MachineLearningUtils.slice(stratify, hyperSplit.testRows);

        // set the number of threads

        final Map<String, Object> classifierParameters = new HashMap<> (XGBoostUtils.DEFAULT_CLASSIFIER_PARAMETERS);
        classifierParameters.put(XGBoostUtils.NUM_THREADS_KEY, nthread);
        classifierParameters.put(XGBoostUtils.EVAL_METRIC_KEY, evalMetric);

        final XGBoostUtils.GATKXGBooster classifier = new XGBoostUtils.GATKXGBooster();
        classifier.chooseNumThreads(classifierParameters, XGBoostUtils.NUM_THREADS_KEY, tuneMatrix);
        localLogger.info("Chose " + classifierParameters.get(XGBoostUtils.NUM_THREADS_KEY) + " threads");

        localLogger.info("Tuning hyperparameters");
        final Map<String, Object> bestClassifierParameters = classifier.tuneClassifierParameters(
                classifierParameters, XGBoostUtils.DEFAULT_TUNING_PARAMETERS, classifierTuningStrategy,
                tuneMatrix, random, tuneStratify, numCrossvalidationFolds, maxTrainingRounds, earlyStoppingRounds,
                numTuningRounds
        );
        localLogger.info("bestClassifierParameters: " + bestClassifierParameters.toString());

        localLogger.info("Cross-val predicting");
        final int[] predictedTestLabels = classifier.crossvalidatePredict(
                validateMatrix, bestClassifierParameters, random, validateStratify, numCrossvalidationFolds
        );

        final double accuracy = MachineLearningUtils.getPredictionAccuracy(
                predictedTestLabels, MachineLearningUtils.getClassLabels(validateMatrix)
        );
        localLogger.info("Crossvalidated accuracy = " + String.format("%.1f%%", 100.0 * accuracy));


        localLogger.info("Training final classifier");
        classifier.train(bestClassifierParameters, dataMatrix);

        localLogger.info("Saving final classifier to " + classifierModelFile);
        try {
            classifier.save(classifierModelFile);
        } catch(IOException err) {
            throw new GATKException(err.getClass() + ": " + err.getMessage());
        }
        return null;
    }

}
