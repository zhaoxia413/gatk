package org.broadinstitute.hellbender.tools.spark.sv.utils;


import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;
import org.apache.commons.collections4.map.HashedMap;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.*;


public class XGBoostUtils extends MachineLearningUtils {
    public static final double GUESS_OPTIMAL_NUM_THREADS_PROPORTION = 0.85;

    public static final String LEARNING_RATE_KEY = "eta";
    public static final double DEFAULT_LEARNING_RATE = 0.5;
    public static final String MAX_DEPTH_KEY = "max_depth";
    public static final String GAMMA_KEY = "gamma";
    public static final String MIN_CHILD_WEIGHT_KEY = "min_child_weight";
    public static final String SUBSAMPLE_KEY = "subsample";
    public static final String COLSAMPLE_BY_TREE_KEY = "colsample_by_tree";
    public static final String COLSAMPLE_BY_LEVEL_KEY = "colsample_by_level";
    public static final String MAX_DELTA_STEP_KEY = "max_delta_step";

    public static final int DEFAULT_MAX_DEPTH = 6;
    public static final String SILENT_KEY = "silent";
    public static final int DEFAULT_SILENT = 1;
    public static final String OBJECTIVE_KEY = "objective";
    public static final String DEFAULT_OBJECTIVE_VALUE = "binary:logistic";
    public static final String EVAL_METRIC_KEY = "eval_metric";
    public static final String DEFAULT_EVAL_METRIC = "auc";
    public static final String NUM_THREADS_KEY = "nthread";
    public static final String SCALE_POS_WEIGHT_KEY = "scale_pos_weight";
    public static final String SEED_KEY = "seed";
    public static final long DEFAULT_SEED = 0L;
    public static final int DEFAULT_NUM_TRAINING_ROUNDS = 1000;
    public static final int DEFAULT_EARLY_STOPPING_ROUNDS = 50;
    @SuppressWarnings("serial")
    public static final Map<String, Object> DEFAULT_CLASSIFIER_PARAMETERS = new HashMap<String, Object>() {
        {
            put(MachineLearningUtils.NUM_TRAINING_ROUNDS_KEY, DEFAULT_NUM_TRAINING_ROUNDS);

            put(LEARNING_RATE_KEY, DEFAULT_LEARNING_RATE);
            put(MAX_DEPTH_KEY, DEFAULT_MAX_DEPTH);

            put(OBJECTIVE_KEY, DEFAULT_OBJECTIVE_VALUE);
            put(EVAL_METRIC_KEY, DEFAULT_EVAL_METRIC);

            put(SILENT_KEY, DEFAULT_SILENT);
            put(SEED_KEY, DEFAULT_SEED);
        }
    };
    @SuppressWarnings("serial")
    public static final Map<String, MachineLearningUtils.ClassifierParamRange<?>> DEFAULT_TUNING_PARAMETERS
            = new HashMap<String, MachineLearningUtils.ClassifierParamRange<?>>() {
        {
            put(LEARNING_RATE_KEY, new MachineLearningUtils.ClassifierLogParamRange(0.01, 10.0));
            put(MAX_DEPTH_KEY, new MachineLearningUtils.ClassifierIntegerLinearParamRange(2, 20));
            put(GAMMA_KEY, new MachineLearningUtils.ClassifierLinearParamRange(0.0, 20.0));
            put(MIN_CHILD_WEIGHT_KEY, new MachineLearningUtils.ClassifierLogParamRange(1.0, 100.0));
            put(SUBSAMPLE_KEY, new MachineLearningUtils.ClassifierLinearParamRange(0.5, 1.0));
            put(COLSAMPLE_BY_TREE_KEY, new MachineLearningUtils.ClassifierLinearParamRange(0.5, 1.0));
            put(COLSAMPLE_BY_LEVEL_KEY, new MachineLearningUtils.ClassifierLinearParamRange(0.5, 1.0));
            put(MAX_DELTA_STEP_KEY, new MachineLearningUtils.ClassifierLinearParamRange(0.0, 10.0));
        }
    };
    @SuppressWarnings("serial")
    public static final Map<String, Boolean> GET_BUILTIN_MAXIMIZE_EVAL_METRIC = new HashMap<String, Boolean>() {
        {
            put("rmse", false);
            put("mae", false);
            put("logloss", false);
            put("error", false);
            put("merror", false);
            put("auc", true);
            put("ndcg", true);
            put("ndcg-", true);
            put("map", true);
            put("map-", true);
            put("poisson-nloglik", false);
            put("gamma-nloglik", false);
            put("cox-nloglik", false);
            put("gamma-deviance", false);
            put("tweedie-nloglik", false);
        }
    };

    public static Booster getBooster() {
        final Booster booster;
        try {
            final DMatrix dMatrix = new DMatrix(new float[0], 0, 0);
            final Map<String, Object> params = new HashMap<>();
            final Map<String, DMatrix> watches = new HashedMap<>();
            booster = XGBoost.train(dMatrix, params, 0, watches, null, null);
        } catch(XGBoostError err) {
            throw new GATKException(err.getMessage());
        }
        return booster;
    }

    static DMatrix realMatrixToDMatrix(final RealMatrix realMatrix) {
        final int numRows = realMatrix.getRowDimension();
        final int numColumns = realMatrix.getColumnDimension();
        final int numDataColumns = numColumns  - 1;
        final int numData = numRows * numDataColumns;
        final float[] values = new float[numData];
        final float[] classLabels = new float[numRows];
        int index = 0;
        for(int row = 0; row < numRows; ++row) {
            for(int column = 0; column < numColumns; ++column) {
                if(column == MachineLearningUtils.CLASS_LABEL_COLUMN) {
                    classLabels[row] = (float) realMatrix.getEntry(row, column);
                } else {
                    values[index] = (float) realMatrix.getEntry(row, column);
                    ++index;
                }
            }
        }
        try {
            DMatrix dMatrix = new DMatrix(values, numRows, numDataColumns, Float.NaN);
            dMatrix.setLabel(classLabels);
            return dMatrix;
        } catch(XGBoostError err) {
            throw new GATKException(err.getMessage());
        }
    }

    public static class GATKXGBooster extends GATKClassifier {
        private static final long serialVersionUID = 1L;
        private Booster booster;

        public GATKXGBooster() {
            this.booster = null;
        }

        public void write(Kryo kryo, Output output) {
            final byte[] bytearray;
            try {
                bytearray = booster == null ? null : booster.toByteArray();
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
            }
            kryo.writeObjectOrNull(output, bytearray, byte[].class);
        }

        public void read(Kryo kryo, Input input) {
            final byte[] bytearray = kryo.readObjectOrNull(input, byte[].class);
            try {
                booster = bytearray == null ? null : XGBoost.loadModel(new ByteArrayInputStream(bytearray));
            } catch(XGBoostError | IOException err) {
                throw new GATKException(err.getMessage());
            }
        }

        @Override
        public GATKClassifier train(final Map<String, Object> classifierParameters, final RealMatrix trainingMatrix) {
            // specify data sets to evaluate
            Map<String, DMatrix> watches = new HashMap<>();

            final Map<String, Object> trainingParams = new HashMap<>(classifierParameters);
            if(!trainingParams.containsKey(NUM_TRAINING_ROUNDS_KEY)) {
                throw new IllegalArgumentException(
                        "classifierParameters must contain key MachineLearningUtils.NUM_TRAINING_ROUNDS_KEY (\""
                                + NUM_TRAINING_ROUNDS_KEY + "\") with int value."
                );
            }
            final int numTrainingRounds = (int)trainingParams.get(NUM_TRAINING_ROUNDS_KEY);
            trainingParams.remove(NUM_TRAINING_ROUNDS_KEY);

            if(!trainingParams.containsKey(SCALE_POS_WEIGHT_KEY)) {
                final int numPositiveClass = Arrays.stream(getClassLabels(trainingMatrix)).sum();
                trainingParams.put(SCALE_POS_WEIGHT_KEY,
                                   (trainingMatrix.getRowDimension() - numPositiveClass) / (float)numPositiveClass);
            }
            try {
                booster = XGBoost.train(realMatrixToDMatrix(trainingMatrix), trainingParams, numTrainingRounds,
                                        watches, null, null);
            } catch (XGBoostError err) {
                throw new GATKException(err.getMessage());
            }
            return this;
        }

        @Override
        public boolean getMaximizeEvalMetric(final Map<String, Object> classifierParameters) {
            final String evalMetric = (String)classifierParameters.get(EVAL_METRIC_KEY);
            final String key = evalMetric.contains("@") ? evalMetric.substring(0, evalMetric.indexOf('@')) : evalMetric;
            if(!GET_BUILTIN_MAXIMIZE_EVAL_METRIC.containsKey(key)) {
                throw new IllegalArgumentException(
                        evalMetric + " not in default GET_BUILTIN_MAXIMIZE_EVAL_METRIC. To use a custom evaluation"
                                + " metric, a boolean must be added to that map specifying if the metric should be maximized"
                                + " (or minimized)"
                );
            }
            return GET_BUILTIN_MAXIMIZE_EVAL_METRIC.get(key);
        }


        @Override
        public double[] trainAndReturnQualityTrace(
                final Map<String, Object> classifierParameters, final RealMatrix trainingMatrix,
                final RealMatrix evaluationMatrix, final int maxTrainingRounds, final int earlyStoppingRounds
        ) {
            final boolean maximizeEvalMetric = getMaximizeEvalMetric(classifierParameters);
            final DMatrix trainingDMatrix = realMatrixToDMatrix(trainingMatrix);
            final DMatrix[] evalMatrices = {realMatrixToDMatrix(evaluationMatrix)};
            final String[] evalNames = {"test"};
            final float[] metricsOut = new float[1];
            final double[] trainingTrace = new double [maxTrainingRounds];
            double bestTraceValue;
            int stopRound = earlyStoppingRounds;
            try {
                Map<String, DMatrix> watches = new HashMap<> ();
                booster = XGBoost.train(trainingDMatrix, classifierParameters, 1, watches, null, null);
                booster.evalSet(evalMatrices, evalNames, 0, metricsOut);
                trainingTrace[0] = metricsOut[0];
                bestTraceValue = trainingTrace[0];
                for(int trainingRound = 1; trainingRound < maxTrainingRounds; ++trainingRound) {
                    booster.update(trainingDMatrix, trainingRound);
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

        @Override
        public double[][] predictProbability(final RealMatrix matrix) {
            final float [][] floatProbabilities;
            try {
                floatProbabilities = booster.predict(realMatrixToDMatrix(matrix));
            } catch(XGBoostError err) {
                throw new GATKException(err.getClass() + ": " + err.getMessage());
            }
            // convert from float[][] to double[][]
            final int numRows = floatProbabilities.length;
            final int numCols = numRows > 0 ? floatProbabilities[0].length : 0;
            final double[][] doubleProbabilities = new double [numRows][numCols];

            for(int row = 0; row < numRows; ++row) {
                final double[] newRow = doubleProbabilities[row];
                final float[] oldRow = floatProbabilities[row];
                for(int column = 0; column < numCols; ++column) {
                    newRow[column] = oldRow[column];
                }
            }
            return doubleProbabilities;
        }
    }
}
