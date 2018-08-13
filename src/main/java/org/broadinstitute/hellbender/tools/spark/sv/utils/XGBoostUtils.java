package org.broadinstitute.hellbender.tools.spark.sv.utils;


import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;
import ml.dmlc.xgboost4j.java.XGBoostError;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.HashMap;
import java.util.Map;


public class XGBoostUtils extends MachineLearningUtils {
    public static GATKDMatrix loadSvmFile(final String fileName) {
        try {
            return new GATKDMatrix(new DMatrix(fileName));
        } catch(XGBoostError err) {
            throw new GATKException(err.getMessage());
        }
    }

    public static class GATKDMatrix implements GATKMatrix {
        public final DMatrix dMatrix;

        GATKDMatrix(final DMatrix dMatrix) {
            this.dMatrix = dMatrix;
        }

        @Override
        public GATKDMatrix sliceRows(final int[] rowIndices) {
            try {
                return new GATKDMatrix(dMatrix.slice(rowIndices));
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
            }
        }

        @Override
        public int getNumRows() {
            try {
                return (int)dMatrix.rowNum();
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
            }
        }

        @Override
        public int[] getClassLabels() {
            final float[] floatClassLabels;
            try {
                floatClassLabels = dMatrix.getLabel();
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
            }
            final int[] classLabels = new int[floatClassLabels.length];
            for(int i = 0; i < classLabels.length; ++i) {
                classLabels[i] = (int)floatClassLabels[i];
            }
            return classLabels;
        }
    }

    public static class GATKXGBooster extends GATKClassifier {
        private Booster booster;

        public GATKXGBooster() {
            this.booster = null;
        }

        public GATKXGBooster(final Booster booster) {
            this.booster = booster;
        }

        @Override
        public GATKClassifier train(final Map<String, Object> classifierParameters, final GATKMatrix trainingMatrix) {
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

            try {
                booster = XGBoost.train(((GATKDMatrix)trainingMatrix).dMatrix, trainingParams, numTrainingRounds,
                                        watches, null, null);
            } catch (XGBoostError err) {
                throw new GATKException(err.getMessage());
            }
            return this;
        }

        @Override
        public float[] trainAndReturnQualityTrace(
                final Map<String, Object> classifierParameters, final GATKMatrix trainingMatrix,
                final GATKMatrix evaluationMatrix, final int maxTrainingRounds, final int earlyStoppingRounds,
                final boolean maximizeEvalMetric
        ) {
            final DMatrix trainingDMatrix = ((GATKDMatrix)trainingMatrix).dMatrix;
            final DMatrix[] evalMatrices = {((GATKDMatrix)evaluationMatrix).dMatrix};
            final String[] evalNames = {"test"};
            final float[] metricsOut = new float[1];
            final float[] trainingTrace = new float [maxTrainingRounds];
            float bestTraceValue;
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
        public float[][] predictProbability(final GATKMatrix matrix) {
            try {
                return booster.predict(((GATKDMatrix) matrix).dMatrix);
            } catch(XGBoostError err) {
                throw new GATKException(err.getMessage());
            }
        }
    }
}
