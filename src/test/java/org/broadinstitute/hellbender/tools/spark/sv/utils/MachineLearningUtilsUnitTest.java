package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import javax.crypto.Mac;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class MachineLearningUtilsUnitTest extends GATKBaseTest {
    private static final Random random = new Random();
    // arbitrary number to test some of the array-building / manipulating routines that should not expected to be
    // brittle to differences between positive numbers. Choosing a large one to verify that it scales to reasonable size
    // data sets.
    private static final int NUM_ELEMENTS_TEST = 1001001;
    private static final double TRAINING_FRACTION = 0.3;
    private static final int NUM_CROSSVALIDATION_FOLDS = 5;
    // choose an arbitrary large value in between 1 and NUM_ELEMENTS_TEST
    private static final int NUM_STRATIFY_CLASSES = (int)Math.round(Math.sqrt(NUM_ELEMENTS_TEST));

    @Test(groups = "sv")
    protected void testGetRange() {
        try {
            MachineLearningUtils.getRange(-1);
            Assert.fail("getRange(negativeInteger) should throw IllegalArgumentException");
        } catch (IllegalArgumentException ignored) {
        }
        Integer negativeInteger = -1;
        try {
            MachineLearningUtils.getRange(negativeInteger);
            Assert.fail("getRange(negativeInteger) should throw IllegalArgumentException");
        } catch (IllegalArgumentException ignored) {
        }
        assertEmpty(MachineLearningUtils.getRange(0), "getRange(0) should be empty");
        assertEmpty(MachineLearningUtils.getRange((Integer) 0), "getRange((Integer)0) should be empty");

        // if these work for any given positive integer they'll work for the rest. Just test one.
        final int[] intRange = MachineLearningUtils.getRange(NUM_ELEMENTS_TEST);
        final Integer[] integerRange = MachineLearningUtils.getRange((Integer) NUM_ELEMENTS_TEST);
        Assert.assertEquals(intRange.length, NUM_ELEMENTS_TEST, "getRange(int) returned wrong number of elements");
        Assert.assertEquals(integerRange.length, NUM_ELEMENTS_TEST, "getRange(Integer) returned wrong number of elements");
        for (int j = 0; j < NUM_ELEMENTS_TEST; ++j) {
            Assert.assertEquals(intRange[j], j, "getRange(int) has wrong value");
            Assert.assertEquals(integerRange[j], (Integer) j, "getRange(Integer) has wrong value");
        }
    }

    @Test(groups = "sv")
    protected void testArrayOrderManipulation() {
        try {
            MachineLearningUtils.getRandomPermutation(random, -1);
            Assert.fail("getRandomPermutation(random, negativeInteger) should throw IllegalArgumentException");
        } catch (IllegalArgumentException ignored) {
        }
        assertEmpty(MachineLearningUtils.getRandomPermutation(random, 0), "getRandomPermutation(random, 0) should be empty");

        // if these work for any given positive integer they'll work for the rest. Just test one.
        final int[] permutation = MachineLearningUtils.getRandomPermutation(random, NUM_ELEMENTS_TEST);
        Assert.assertEquals(permutation.length, NUM_ELEMENTS_TEST, "getRandomPermutation() returned wrong number of elements");
        // assert the permutation is non-trivial. **Technically** it's possible this will fail, but it's so unlikely it
        // should never happen during the epoch where the universe is capable of hosting life...
        final int[] range = MachineLearningUtils.getRange(NUM_ELEMENTS_TEST);
        assertArraysNotEqual(permutation, range,"getPermutation returned trivial permutaiton.");

        // use argsort to find indices that invert the permutation.
        final int[] sortedInds = MachineLearningUtils.argsort(permutation);
        // it's easy to use sliceAssign to get the correct inverse permutation
        final int[] inverseInds = new int[NUM_ELEMENTS_TEST];
        MachineLearningUtils.sliceAssign(inverseInds, permutation, range);
        assertArrayEquals(sortedInds, inverseInds, "sortedInds don't equal inverseInds");
        // using slice should put the permutation back into order
        Assert.assertEquals(MachineLearningUtils.slice(permutation, sortedInds), range,
                "argsort did not correctly unscramble the permutation.");
        // check corner-cases
        assertEmpty(MachineLearningUtils.slice(permutation, new int[0]), "empty slice should return empty result");
        assertEmpty(MachineLearningUtils.argsort(new int[0]), "argsort of empty array should return empty result");
        try {
            MachineLearningUtils.slice(permutation, new int[] {NUM_ELEMENTS_TEST});
            Assert.fail("slice should throw ArrayIndexOutOfBoundsException when passed a bad index");
        } catch (ArrayIndexOutOfBoundsException ignored) {
        }
        try {
            MachineLearningUtils.slice(permutation, new int[] {-1});
            Assert.fail("slice should throw ArrayIndexOutOfBoundsException when passed a bad index");
        } catch (ArrayIndexOutOfBoundsException ignored) {
        }

        // get two differently-sized arbitrary arrays of indices that are a subset of valid indices into permutation.
        final int[] sliceArr1 = MachineLearningUtils.slice(
                permutation, MachineLearningUtils.getRange(NUM_ELEMENTS_TEST / 2)
        );
        final int[] sliceArr2 = MachineLearningUtils.slice(
                permutation, MachineLearningUtils.getRange(NUM_ELEMENTS_TEST / 4)
        );
        try {
            MachineLearningUtils.sliceAssign(permutation, sliceArr1, sliceArr2);
            Assert.fail("sliceAssign should throw IllegalArgumentException when slice indices don't have same size as slice values");
        } catch (IllegalArgumentException ignored) {
        }
        try {
            MachineLearningUtils.sliceAssign(permutation, sliceArr2, sliceArr1);
            Assert.fail("sliceAssign should throw IllegalArgumentException when slice indices don't have same size as slice values");
        } catch (IllegalArgumentException ignored) {
        }
        // sliceAssign should throw IndexOutOfBoundsException when asked for a bad index
        try {
            MachineLearningUtils.sliceAssign(permutation, new int[] {NUM_ELEMENTS_TEST}, new int[] {42});
            Assert.fail("sliceAssign should throw ArrayIndexOutOfBoundsException when passed a bad index");
        } catch (ArrayIndexOutOfBoundsException ignored) {
        }
        try {
            MachineLearningUtils.sliceAssign(permutation, new int[] {-1}, new int[] {42});
            Assert.fail("sliceAssign should throw ArrayIndexOutOfBoundsException when passed a bad index");
        } catch (ArrayIndexOutOfBoundsException ignored) {
        }
    }

    @Test(groups = "sv")
    /**
     * Test expected behavior of TrainTestSplit factory functions.
     * Basic functionality is tested in actual use of tuning / cross-validating. Here check that the details match
     * expectations (data set is divided evenly according to specification).
     */
    protected void testTrainTestSplit() {
        // Create stratify array, random array of class membership indexes used to stratify splits. Make it class
        // populations uneven to further stress-test algorithms.
        final int[] stratify = Arrays.stream(
                new MachineLearningUtils.ClassifierIntegerLogParamRange(1, NUM_STRATIFY_CLASSES)
                .getRandomSamples(random, NUM_ELEMENTS_TEST)
        ).mapToInt(i->i).toArray();
        Assert.assertEquals(stratify.length, NUM_ELEMENTS_TEST, "stratify was generated improperly");

        // Test stratified and un-stratified train-test split
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        TRAINING_FRACTION, NUM_ELEMENTS_TEST, random, stratify
                ),
                NUM_ELEMENTS_TEST, TRAINING_FRACTION, stratify
        );

        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        TRAINING_FRACTION, NUM_ELEMENTS_TEST, random, null
                ),
                NUM_ELEMENTS_TEST, TRAINING_FRACTION, null
        );


        // Test stratified and un-stratified cross-validation splits
        final double cvTrainingFraction = 1.0 - 1.0 / (double)NUM_CROSSVALIDATION_FOLDS;

        final Iterator<MachineLearningUtils.TrainTestSplit> stratifiedCvSplits
                = MachineLearningUtils.TrainTestSplit.getCrossvalidationSplits(
                        NUM_CROSSVALIDATION_FOLDS, NUM_ELEMENTS_TEST, random, stratify
            );
        while(stratifiedCvSplits.hasNext()) {
            final MachineLearningUtils.TrainTestSplit stratifiedCvSplit = stratifiedCvSplits.next();
            assertGoodSplit(stratifiedCvSplit, NUM_ELEMENTS_TEST, cvTrainingFraction, stratify);
        }

        final Iterator<MachineLearningUtils.TrainTestSplit> flatCvSplits
                = MachineLearningUtils.TrainTestSplit.getCrossvalidationSplits(
                NUM_CROSSVALIDATION_FOLDS, NUM_ELEMENTS_TEST, random, null
        );
        while(flatCvSplits.hasNext()) {
            final MachineLearningUtils.TrainTestSplit flatCvSplit = flatCvSplits.next();
            assertGoodSplit(flatCvSplit, NUM_ELEMENTS_TEST, cvTrainingFraction, null);
        }

        // check corner-cases
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        0.0, NUM_ELEMENTS_TEST, random, stratify
                ),
                NUM_ELEMENTS_TEST, 0.0, stratify
        );
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        1.0, NUM_ELEMENTS_TEST, random, stratify
                ),
                NUM_ELEMENTS_TEST, 1.0, stratify
        );
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        0.0, NUM_ELEMENTS_TEST, random, null
                ),
                NUM_ELEMENTS_TEST, 0.0, null
        );
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        1.0, NUM_ELEMENTS_TEST, random, null
                ),
                NUM_ELEMENTS_TEST, 1.0, null
        );

        // ensure errors are thrown when crazy values are passed
        try {
            MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                    -0.01, NUM_ELEMENTS_TEST, random, stratify
            );
            Assert.fail("getTrainTestSplit should throw IllegalArgmentException when trainingFraction is not in range [0, 1]");
        } catch (IllegalArgumentException ignored) {
        }
        try {
            MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                    1.01, NUM_ELEMENTS_TEST, random, stratify
            );
            Assert.fail("getTrainTestSplit should throw IllegalArgmentException when trainingFraction is not in range [0, 1]");
        } catch (IllegalArgumentException ignored) {
        }
        try {
            MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                    -0.01, NUM_ELEMENTS_TEST, random, null
            );
            Assert.fail("getTrainTestSplit should throw IllegalArgmentException when trainingFraction is not in range [0, 1]");
        } catch (IllegalArgumentException ignored) {
        }
        try {
            MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                    1.01, NUM_ELEMENTS_TEST, random, null
            );
            Assert.fail("getTrainTestSplit should throw IllegalArgmentException when trainingFraction is not in range [0, 1]");
        } catch (IllegalArgumentException ignored) {
        }

        try {
            MachineLearningUtils.TrainTestSplit.getCrossvalidationSplits(
                    1, NUM_ELEMENTS_TEST, random, stratify
            );
            Assert.fail("getCrossvalidationSplits should throw IllegalArgmentException when numCrossvalidationFolds < 2");
        } catch (IllegalArgumentException ignored) {
        }
    }

    private static void assertEmpty(final int[] arr, final String message) {
        if(arr.length != 0) {
            Assert.fail(message);
        }
    }

    private static <T> void assertEmpty(final T[] arr, final String message) {
        if(arr.length != 0) {
            Assert.fail(message);
        }
    }

    private static void assertArraysNotEqual(final int[] actuals, final int[] expecteds, final String message) {
        if(actuals.length != expecteds.length) {
            return;
        }
        for (int index = 0; index < expecteds.length; ++index) {
            if(actuals[index] != expecteds[index]) {
                return;
            }
        }
        Assert.fail(message);
    }

    private static void assertArrayEquals(final int[] actuals, final int[] expecteds, final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for (int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], "at index=" + index + ": " + message);
        }
    }

    private static void assertArrayEquals(final int[] actuals, final double[] expecteds, final double tol,
                                          final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for (int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], tol, "at index=" + index + ": " + message);
        }
    }

    private static void assertArraysDisjoint(final int[] arr1, final int[] arr2, final String message) {
        final Set<Integer> arr1Values = Arrays.stream(arr1).boxed().collect(Collectors.toSet());
        for(final int val2 : arr2) {
            if(arr1Values.contains(val2)) {
                Assert.fail(message + ": " + val2 + " in both arrays");
            }
        }
    }

    private static void assertGoodSplit(final MachineLearningUtils.TrainTestSplit split, final int numElements,
                                        final double trainingFraction, final int[] stratify) {
        // first check if the split is trivial
        if(trainingFraction == 0.0) {
            assertEmpty(split.trainRows, "trainRows should be empty");
            assertArrayEquals(split.testRows, MachineLearningUtils.getRange(numElements),
                    "testRows should be all elements");
            return;
        } else if(trainingFraction == 1.0) {
            assertEmpty(split.testRows, "testRows should be empty");
            assertArrayEquals(split.trainRows, MachineLearningUtils.getRange(numElements),
                    "trainRows should be all elements");
            return;
        }

        // next check that it *is* a split
        Assert.assertEquals(split.testRows.length + split.trainRows.length, numElements,
                "number of trainRows + number of testRows != number of elements");
        assertArraysDisjoint(split.testRows, split.trainRows, "trainRows and testRows have overlap");
        final IntSummaryStatistics trainStats = Arrays.stream(split.trainRows).summaryStatistics();
        final IntSummaryStatistics testStats = Arrays.stream(split.testRows).summaryStatistics();
        final int minIndex = Math.min(trainStats.getMin(), testStats.getMin());
        final int maxIndex = Math.max(trainStats.getMax(), testStats.getMax());
        Assert.assertEquals(minIndex, 0, "split contains indices less than 0");
        Assert.assertEquals(maxIndex, numElements - 1, "split contains indices past end of data");

        // finally check balance of division
        Assert.assertEquals(split.trainRows.length, trainingFraction * numElements, 2,
                "Number of training rows differs from expected by more than 2");
        if(stratify != null) {
            // check balance of division for each stratify value
            final double[] expectedStratifyCounts = Arrays.stream(
                    getStratifyCounts(stratify)
            ).mapToDouble(c -> c * trainingFraction).toArray();
            final int[] trainStratifyCounts = getStratifyCounts(MachineLearningUtils.slice(stratify, split.trainRows));
            assertArrayEquals(trainStratifyCounts, expectedStratifyCounts, 2,
                    "split wasn't properly balanced for stratify");
        }
    }

    private static int[] getStratifyCounts(final int[] stratify) {
        final Map<Integer, Integer> stratifyCountsMap = new HashMap<>();
        for(final Integer value: stratify) {
            stratifyCountsMap.put(value, stratifyCountsMap.getOrDefault(value, 0) + 1);
        }
        if(stratifyCountsMap.isEmpty()) {
            return new int[0];
        }
        final int maxValue = stratifyCountsMap.keySet().stream().mapToInt(Integer::intValue).max().getAsInt();
        final int[] stratifyCounts = new int[maxValue + 1];
        for(final Map.Entry entry : stratifyCountsMap.entrySet()) {
            stratifyCounts[(int)(Integer)entry.getKey()] = (int)(Integer)entry.getValue();
        }
        return stratifyCounts;
    }
}
