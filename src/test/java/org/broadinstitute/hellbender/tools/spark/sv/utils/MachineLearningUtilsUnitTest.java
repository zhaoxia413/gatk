package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import javax.crypto.Mac;
import java.io.IOException;
import java.util.Random;


public class MachineLearningUtilsUnitTest extends GATKBaseTest {
    private static final Random random = new Random();
    // arbitrary number to test some of the array-building / manipulating routines that should not expected to be
    // brittle to differences between positive numbers. Choosing a large one to verify that it scales to reasonable size
    // data sets.
    private static final int NUM_ELEMENTS_TEST = 11001001;

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

        // get two differently-sied arbitrary arrays of indices that are a subset of valid indices into permutation.
        final int[] sliceArr1 = MachineLearningUtils.slice(permutation, MachineLearningUtils.getRange(NUM_ELEMENTS_TEST / 2));
        final int[] sliceArr2 = MachineLearningUtils.slice(permutation, MachineLearningUtils.getRange(NUM_ELEMENTS_TEST / 4));
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
}
