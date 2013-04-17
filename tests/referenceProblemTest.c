/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stReferenceProblem.h"
#include "stMatchingAlgorithms.h"
#include "stCheckEdges.h"

static stList *stubs;
static stList *chains;
static float *zMatrix;
static int64_t nodeNumber;

static void teardown() {
    if (nodeNumber != -1) {
        stList_destruct(stubs);
        stList_destruct(chains);
        free(zMatrix);
        nodeNumber = -1;
    }
}

static void setup() {
    teardown();
    assert(nodeNumber == -1);
    while(nodeNumber % 2 != 0) {
        nodeNumber = st_randomInt(0, 100);
    }
    assert(nodeNumber >= 0);
    assert(nodeNumber % 2 == 0);
    stubs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    chains = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int64_t i=0; i<nodeNumber/2; i++) {
        assert(nodeNumber/2 > 0);
        stIntTuple *edge = stIntTuple_construct2(i, nodeNumber/2 + i);
        if(stList_length(stubs) == 0 || st_random() > 0.9) {
            stList_append(stubs, edge);
        }
        else {
            stList_append(chains, edge);
        }
    }
    zMatrix = st_calloc(nodeNumber*nodeNumber, sizeof(float));
    for(int64_t i=0; i<nodeNumber; i++) {
        for(int64_t j=i+1; j<nodeNumber; j++) {
            double score = st_random();
            zMatrix[i * nodeNumber + j] = score;
            zMatrix[j * nodeNumber + i] = score;
        }
    }
    st_logDebug("To test the adjacency problem we've created a problem with %" PRIi64 " nodes %" PRIi64 " stubs and %" PRIi64 " chains\n", nodeNumber, stList_length(stubs), stList_length(chains));
}

static void checkIsValidReference(CuTest *testCase, stList *reference,
        double totalScore) {
    stList *chosenEdges = convertReferenceToAdjacencyEdges(reference);
    //Check that everyone has a partner.
    CuAssertIntEquals(testCase, nodeNumber, stList_length(chosenEdges) * 2);
    stSortedSet *nodes = stSortedSet_construct3((int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < nodeNumber; i++) {
        stSortedSet_insert(nodes, stIntTuple_construct1( i));
    }
    checkEdges(chosenEdges, nodes, 1, 0);
    //Check that the score is correct
    double totalScore2 = calculateZScoreOfReference(reference, nodeNumber, zMatrix);
    CuAssertDblEquals(testCase, totalScore2, totalScore, 0.000001);
    //Check that the stubs are properly connected.
    stList *allEdges = stList_copy(chosenEdges, NULL);
    stList_appendAll(allEdges, stubs);
    stList_appendAll(allEdges, chains);
    stList *components = getComponents(allEdges);
    CuAssertIntEquals(testCase, stList_length(stubs), stList_length(reference));
    CuAssertIntEquals(testCase, stList_length(stubs), stList_length(components));
    //Cleanup
    stList_destruct(components);
    stSortedSet_destruct(nodes);
    stList_destruct(allEdges);
    stList_destruct(chosenEdges);
}

static void testMakeReferenceGreedily(CuTest *testCase) {
    for (int64_t i = 0; i < 100; i++) {
        setup();
        double totalScore;
        time_t startTime = time(NULL);
        bool fast = 1; //st_random() > 0.5;
        stList *reference = makeReferenceGreedily(stubs, chains, zMatrix,
                nodeNumber, &totalScore, fast);
        checkIsValidReference(testCase, reference, totalScore);
        time_t totalTime = time(NULL) - startTime;
        char *cA = stString_print("Greedy (fast:%" PRIi64 "), %" PRIi64 " nodes, it took %" PRIi64 " seconds, score: %f\n", fast, nodeNumber, totalTime, totalScore);
        st_logInfo(cA);
        logReference(reference, nodeNumber, zMatrix, totalScore, cA);
        free(cA);
        teardown();
    }
}

static void testGibbsSamplingWithSimulatedAnnealing(CuTest *testCase,
        double(*temperatureFn)(double), bool pureGreedy, int64_t maxNumberOfChainsBeforeSwitchingToFast) {
    for (int64_t i = 0; i < 100; i++) {
        setup();
        double totalScore;
        time_t startTime = time(NULL);
        bool fast = maxNumberOfChainsBeforeSwitchingToFast > nodeNumber;
        stList *reference = makeReferenceGreedily(stubs, chains, zMatrix,
                nodeNumber, &totalScore, fast);
        checkIsValidReference(testCase, reference, totalScore);
        time_t totalTime = time(NULL) - startTime;
        char *cA = stString_print("Pre-annealing greedy (fast:%" PRIi64 "), %" PRIi64 " nodes, it took %" PRIi64 " seconds, score: %f\n", fast, nodeNumber, totalTime, totalScore);
        st_logInfo(cA);
        logReference(reference, nodeNumber, zMatrix, totalScore, cA);
        free(cA);
        startTime = time(NULL);
        int64_t permutations = st_randomInt(0, 100);
        if(pureGreedy) {
            greedyImprovement(reference, chains, zMatrix, permutations);
        }
        else {
            gibbsSamplingWithSimulatedAnnealing(reference, chains, zMatrix,
                    permutations, temperatureFn);
        }
        double totalScoreAfterAnnealing = calculateZScoreOfReference(reference, nodeNumber, zMatrix);
        checkIsValidReference(testCase, reference, totalScoreAfterAnnealing);
        totalTime = time(NULL) - startTime;
        cA = stString_print("Post-annealing permutations:%" PRIi64 ", %" PRIi64 " nodes, it took %" PRIi64 " seconds, score: %f\n", permutations, nodeNumber, totalTime, totalScore);
        st_logInfo(cA);
        logReference(reference, nodeNumber, zMatrix, totalScoreAfterAnnealing,
                cA);
        if(pureGreedy) {
            CuAssertTrue(testCase, totalScoreAfterAnnealing + 0.001 >= totalScore);
        }
        teardown();
    }
}

static void testGibbsSamplingWithSimulatedAnnealing_NoExponentiation_Greedy_Fast(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, but greedy sampling\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase,
            NULL, 1, 0);
}

static void testGibbsSamplingWithSimulatedAnnealing_NoExponentiation_Greedy(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, but greedy sampling\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase,
            NULL, 1, INT64_MAX);
}

static void testGibbsSamplingWithSimulatedAnnealing_NoExponentiation(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, but no exponentiation\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase,
            NULL, 0, INT64_MAX);
}

static void testGibbsSamplingWithSimulatedAnnealing_ConstantTemperature(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, but with constant temperature\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase, constantTemperatureFn, 0, INT64_MAX);
}

static void testGibbsSamplingWithSimulatedAnnealing_WithCooling(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, with exponentially decreasing temperature function\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase,
            exponentiallyDecreasingTemperatureFn, 0, INT64_MAX);
}

static double calculateZScoreSlow(int64_t n, int64_t m, int64_t k, double theta) {
    double score = 0.0;
    for(int64_t i=0; i<n; i++) {
        for(int64_t j=0; j<m; j++) {
            score += pow(1.0 - theta, k + i + j);
        }
    }
    return score;
}

static void testCalculateZScore(CuTest *testCase) {
    for(int64_t i=0; i<100; i++) {
        int64_t n = st_randomInt(0, 100);
        int64_t m = st_randomInt(0, 100);
        int64_t k = st_randomInt(1, 10000);
        double theta = st_random() > 0.05 ? st_random() : 0.0;
        double zScore = calculateZScore(n, m, k, theta);
        double zScoreSlow = calculateZScoreSlow(n, m, k, theta);
        st_logDebug("The slow computed score: %f the fast computed score: %f, n: %" PRIi64 " m: %" PRIi64 " k: %" PRIi64 ", theta: %lf\n", zScoreSlow, zScore, n, m, k, theta);
        CuAssertDblEquals(testCase, zScoreSlow, zScore, 0.000001);
    }
}

CuSuite* referenceProblemTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testMakeReferenceGreedily);
    SUITE_ADD_TEST(suite,
                        testGibbsSamplingWithSimulatedAnnealing_NoExponentiation_Greedy_Fast);
    SUITE_ADD_TEST(suite,
                    testGibbsSamplingWithSimulatedAnnealing_NoExponentiation_Greedy);
    SUITE_ADD_TEST(suite,
                testGibbsSamplingWithSimulatedAnnealing_NoExponentiation);
    SUITE_ADD_TEST(suite,
            testGibbsSamplingWithSimulatedAnnealing_ConstantTemperature);
    SUITE_ADD_TEST(suite, testGibbsSamplingWithSimulatedAnnealing_WithCooling);
    SUITE_ADD_TEST(suite, testCalculateZScore);
    return suite;
}
