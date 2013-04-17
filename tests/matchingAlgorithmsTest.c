/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stMatchingAlgorithms.h"
#include "shared.h"

static int64_t nodeNumber = 0;
static stSortedSet *edges = NULL;
static stList *edgesList = NULL;

static void teardown() {
    //Gets rid of the random flower.
    if(edges != NULL) {
        stSortedSet_destruct(edges);
        assert(edgesList != NULL);
        stList_destruct(edgesList);
        edgesList = NULL;
        edges = NULL;
        nodeNumber = 0;
    }
}

static void setup() {
    teardown();
    do {
        nodeNumber = st_randomInt(0, 100);
    } while(nodeNumber % 2 != 0);
    edgesList = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    if(nodeNumber > 0) {
        int64_t edgeNumber = st_randomInt(0, nodeNumber * 10);
        stSortedSet *seen = getEmptyNodeOrEdgeSetWithCleanup();
        for(int64_t i=0; i<edgeNumber; i++) {
            int64_t from = st_randomInt(0, nodeNumber);
            int64_t to = st_randomInt(0, nodeNumber);
            if(from != to) {
                if(!edgeInSet(seen, from, to)) {
                    addEdgeToSet(seen, from, to);
                    addWeightedEdgeToList(from, to, st_randomInt(0, 100), edgesList);
                }
            }
        }
        stSortedSet_destruct(seen);
    }
    edges = stList_getSortedSet(edgesList, (int (*)(const void *, const void *))stIntTuple_cmpFn);
}

static void checkMatching(CuTest *testCase, stList *matching, bool perfectMatching) {
    /*
     * Checks the matching is valid.
     */
    stSortedSet *seen = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct);
    for(int64_t i=0; i<stList_length(matching); i++) {
        stIntTuple *edge = stList_get(matching, i);
        int64_t from = stIntTuple_get(edge, 0);
        int64_t to = stIntTuple_get(edge, 1);
        int64_t weight = stIntTuple_get(edge, 2);

        /*
         * Check bounds are valid.
         */
        CuAssertTrue(testCase, from != to);
        CuAssertTrue(testCase, from >= 0 && from < nodeNumber);
        CuAssertTrue(testCase, to >= 0 && to < nodeNumber);
        CuAssertTrue(testCase, weight >= 0 && weight < 100);

        /*
         * Check the matching is valid.
         */
        stIntTuple *fromTuple = stIntTuple_construct1( from);
        CuAssertTrue(testCase, stSortedSet_search(seen, fromTuple) == NULL);
        stSortedSet_insert(seen, fromTuple);
        stIntTuple *toTuple = stIntTuple_construct1( to);
        CuAssertTrue(testCase, stSortedSet_search(seen, toTuple) == NULL);
        stSortedSet_insert(seen, toTuple);

        /*
         * Check that any edge with non-zero weight is in our original set.
         */
        if(weight > 0.0) {
            CuAssertTrue(testCase, stSortedSet_search(edges, edge) != NULL);
        }
    }
    if(perfectMatching) {
        CuAssertTrue(testCase, stList_length(matching)*2 == nodeNumber);
        CuAssertTrue(testCase, stSortedSet_size(seen) == nodeNumber);
    }
    stSortedSet_destruct(seen);
}

static void testGreedy(CuTest *testCase) {
    for(int64_t i=0; i<100; i++) {
        setup();
        stList *matching = chooseMatching_greedy(edgesList, nodeNumber);
        checkMatching(testCase, matching, 0);
        int64_t totalWeight = matchingWeight(matching);
        st_logInfo("The total weight of the greedy matching is %" PRIi64 "\n", totalWeight);
        stList_destruct(matching);
        teardown();
    }
}

static void testMaximumWeight(CuTest *testCase) {
    /*
     * Creates random graphs, constructs matchings with the blossom5 and maximumWeight algorithms
     * and checks that they have equal weight and higher or equal weight to the greedy matching
     * as well as sanity checking the matchings (the blossom matching is perfect).
     */
    for(int64_t i=0; i<100; i++) {
        setup();
        stList *greedyMatching = chooseMatching_greedy(edgesList, nodeNumber);
        stList *blossomMatching = chooseMatching_blossom5(edgesList, nodeNumber);
        stList *maximumWeightMatching = chooseMatching_maximumWeightMatching(edgesList, nodeNumber);
        checkMatching(testCase, greedyMatching, 0);
        checkMatching(testCase, blossomMatching, 0);
        checkMatching(testCase, maximumWeightMatching, 0);
        int64_t totalGreedyWeight = matchingWeight(greedyMatching);
        int64_t totalBlossumWeight = matchingWeight(blossomMatching);
        int64_t totalMaximumWeightWeight = matchingWeight(maximumWeightMatching);
        st_logInfo("The total weight of the greedy matching is %" PRIi64 ", the total weight of the blossom5 matching is %" PRIi64 ", the total weight of the maximum weight matching is %" PRIi64 "\n",
                totalGreedyWeight, totalBlossumWeight, totalMaximumWeightWeight);
        st_logInfo("The total cardinality of the greedy matching is %" PRIi64 ", the total cardinality of the blossom5  matching is %" PRIi64 ", the total cardinality of the maximum weight matching is %" PRIi64 "\n",
                        stList_length(greedyMatching), stList_length(blossomMatching), stList_length(maximumWeightMatching));
        CuAssertTrue(testCase, totalGreedyWeight <= totalBlossumWeight);
        CuAssertTrue(testCase, totalGreedyWeight <= totalMaximumWeightWeight);
        CuAssertTrue(testCase, totalBlossumWeight == totalMaximumWeightWeight);
        //CuAssertTrue(testCase, stList_length(greedyMatching) <= stList_length(blossomMatching));
        stList_destruct(greedyMatching);
        stList_destruct(blossomMatching);
        stList_destruct(maximumWeightMatching);
        teardown();
    }
}

static void testMaximumCardinality(CuTest *testCase) {
    /*
     * Tests a maximum (cardinality) matching algorithm, checking that it has higher or equal
     * cardinality to the greedy algorithm.
     */
    for(int64_t i=0; i<100; i++) {
        setup();
        stList *greedyMatching = chooseMatching_greedy(edgesList, nodeNumber);
        stList *edmondsMatching = chooseMatching_maximumCardinalityMatching(edgesList, nodeNumber);
        checkMatching(testCase, greedyMatching, 0);
        checkMatching(testCase, edmondsMatching, 0);

        int64_t totalGreedyWeight = matchingWeight(greedyMatching);
        int64_t totalEdmondsWeight = matchingWeight(edmondsMatching);
        st_logInfo("The total weight of the greedy matching is %" PRIi64 ", the total weight of the edmonds matching is %" PRIi64 "\n",
                totalGreedyWeight, totalEdmondsWeight);
        st_logInfo("The total cardinality of the greedy matching is %" PRIi64 ", the total cardinality of the edmonds matching is %" PRIi64 "\n",
                stList_length(greedyMatching), stList_length(edmondsMatching));
        CuAssertTrue(testCase, stList_length(greedyMatching) <= stList_length(edmondsMatching));
        stList_destruct(greedyMatching);
        stList_destruct(edmondsMatching);
        teardown();
    }
}

CuSuite* matchingAlgorithmsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testGreedy);
    SUITE_ADD_TEST(suite, testMaximumWeight);
    SUITE_ADD_TEST(suite, testMaximumCardinality);

    return suite;
}
