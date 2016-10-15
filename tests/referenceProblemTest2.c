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
#include "stReferenceProblem2.h"
#include "stReferenceProblem.h"

static refAdjList *aL = NULL;
static refAdjList *dAL = NULL;
static reference *ref;
static int64_t nodeNumber;
static int64_t intervalNumber;
static int64_t weightNumber;
static int64_t testNumber = 20;

static void teardown() {
    if (aL != NULL) {
        refAdjList_destruct(aL);
        refAdjList_destruct(dAL);
        reference_destruct(ref);
        aL = NULL;
    }
}

int64_t getRandomNode(int64_t nodeNumber) {
    assert(nodeNumber > 0);
    while (1) {
        int64_t n = st_randomInt(-nodeNumber, nodeNumber + 1);
        if (n != 0) {
            return n;
        }
    }
    return 0;
}

static void setup() {
    teardown();
    nodeNumber = st_randomInt(1, 100) * 2;
    intervalNumber = nodeNumber > 0 ? (nodeNumber > 2 ? st_randomInt(1, nodeNumber / 2) : 1) : 0;
    weightNumber = nodeNumber >= 1 ? st_randomInt(0, nodeNumber * nodeNumber * 4) : 0;

    aL = refAdjList_construct(nodeNumber);
    dAL = refAdjList_construct(nodeNumber);
    ref = reference_construct(0);
    for (int64_t i = 0; i < intervalNumber; i++) {
        reference_makeNewInterval(ref, 2 * i + 1, 2 * i + 2);
    }

    for (int64_t i = 0; i < weightNumber; i++) {
        int64_t node1 = getRandomNode(nodeNumber);
        int64_t node2 = getRandomNode(nodeNumber);
        double score = st_random();
        refAdjList_addToWeight(aL, node1, node2, score);
        if (score > 0.8 && refAdjList_getWeight(dAL, node1, node2) == 0.0) {
            refAdjList_addToWeight(dAL, node1, node2, score); //0.8 thresholds ensures that nudging algorithms always results in fewer adjacencies.
        }
    }
    st_logDebug(
            "To test the adjacency problem we've created a problem with %" PRIi64 " nodes, %" PRIi64 " intervals and %" PRIi64 " weights\n",
            nodeNumber, intervalNumber, weightNumber);
}

static void checkIsValidReference(CuTest *testCase) {
    bool *nodes = st_calloc(nodeNumber, sizeof(bool));
    for (int64_t i = 0; i < nodeNumber; i++) {
        nodes[i] = 0;
    }
    for (int64_t i = 0; i < intervalNumber; i++) {
        int64_t n = reference_getFirstOfInterval(ref, i);
        int64_t first = reference_getFirst(ref, n);
        int64_t last = reference_getLast(ref, n);
        //CuAssertIntEquals(testCase, 2*i+1, n);
        while (reference_getNext(ref, n) != INT64_MAX) {
            CuAssertTrue(testCase, n <= nodeNumber);
            CuAssertTrue(testCase, n >= -nodeNumber);
            CuAssertTrue(testCase, n != 0);
            CuAssertTrue(testCase, nodes[llabs(n) - 1] == 0);
            CuAssertIntEquals(testCase, first, reference_getFirst(ref, n));
            CuAssertIntEquals(testCase, last, reference_getLast(ref, n));
            nodes[llabs(n) - 1] = 1;
            n = reference_getNext(ref, n);
        }
        CuAssertTrue(testCase, nodes[llabs(n) - 1] == 0);
        nodes[llabs(n) - 1] = 1;
        //CuAssertIntEquals(testCase, 2*i+2, n);
    }
    for (int64_t i = 0; i < nodeNumber; i++) {
        CuAssertTrue(testCase, nodes[i] == 1);
    }
    free(nodes);
}

static void testEdge(CuTest *testCase) {
    refEdge e = refEdge_construct(1, 5.0);
    CuAssertIntEquals(testCase, refEdge_to(&e), 1);
    CuAssertDblEquals(testCase, refEdge_weight(&e), 5.0, 0.0);
}

static void testAdjList(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        CuAssertIntEquals(testCase, refAdjList_getNodeNumber(aL), nodeNumber);
        if (nodeNumber > 0) {
            for (int64_t j = 0; j < 100; j++) {
                int64_t n1 = getRandomNode(nodeNumber), n2 = getRandomNode(nodeNumber);
                float w = st_random();
                refAdjList_setWeight(aL, n1, n2, w);
                CuAssertDblEquals(testCase, refAdjList_getWeight(aL, n1, n2), w, 0.0);
            }
        }
        long double totalWeight = 0.0;
        for (int64_t n1 = -nodeNumber; n1 <= nodeNumber; n1++) {
            if (n1 != 0) {
                for (int64_t n2 = n1; n2 <= nodeNumber; n2++) {
                    if (n2 != 0) {
                        double w = st_random();
                        refAdjList_setWeight(aL, n1, n2, w);
                        totalWeight += w;
                    }
                }
            }
        }
        st_logDebug("The weights were %f %f\n", refAdjList_getMaxPossibleScore(aL), totalWeight);
        CuAssertDblEquals(testCase, refAdjList_getMaxPossibleScore(aL), totalWeight, 0.0001);
        teardown();
    }
}

static void testReference(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        CuAssertIntEquals(testCase, reference_getIntervalNumber(ref), intervalNumber);
        for (int64_t j = 0; j < intervalNumber; j++) {
            CuAssertIntEquals(testCase, reference_getFirstOfInterval(ref, j), j * 2 + 1);
            CuAssertIntEquals(testCase, reference_getNext(ref, j * 2 + 1), j * 2 + 2);
            CuAssertTrue(testCase, reference_getNext(ref, j*2 + 2) == INT64_MAX);
            CuAssertTrue(testCase, reference_getLast(ref, j * 2 + 1) == j * 2 + 2);
            CuAssertTrue(testCase, reference_getPrevious(ref, j*2 + 1) == INT64_MAX);
            CuAssertTrue(testCase, reference_getPrevious(ref, j * 2 + 2) == j * 2 + 1);
            CuAssertTrue(testCase, reference_inGraph(ref, j * 2 + 1));
            CuAssertTrue(testCase, reference_inGraph(ref, j * 2 + 2));
            CuAssertTrue(testCase, reference_cmp(ref, j * 2 + 1, j * 2 + 2) == -1);
            CuAssertTrue(testCase, reference_cmp(ref, j * 2 + 2, j * 2 + 1) == 1);
            CuAssertTrue(testCase, reference_cmp(ref, j * 2 + 2, j * 2 + 2) == 0);
            CuAssertTrue(testCase, reference_getOrientation(ref, (j * 2 + 1)));
            CuAssertTrue(testCase, !reference_getOrientation(ref, -(j * 2 + 1)));
            CuAssertTrue(testCase, reference_getOrientation(ref, (j * 2 + 2)));
            CuAssertTrue(testCase, !reference_getOrientation(ref, -(j * 2 + 2)));
            CuAssertTrue(testCase, reference_getFirst(ref, j * 2 + 1) == j * 2 + 1);
            CuAssertTrue(testCase, reference_getFirst(ref, j * 2 + 2) == j * 2 + 1);
        }
        st_logInfo("The reference for the %" PRIi64 " th test\n", i);
        reference_log(ref);
        teardown();
    }
}

static void fillReference() {
    for (int64_t n = 2 * intervalNumber + 1; n <= nodeNumber; n++) {
        reference_insertNode(ref, n - 1 > 2 * intervalNumber ? n - 1 : 1, n);
    }
}

static void testReferenceRandom(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        time_t startTime = time(NULL);
        fillReference();
        st_logInfo("Random it took %" PRIi64 " seconds, score: %f of possible: %f\n", time(NULL) - startTime,
                getReferenceScore(aL, ref), refAdjList_getMaxPossibleScore(aL));
        st_logInfo("The reference for the %" PRIi64 " th test\n", i);
        reference_log(ref);
        checkIsValidReference(testCase);
        teardown();
    }
}

void testReference_removeIntervals(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        fillReference();
        checkIsValidReference(testCase);
        while(intervalNumber > 0) {
            //Pick a random interval and remove it.
            stSortedSet *intervalsToRemove = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct);
            int64_t intervalToRemove = st_randomInt(0, intervalNumber);
            int64_t node = reference_getFirstOfInterval(ref, intervalToRemove);
            CuAssertTrue(testCase, reference_inGraph(ref, node));
            stSortedSet_insert(intervalsToRemove, stIntTuple_construct1(node));
            //Get nodes to remove
            stList *nodes = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
            while(node != INT64_MAX) {
                stList_append(nodes, stIntTuple_construct1(node));
                node = reference_getNext(ref, node);
            }
            CuAssertTrue(testCase, reference_inGraph(ref, reference_getFirstOfInterval(ref, intervalToRemove)));
            reference_removeIntervals(ref, intervalsToRemove);
            stSortedSet_destruct(intervalsToRemove);
            //Check nodes no longer in the graph
            for(int64_t i=0; i<stList_length(nodes); i++) {
                CuAssertTrue(testCase, !reference_inGraph(ref, stIntTuple_get(stList_get(nodes, i), 0)));
            }
            //Update recording of reference and check validity
            intervalNumber--;
            nodeNumber -= stList_length(nodes);
            //checkIsValidReference(testCase); This doesn't work because it assumes all the nodes in the sequence 1, 2, ..., nodeNumber-1 are in the graph.
            CuAssertIntEquals(testCase, reference_getIntervalNumber(ref), intervalNumber);
        }
        CuAssertIntEquals(testCase, nodeNumber, 0);
        teardown();
    }
}

void testReference_getMaximumNode(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        fillReference();
        int64_t maxNode = reference_getMaximumNode(ref);
        CuAssertIntEquals(testCase, reference_getMaximumNode(ref), nodeNumber);
        CuAssertTrue(testCase, reference_getMaximumNode(ref) >= 0);
        CuAssertTrue(testCase, maxNode <= reference_getMaximumNode(ref));
        checkIsValidReference(testCase);
        teardown();
    }
}

static void testReference_splitInterval(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        fillReference();
        checkIsValidReference(testCase);
        if (reference_getIntervalNumber(ref) > 0) {
            while (st_random() > 0.1) {
                int64_t n = reference_getFirstOfInterval(ref, st_randomInt(0, reference_getIntervalNumber(ref)));
                while (reference_getNext(ref, n) != INT64_MAX && st_random() > 0.5) {
                    n = reference_getNext(ref, n);
                }
                int64_t m = reference_getNext(ref, n);
                if (m != INT64_MAX) {
                    int64_t i = reference_getIntervalNumber(ref);
                    int64_t first = reference_getFirst(ref, m);
                    assert(reference_getFirst(ref, m) == reference_getFirst(ref, n));
                    int64_t last = reference_getLast(ref, m);
                    assert(reference_getLast(ref, m) == reference_getLast(ref, n));
                    int64_t stub1 = nodeNumber + 1, stub2 = nodeNumber + 2;
                    reference_splitInterval(ref, n, stub1, stub2);
                    CuAssertIntEquals(testCase, i + 1, reference_getIntervalNumber(ref));
                    CuAssertIntEquals(testCase, reference_getNext(ref, n), stub1);
                    CuAssertIntEquals(testCase, reference_getPrevious(ref, stub1), n);
                    CuAssertIntEquals(testCase, reference_getNext(ref, stub2), m);
                    CuAssertIntEquals(testCase, reference_getPrevious(ref, m), stub2);

                    CuAssertIntEquals(testCase, reference_getFirst(ref, stub1), first);
                    CuAssertIntEquals(testCase, reference_getFirst(ref, n), first);
                    CuAssertIntEquals(testCase, reference_getLast(ref, stub1), stub1);
                    CuAssertIntEquals(testCase, reference_getLast(ref, n), stub1);

                    CuAssertIntEquals(testCase, reference_getFirst(ref, stub2), stub2);
                    CuAssertIntEquals(testCase, reference_getFirst(ref, m), stub2);
                    CuAssertIntEquals(testCase, reference_getLast(ref, stub2), last);
                    CuAssertIntEquals(testCase, reference_getLast(ref, m), last);
                    nodeNumber += 2;
                    intervalNumber++;
                    checkIsValidReference(testCase);
                }
            }
        }
        reference_log(ref);
        checkIsValidReference(testCase);
        teardown();
    }
}

bool toySplitFn(int64_t pNode, reference *ref, void *extraArgs) {
    return st_random() > 0.5;
}

static void testMakeReferenceGreedily(CuTest *testCase) {
    long double maxScore = 0, achievedScore = 0;
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        time_t startTime = time(NULL);
        makeReferenceGreedily2(aL, dAL, ref, 0.99);
        int64_t badAdjacencyCount = getBadAdjacencyCount(dAL, ref);
        st_logInfo(
                "Greedy it took %" PRIi64 " seconds, score: %Lf of possible: %Lf, bad adjacency count: %" PRIi64 "\n",
                time(NULL) - startTime, getReferenceScore(aL, ref), refAdjList_getMaxPossibleScore(aL),
                badAdjacencyCount);
        checkIsValidReference(testCase);
        updateReferenceGreedily(aL, dAL, ref, 10);
        long double greedyPermutationScore = getReferenceScore(aL, ref);
        int64_t badAdjacencyCountGreedyPermutations = getBadAdjacencyCount(dAL, ref);
        st_logInfo(
                "Greedy with update permutations, it took %" PRIi64 " seconds, score: %Lf of possible: %Lf, bad adjacency count: %" PRIi64 "\n",
                time(NULL) - startTime, greedyPermutationScore, refAdjList_getMaxPossibleScore(aL),
                badAdjacencyCountGreedyPermutations);
        checkIsValidReference(testCase);
        reorderReferenceToAvoidBreakpoints(aL, ref);
        long double topologicalReorderedScore = getReferenceScore(aL, ref);
        checkIsValidReference(testCase);
        int64_t topologicalBadAdjacencyCount = getBadAdjacencyCount(dAL, ref);
        st_logInfo(
                "Reordered score, it took %" PRIi64 " seconds, score: %Lf of possible: %Lf, bad adjacency count: %" PRIi64 "\n",
                time(NULL) - startTime, topologicalReorderedScore, refAdjList_getMaxPossibleScore(aL),
                topologicalBadAdjacencyCount);
        CuAssertTrue(testCase, topologicalReorderedScore >= greedyPermutationScore);
        //CuAssertTrue(testCase, getBadAdjacencyCount(dAL, ref) <= badAdjacencyCountGreedyPermutations);
        nudgeGreedily(dAL, aL, ref, 10, 100);
        long double nudgeScore = getReferenceScore(aL, ref);
        checkIsValidReference(testCase);
        st_logInfo(
                "Nudge score, it took %" PRIi64 " seconds, score: %Lf of possible: %Lf, bad adjacency count: %" PRIi64 "\n",
                time(NULL) - startTime, nudgeScore, refAdjList_getMaxPossibleScore(aL), getBadAdjacencyCount(dAL, ref));
        CuAssertTrue(testCase, nudgeScore >= topologicalReorderedScore);
        CuAssertTrue(testCase, getBadAdjacencyCount(aL, ref) <= topologicalBadAdjacencyCount);
        stList *splitNodes = splitReferenceAtIndicatedLocations(ref, toySplitFn, NULL);
        CuAssertTrue(testCase, stList_length(splitNodes) % 2 == 0);
        nodeNumber += stList_length(splitNodes);
        intervalNumber += stList_length(splitNodes) / 2;
        checkIsValidReference(testCase);
        //
        reference_log(ref);
        maxScore += refAdjList_getMaxPossibleScore(aL);
        achievedScore += nudgeScore;
        stList_destruct(splitNodes);
        teardown();
    }
    st_logInfo("Got %Lf of possible %Lf score\n", achievedScore, maxScore);
}

static void fn(double theta, int64_t node1, int64_t node2, int64_t adjacencyLength, int64_t node1Length,
        int64_t node2Length, int64_t degree) {
    double d = degree * calculateZScore(node1Length, node2Length, adjacencyLength, theta);
    assert(node1 != node2);
    refAdjList_setWeight(aL, node1, node2, refAdjList_getWeight(aL, node1, node2) + d);
}

static void testADBDCExample(CuTest *testCase) {
    /*
     * Tests example from paper.
     */
    //Nodes
    int64_t A = 1;
    int64_t AL = 2;
    int64_t C = 2;
    int64_t CL = 2;
    int64_t _5B = 3;
    int64_t _3B = -3;
    int64_t BL = 2;
    int64_t _5D = 4;
    int64_t _3D = -4;
    int64_t DL = 2;
    nodeNumber = 4;
    aL = refAdjList_construct(nodeNumber);
    ref = reference_construct(nodeNumber);

    int64_t adjacencyLength = 1;
    int64_t n = 100;
    double theta = 0.0;

    reference_makeNewInterval(ref, C, A);
    //stList_append(chains, stIntTuple_construct(2, _5B, _3B));
    //stList_append(chains, stIntTuple_construct(2, _5D, _3D));

    fn(theta, A, _5B, 2 * adjacencyLength + DL, AL, BL, n - 1);
    fn(theta, A, _3B, 2 * adjacencyLength + DL, AL, BL, 1);
    fn(theta, A, _5D, adjacencyLength, AL, DL, n);
    fn(theta, A, _3D, 3 * adjacencyLength + DL + BL, AL, DL, n);

    fn(theta, A, C, 4 * adjacencyLength + 2 * DL + BL, AL, CL, n);

    fn(theta, C, _5B, 2 * adjacencyLength + DL, CL, BL, 1);
    fn(theta, C, _3B, 2 * adjacencyLength + DL, CL, BL, n - 1);
    fn(theta, C, _5D, adjacencyLength + DL, CL, DL, n);
    fn(theta, C, _3D, 3 * adjacencyLength + DL + BL, CL, DL, n);

    fn(theta, _3D, _5B, adjacencyLength, DL, BL, n);
    fn(theta, _3D, _3B, adjacencyLength, DL, BL, n);

    makeReferenceGreedily2(aL, aL, ref, 0.99);
    updateReferenceGreedily(aL, aL, ref, 100);
    st_logInfo("Running reference example problem, score: %f of possible: %f\n", getReferenceScore(aL, ref),
            refAdjList_getMaxPossibleScore(aL));
    reference_log(ref);
}

CuSuite* referenceProblem2TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testEdge);
    SUITE_ADD_TEST(suite, testAdjList);
    SUITE_ADD_TEST(suite, testReference);
    SUITE_ADD_TEST(suite, testReferenceRandom);
    SUITE_ADD_TEST(suite, testMakeReferenceGreedily);
    SUITE_ADD_TEST(suite, testReference_splitInterval);
    SUITE_ADD_TEST(suite, testReference_getMaximumNode);
    SUITE_ADD_TEST(suite, testReference_removeIntervals);
    SUITE_ADD_TEST(suite, testADBDCExample);
    return suite;
}
