/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stReferenceProblem.h"
#include "stMatchingAlgorithms.h"
#include "stCheckEdges.h"

static void fn(float *zMatrix, int64_t nodeNumber,  double theta, int64_t node1, int64_t node2,
        int64_t adjacencyLength, int64_t node1Length, int64_t node2Length, int64_t degree) {
    double d = degree * calculateZScore(node1Length, node2Length, adjacencyLength, theta);
    assert(node1 != node2);
    zMatrix[nodeNumber * node1 + node2] += d;
    zMatrix[nodeNumber * node2 + node1] += d;
}

static void testADBDCExample(
        CuTest *testCase) {
    /*
     * Tests example from paper.
     */
    //Nodes
    int64_t A=0;
    int64_t AL = 2;
    int64_t C=1;
    int64_t CL = 2;
    int64_t _5B = 2;
    int64_t _3B = 3;
    int64_t BL = 2;
    int64_t _5D = 4;
    int64_t _3D = 5;
    int64_t DL = 2;
    int64_t nodeNumber = 6;

    int64_t adjacencyLength = 1;
    int64_t n = 100;
    double theta = 0.0;

    stList *stubs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    stList_append(stubs, stIntTuple_construct2(C, A));

    stList *chains = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    stList_append(chains, stIntTuple_construct2(_5B, _3B));
    stList_append(chains, stIntTuple_construct2(_5D, _3D));

    float *zMatrix = st_calloc(nodeNumber*nodeNumber, sizeof(float));

    fn(zMatrix, nodeNumber, theta, A, _5B, 2*adjacencyLength + DL, AL, BL, n-1);
    fn(zMatrix, nodeNumber, theta, A, _3B, 2*adjacencyLength + DL, AL, BL, 1);
    fn(zMatrix, nodeNumber, theta, A, _5D, adjacencyLength, AL, DL, n);
    fn(zMatrix, nodeNumber, theta, A, _3D, 3*adjacencyLength + DL + BL, AL, DL, n);

    fn(zMatrix, nodeNumber, theta, A, C, 4*adjacencyLength + 2*DL + BL, AL, CL, n);

    fn(zMatrix, nodeNumber, theta, C, _5B, 2*adjacencyLength + DL, CL, BL, 1);
    fn(zMatrix, nodeNumber, theta, C, _3B, 2*adjacencyLength + DL, CL, BL, n-1);
    fn(zMatrix, nodeNumber, theta, C, _5D, adjacencyLength + DL, CL, DL, n);
    fn(zMatrix, nodeNumber, theta, C, _3D, 3*adjacencyLength + DL + BL, CL, DL, n);

    fn(zMatrix, nodeNumber, theta, _3D, _5B, adjacencyLength, DL, BL, n);
    fn(zMatrix, nodeNumber, theta, _3D, _3B, adjacencyLength, DL, BL, n);

    double totalScore;
    stList *reference = makeReferenceGreedily(stubs, chains, zMatrix, nodeNumber, &totalScore, INT64_MAX);

    greedyImprovement(reference, chains, zMatrix, 100);
    //(reference, chains, zMatrix, 100, NULL);
    logZScoreOfReference(reference, nodeNumber, zMatrix);

    logReference(reference, nodeNumber, zMatrix, totalScore,
                    "result");
    stList *chosenEdges = convertReferenceToAdjacencyEdges(reference);
    stSortedSet *chosenEdgesSet = stList_getSortedSet(chosenEdges, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    for(int64_t i=0; i<stList_length(chosenEdges); i++) {
        stIntTuple *edge = stList_get(chosenEdges, i);
        st_logDebug("I got the edge %" PRIi64 " %" PRIi64 "\n", stIntTuple_get(edge, 0), stIntTuple_get(edge, 1));
    }

    CuAssertTrue(testCase, stList_length(reference) == 1);
    CuAssertTrue(testCase, stList_length(chosenEdges) == 3);
    CuAssertTrue(testCase, stSortedSet_size(chosenEdgesSet) == 3);
    stSortedSet *config1 = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
    stSortedSet_insert(config1, constructEdge(0, 2));
    stSortedSet_insert(config1, constructEdge(3, 5));
    stSortedSet_insert(config1, constructEdge(1, 4));
    stSortedSet *config2 = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
    stSortedSet_insert(config2, constructEdge(0, 4));
    stSortedSet_insert(config2, constructEdge(2, 5));
    stSortedSet_insert(config2, constructEdge(1, 3));
    CuAssertTrue(testCase, stSortedSet_equals(config1, chosenEdgesSet) || stSortedSet_equals(config2, chosenEdgesSet));
}

CuSuite* referenceProblemExamplesTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testADBDCExample);
    return suite;
}
