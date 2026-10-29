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
     * Creates random graphs, constructs matchings with the edmonds and maximumWeight algorithms
     * and checks that they have equal weight and higher or equal weight to the greedy matching
     * as well as sanity checking the matchings.
     */
    for(int64_t i=0; i<100; i++) {
        setup();
        stList *greedyMatching = chooseMatching_greedy(edgesList, nodeNumber);
        stList *edmondsMatching = chooseMatching_edmonds(edgesList, nodeNumber);
        stList *maximumWeightMatching = chooseMatching_maximumWeightMatching(edgesList, nodeNumber);
        checkMatching(testCase, greedyMatching, 0);
        checkMatching(testCase, edmondsMatching, 0);
        checkMatching(testCase, maximumWeightMatching, 0);
        int64_t totalGreedyWeight = matchingWeight(greedyMatching);
        int64_t totalEdmondsWeight = matchingWeight(edmondsMatching);
        int64_t totalMaximumWeightWeight = matchingWeight(maximumWeightMatching);
        st_logInfo("The total weight of the greedy matching is %" PRIi64 ", the total weight of the edmonds matching is %" PRIi64 ", the total weight of the maximum weight matching is %" PRIi64 "\n",
                totalGreedyWeight, totalEdmondsWeight, totalMaximumWeightWeight);
        st_logInfo("The total cardinality of the greedy matching is %" PRIi64 ", the total cardinality of the edmonds matching is %" PRIi64 ", the total cardinality of the maximum weight matching is %" PRIi64 "\n",
                        stList_length(greedyMatching), stList_length(edmondsMatching), stList_length(maximumWeightMatching));
        CuAssertTrue(testCase, totalGreedyWeight <= totalEdmondsWeight);
        CuAssertTrue(testCase, totalGreedyWeight <= totalMaximumWeightWeight);
        CuAssertTrue(testCase, totalEdmondsWeight == totalMaximumWeightWeight);
        stList_destruct(greedyMatching);
        stList_destruct(edmondsMatching);
        stList_destruct(maximumWeightMatching);
        teardown();
    }
}

/*
 * Tests of chooseMatching_edmonds on its own.
 */

static int64_t bruteForceMatchingWeight(stList *edges, int64_t nodeNumber) {
    /*
     * The weight of a maximum weight matching, by exhaustive search: best[mask] is the weight of a
     * maximum weight matching among the nodes in mask.
     */
    int64_t *weights = st_calloc(nodeNumber * nodeNumber + 1, sizeof(int64_t));
    for(int64_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        int64_t x = stIntTuple_get(edge, 0), y = stIntTuple_get(edge, 1), w = stIntTuple_get(edge, 2);
        if(w > weights[x * nodeNumber + y]) {
            weights[x * nodeNumber + y] = weights[y * nodeNumber + x] = w;
        }
    }
    int64_t *best = st_calloc(((int64_t)1) << nodeNumber, sizeof(int64_t));
    for(int64_t mask=1; mask < (((int64_t)1) << nodeNumber); mask++) {
        int64_t i = 0;
        while(!(mask >> i & 1)) {
            i++;
        }
        int64_t rest = mask & ~(((int64_t)1) << i);
        int64_t b = best[rest]; //Leave node i unmatched
        for(int64_t j=i+1; j<nodeNumber; j++) {
            int64_t w = weights[i * nodeNumber + j];
            if((rest >> j & 1) && w > 0 && w + best[rest & ~(((int64_t)1) << j)] > b) {
                b = w + best[rest & ~(((int64_t)1) << j)];
            }
        }
        best[mask] = b;
    }
    int64_t weight = best[(((int64_t)1) << nodeNumber) - 1];
    free(weights);
    free(best);
    return weight;
}

static void checkEdmondsMatching(CuTest *testCase, stList *edges, int64_t nodeNumber, stList *matching) {
    /*
     * Checks the chosen edges are input edges of positive weight, that no node is used twice and
     * that they are listed in increasing order of their lesser node.
     */
    stSet *inputEdges = stList_getSet(edges);
    bool *used = st_calloc(nodeNumber + 1, sizeof(bool));
    int64_t previousNode = -1;
    for(int64_t i=0; i<stList_length(matching); i++) {
        stIntTuple *edge = stList_get(matching, i);
        int64_t x = stIntTuple_get(edge, 0), y = stIntTuple_get(edge, 1);
        CuAssertTrue(testCase, stSet_search(inputEdges, edge) != NULL);
        CuAssertTrue(testCase, stIntTuple_get(edge, 2) > 0);
        CuAssertTrue(testCase, x != y && x >= 0 && x < nodeNumber && y >= 0 && y < nodeNumber);
        CuAssertTrue(testCase, !used[x] && !used[y]);
        used[x] = used[y] = 1;
        int64_t lesserNode = x < y ? x : y;
        CuAssertTrue(testCase, lesserNode > previousNode);
        previousNode = lesserNode;
    }
    free(used);
    stSet_destruct(inputEdges);
}

/*
 * Graphs whose matchings make, relabel, expand and augment through blossoms in the ways that are
 * hard to get right, with the weights of their maximum weight matchings. They come from the unit
 * tests of Joris van Rantwijk's maximum-weight-matching (MIT license), which
 * impl/maximumWeightMatching.c is translated from.
 */
typedef struct _knownGraph {
    int64_t weight;
    int64_t edgeNumber;
    int64_t edges[3 * 13];
} KnownGraph;

static const KnownGraph knownGraphs[] = {
    { 1, 1, { 0,1,1 } }, //single edge
    { 11, 2, { 1,2,10, 2,3,11 } },
    { 11, 3, { 1,2,5, 2,3,11, 3,4,5 } },
    { 2, 5, { 1,2,2, 1,3,-2, 2,3,1, 2,4,-1, 3,4,-6 } }, //negative weights
    { 15, 4, { 1,2,8, 1,3,9, 2,3,10, 3,4,7 } }, //S-blossom, used to augment
    { 21, 6, { 1,2,8, 1,3,9, 2,3,10, 3,4,7, 1,6,5, 4,5,6 } },
    { 17, 6, { 1,2,9, 1,3,8, 2,3,10, 1,4,5, 4,5,4, 1,6,3 } }, //S-blossom relabelled T, used to augment
    { 17, 6, { 1,2,9, 1,3,8, 2,3,10, 1,4,5, 4,5,3, 1,6,4 } },
    { 16, 6, { 1,2,9, 1,3,8, 2,3,10, 1,4,5, 4,5,3, 3,6,4 } },
    { 23, 7, { 1,2,9, 1,3,9, 2,3,10, 2,4,8, 3,5,8, 4,5,10, 5,6,6 } }, //nested S-blossom, used to augment
    { 48, 9, { 1,2,10, 1,7,10, 2,3,12, 3,4,20, 3,5,20, 4,5,25, 5,6,10, 6,7,10, 7,8,8 } }, //S-blossom relabelled S, nested
    { 44, 10, { 1,2,8, 1,3,8, 2,3,10, 2,4,12, 3,5,12, 4,5,14, 4,6,12, 5,7,12, 6,7,14, 7,8,12 } }, //nested S-blossom, augment, expand recursively
    { 67, 8, { 1,2,23, 1,5,22, 1,6,15, 2,3,25, 3,4,22, 4,5,25, 4,8,14, 5,7,13 } }, //S-blossom relabelled T, expanded
    { 47, 9, { 1,2,19, 1,3,20, 1,8,8, 2,3,25, 2,4,18, 3,5,18, 4,5,13, 4,7,7, 5,6,7 } }, //nested S-blossom relabelled T, expanded
    { 146, 10, { 1,2,45, 1,5,45, 2,3,50, 3,4,45, 4,5,50, 1,6,30, 3,9,35, 4,8,35, 5,7,26, 9,10,5 } }, //blossom relabelled T in more than one way, expanded, augment
    { 151, 10, { 1,2,45, 1,5,45, 2,3,50, 3,4,45, 4,5,50, 1,6,30, 3,9,35, 4,8,26, 5,7,40, 9,10,5 } },
    { 139, 10, { 1,2,45, 1,5,45, 2,3,50, 3,4,45, 4,5,50, 1,6,30, 3,9,35, 4,8,28, 5,7,26, 9,10,5 } }, //expansion makes a new least slack edge
    { 241, 13, { 1,2,45, 1,7,45, 2,3,50, 3,4,45, 4,5,95, 4,6,94, 5,6,94, 6,7,50, 1,8,30, 3,11,35, 5,9,36, 7,10,26, 11,12,5 } }, //expansion puts an inner blossom on the augmenting path
    { 145, 11, { 1,2,40, 1,3,40, 2,3,60, 2,4,55, 3,5,55, 4,5,50, 1,8,15, 5,7,30, 7,6,10, 8,10,10, 4,9,30 } }, //nested S-blossom relabelled S, expanded recursively
    { 15, 8, { 0,1,2, 0,4,3, 1,2,7, 1,5,2, 2,3,9, 2,5,4, 3,4,8, 3,5,4 } }, //nodes left unmatched
    { 20, 7, { 0,1,8, 0,2,8, 1,2,9, 1,3,6, 2,4,7, 3,4,8, 3,5,5 } }, //nested S-blossom, only the outer expanded
    { 20, 7, { 0,1,7, 0,2,7, 1,2,9, 0,3,7, 0,4,7, 3,4,9, 5,6,2 } }, //nested S-blossom left with an inner zero dual
    { 18, 7, { 0,1,7, 0,4,6, 1,2,9, 2,3,8, 3,4,9, 2,5,1, 4,5,1 } }, //redundant edge dropped when making a blossom
    { 8, 5, { 1,3,7, 1,4,3, 2,4,1, 3,4,5, 3,5,5 } }, //blossom expanded before it is labelled T
    { 34, 8, { 1,2,15, 1,3,10, 1,4,11, 1,6,17, 2,5,12, 2,6,8, 5,7,15, 6,7,7 } },
    { 37, 8, { 1,2,12, 1,4,11, 1,5,11, 2,4,12, 2,6,11, 3,4,14, 3,7,14, 4,7,11 } },
    { 52, 8, { 1,2,19, 1,4,17, 1,5,19, 2,3,15, 2,5,21, 4,6,18, 4,7,11, 5,6,19 } }, //nested blossoms with a common base augmented
    { 46, 7, { 1,2,19, 1,3,19, 1,4,15, 1,5,17, 2,3,21, 4,5,16, 5,6,10 } },
    { 4, 11, { 1,2,1, 1,3,1, 2,3,1, 4,5,1, 4,6,1, 5,6,1, 7,8,1, 7,9,1, 8,9,1, 1,4,1, 5,8,1 } } //three triangles
};

static void testEdmondsKnownGraphs(CuTest *testCase) {
    /*
     * Checks the weight of the matching of each known graph, with its edges in order, in reverse
     * order and with their ends swapped, and against brute force.
     */
    for(int64_t i=0; i<(int64_t)(sizeof(knownGraphs) / sizeof(KnownGraph)); i++) {
        const KnownGraph *graph = &knownGraphs[i];
        int64_t graphNodeNumber = 0;
        for(int64_t j=0; j<graph->edgeNumber; j++) {
            for(int64_t k=0; k<2; k++) {
                if(graph->edges[3 * j + k] + 1 > graphNodeNumber) {
                    graphNodeNumber = graph->edges[3 * j + k] + 1;
                }
            }
        }
        for(int64_t variant=0; variant<3; variant++) {
            stList *graphEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
            for(int64_t j=0; j<graph->edgeNumber; j++) {
                const int64_t *e = &graph->edges[3 * (variant == 1 ? graph->edgeNumber - 1 - j : j)];
                stList_append(graphEdges, variant == 2 ? stIntTuple_construct3(e[1], e[0], e[2]) : stIntTuple_construct3(e[0], e[1], e[2]));
            }
            stList *matching = chooseMatching_edmonds(graphEdges, graphNodeNumber);
            checkEdmondsMatching(testCase, graphEdges, graphNodeNumber, matching);
            CuAssertIntEquals(testCase, graph->weight, matchingWeight(matching));
            CuAssertIntEquals(testCase, graph->weight, bruteForceMatchingWeight(graphEdges, graphNodeNumber));
            stList_destruct(matching);
            stList_destruct(graphEdges);
        }
    }
}

static void testEdmondsBruteForce(CuTest *testCase) {
    /*
     * Compares the weight of the matching with brute force on random graphs of up to 14 nodes,
     * odd and even, with weights that tie often, rarely, and that come close to the INT64_MAX / 6
     * limit (blossom5 got weights of 2^30 or more wrong).
     */
    for(int64_t i=0; i<2000; i++) {
        int64_t weightKind = i % 4;
        int64_t graphNodeNumber = st_randomInt(0, weightKind == 3 ? 9 : 15);
        double edgeProbability = st_random();
        stList *graphEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
        for(int64_t x=0; x<graphNodeNumber; x++) {
            for(int64_t y=x+1; y<graphNodeNumber; y++) {
                if(st_random() < edgeProbability) {
                    int64_t w = weightKind == 0 ? st_randomInt(0, 2) : weightKind == 1 ? st_randomInt(0, 11) :
                            weightKind == 2 ? st_randomInt(-10, 1000) : INT64_MAX / 6 - st_randomInt(0, 1000000) * (INT64_MAX / 6000000);
                    stList_append(graphEdges, st_random() < 0.5 ? stIntTuple_construct3(x, y, w) : stIntTuple_construct3(y, x, w));
                }
            }
        }
        stList_shuffle(graphEdges);
        stList *matching = chooseMatching_edmonds(graphEdges, graphNodeNumber);
        checkEdmondsMatching(testCase, graphEdges, graphNodeNumber, matching);
        CuAssertTrue(testCase, matchingWeight(matching) == bruteForceMatchingWeight(graphEdges, graphNodeNumber));
        stList_destruct(matching);
        stList_destruct(graphEdges);
    }
}

static void testEdmondsDisjointEdges(CuTest *testCase) {
    /*
     * The graphs cactus gives it: an edge between the two stub ends of each thread, of any weight
     * up to INT64_MAX, with the other pairs of weight zero. All of the positive edges are chosen.
     */
    for(int64_t i=0; i<20; i++) {
        int64_t graphNodeNumber = 2 * st_randomInt(0, 500);
        stList *graphEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
        int64_t positiveEdgeNumber = 0;
        for(int64_t x=0; x<graphNodeNumber; x+=2) {
            int64_t w = st_random() < 0.1 ? 0 : INT64_MAX - st_randomInt(0, 1000000) * (INT64_MAX / 1000000);
            positiveEdgeNumber += w > 0 ? 1 : 0;
            stList_append(graphEdges, stIntTuple_construct3(x + 1, x, w));
            if(x > 0) {
                stList_append(graphEdges, stIntTuple_construct3(x - 1, x, 0));
            }
        }
        stList_shuffle(graphEdges);
        stList *matching = chooseMatching_edmonds(graphEdges, graphNodeNumber);
        checkEdmondsMatching(testCase, graphEdges, graphNodeNumber, matching);
        CuAssertIntEquals(testCase, positiveEdgeNumber, stList_length(matching));
        for(int64_t j=0; j<stList_length(matching); j++) {
            stIntTuple *edge = stList_get(matching, j);
            CuAssertTrue(testCase, stIntTuple_get(edge, 0) / 2 == stIntTuple_get(edge, 1) / 2);
        }
        stList_destruct(matching);
        stList_destruct(graphEdges);
    }
}

static void testEdmondsHeavyWeights(CuTest *testCase) {
    /*
     * Weights above INT64_MAX / 6 are divided down to be matched, which can lose a little weight,
     * but not more than 7 per edge.
     */
    stList *graphEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    int64_t w = ((int64_t)1) << 61;
    stList_append(graphEdges, stIntTuple_construct3(0, 1, w));
    stList_append(graphEdges, stIntTuple_construct3(1, 2, 2 * w + 1));
    stList_append(graphEdges, stIntTuple_construct3(2, 3, w));
    stList *matching = chooseMatching_edmonds(graphEdges, 4);
    CuAssertIntEquals(testCase, 1, stList_length(matching));
    CuAssertTrue(testCase, stList_get(matching, 0) == stList_get(graphEdges, 1));
    stList_destruct(matching);
    stList_destruct(graphEdges);
    for(int64_t i=0; i<500; i++) {
        int64_t graphNodeNumber = st_randomInt(0, 7);
        graphEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
        for(int64_t x=0; x<graphNodeNumber; x++) {
            for(int64_t y=x+1; y<graphNodeNumber; y++) {
                if(st_random() < 0.7) {
                    stList_append(graphEdges, stIntTuple_construct3(x, y, INT64_MAX / 4 - st_randomInt(0, 1000)));
                }
            }
        }
        matching = chooseMatching_edmonds(graphEdges, graphNodeNumber);
        checkEdmondsMatching(testCase, graphEdges, graphNodeNumber, matching);
        int64_t bestWeight = bruteForceMatchingWeight(graphEdges, graphNodeNumber);
        CuAssertTrue(testCase, matchingWeight(matching) <= bestWeight);
        CuAssertTrue(testCase, matchingWeight(matching) > bestWeight - 7 * stList_length(matching) - 1);
        stList_destruct(matching);
        stList_destruct(graphEdges);
    }
}

static void testEdmondsSmallCases(CuTest *testCase) {
    /*
     * No nodes, one node, an odd number of nodes and only zero weight edges.
     */
    stList *graphEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int64_t n=0; n<2; n++) {
        stList *matching = chooseMatching_edmonds(graphEdges, n);
        CuAssertIntEquals(testCase, 0, stList_length(matching));
        stList_destruct(matching);
    }
    stList_append(graphEdges, stIntTuple_construct3(0, 1, 5));
    stList_append(graphEdges, stIntTuple_construct3(1, 2, 6));
    stList_append(graphEdges, stIntTuple_construct3(0, 2, 0));
    stList *matching = chooseMatching_edmonds(graphEdges, 3);
    CuAssertIntEquals(testCase, 1, stList_length(matching));
    CuAssertTrue(testCase, stList_get(matching, 0) == stList_get(graphEdges, 1));
    stList_destruct(matching);
    stList_destruct(graphEdges);
    graphEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    stList_append(graphEdges, stIntTuple_construct3(0, 1, 0));
    stList_append(graphEdges, stIntTuple_construct3(2, 3, 0));
    matching = chooseMatching_edmonds(graphEdges, 4);
    CuAssertIntEquals(testCase, 0, stList_length(matching));
    stList_destruct(matching);
    stList_destruct(graphEdges);
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
    SUITE_ADD_TEST(suite, testEdmondsKnownGraphs);
    SUITE_ADD_TEST(suite, testEdmondsBruteForce);
    SUITE_ADD_TEST(suite, testEdmondsDisjointEdges);
    SUITE_ADD_TEST(suite, testEdmondsHeavyWeights);
    SUITE_ADD_TEST(suite, testEdmondsSmallCases);

    return suite;
}
