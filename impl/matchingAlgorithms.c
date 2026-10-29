/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * adjacencyPairsHash.c
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include "commonC.h"
#include "sonLib.h"
#include "shared.h"

const char *MATCHING_EXCEPTION = "MATCHING_EXCEPTION";

/*
 * Functions to assess matching
 */

int64_t matchingCardinality(stList *matching) {
    /*
     * Returns number of edges with weight > 0.
     */
    int64_t totalCardinality = 0;
    for (int64_t i = 0; i < stList_length(matching); i++) {
        stIntTuple *edge = stList_get(matching, i);
        totalCardinality += stIntTuple_get(edge, 2) > 0 ? 1 : 0;
    }
    return totalCardinality;
}

int64_t matchingWeight(stList *matching) {
    /*
     * Returns sum of weights.
     */
    int64_t totalWeight = 0;
    for(int64_t i=0; i<stList_length(matching); i++) {
        stIntTuple *edge = stList_get(matching, i);
        totalWeight += stIntTuple_get(edge, 2);
    }
    return totalWeight;
}

/*
 * Code to talk to matchGraph.py, which computes maximum weight and maximum cardinality matchings.
 */

static void writeGraph(FILE *fileHandle, stList *edges, int64_t nodeNumber) {
    /*
     * Writes out the adjacencies in the format matchGraph.py reads.
     */
    int64_t edgeNumber = stList_length(edges);
    fprintf(fileHandle, "%" PRIi64 " %" PRIi64 "\n", nodeNumber, edgeNumber);
    for(int64_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        int64_t from =  stIntTuple_get(edge, 0);
        int64_t to = stIntTuple_get(edge, 1);
        int64_t weight = stIntTuple_get(edge, 2);
        fprintf(fileHandle, "%" PRIi64 " %" PRIi64 " %" PRIi64 "\n", from, to, weight);
    }
}

static stHash *putEdgesInHash(stList *edges) {
    stHash *intsToEdgesHash = stHash_construct3((uint64_t (*)(const void *))stIntTuple_hashKey, (int (*)(const void *, const void *))stIntTuple_equalsFn, (void (*)(void *))stIntTuple_destruct, NULL);
    for(int64_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        stHash_insert(intsToEdgesHash, constructEdge(stIntTuple_get(edge, 0), stIntTuple_get(edge, 1)), edge);
    }
    return intsToEdgesHash;
}

static stList *readMatching(FILE *fileHandle, stList *originalEdges) {
    /*
     * Reads the matching written by matchGraph.py.
     */
    stHash *originalEdgesHash = putEdgesInHash(originalEdges);
    char *line = stFile_getLineFromFile(fileHandle);
    assert(line != NULL);
    int64_t nodeNumber, edgeNumber;
    int64_t i = sscanf(line, "%" PRIi64 " %" PRIi64 "\n", &nodeNumber, &edgeNumber);
    assert(i == 2);
    free(line);
    stList *chosenEdges = stList_construct();
    for(int64_t j=0; j<edgeNumber; j++) {
        line = stFile_getLineFromFile(fileHandle);
        int64_t node1, node2;
        i = sscanf(line, "%" PRIi64 " %" PRIi64 "", &node1, &node2);
        assert(i == 2);
        free(line);
        assert(node1 >= 0);
        assert(node1 < nodeNumber);
        assert(node2 >= 0);
        assert(node2 < nodeNumber);
        stIntTuple *edge = constructEdge(node1, node2);
        stIntTuple *originalEdge = stHash_search(originalEdgesHash, edge);
        if(originalEdge != NULL) {
            stList_append(chosenEdges, originalEdge);
        }
        stIntTuple_destruct(edge);
    }
    stHash_destruct(originalEdgesHash);
    return chosenEdges;
}

static stList *chooseAdjacencyPairing_externalProgram(stList *edges, int64_t nodeNumber, const char *programName) {
    /*
     * We create temp files to hold stuff.
     */
    if(nodeNumber <= 1) {
        assert(stList_length(edges) == 0);
        return stList_construct();
    }

    char *tempInputFile = getTempFile(), *tempOutputFile = getTempFile();

    /*
     * We write the graph to a temp file.
     */
    FILE *fileHandle = fopen(tempInputFile, "w");
    writeGraph(fileHandle, edges, nodeNumber);
    fclose(fileHandle);

    /*
     * We run the external program.
     */
    char *command = stString_print("%s -e %s -w %s > /dev/null", programName, tempInputFile, tempOutputFile);
    int64_t i = st_system(command);
    if(i != 0) {
        st_errAbort("Something went wrong with the command: %s", command);
        //For some reason this causes a seg fault
        //stThrowNew(MATCHING_EXCEPTION, "Something went wrong with the command: %s", command);
    }
    free(command);

    /*
     * We get back the matching.
     */
    fileHandle = fopen(tempOutputFile, "r");
    stList *matching = readMatching(fileHandle, edges);
    fclose(fileHandle);
    st_logDebug("The adjacency matching for %" PRIi64 " nodes with %" PRIi64 " initial edges contains %" PRIi64 " edges\n", nodeNumber, stList_length(edges), stList_length(matching));

    /*
     * Get rid of the temp files..
     */
    st_system("rm -rf %s %s", tempInputFile, tempOutputFile);
    free(tempInputFile);
    free(tempOutputFile);

    return matching;
}

stList *chooseMatching_maximumCardinalityMatching(stList *edges, int64_t nodeNumber) {
    return chooseAdjacencyPairing_externalProgram(edges, nodeNumber, "matchGraph.py -c");
}

stList *chooseMatching_maximumWeightMatching(stList *edges, int64_t nodeNumber) {
    return chooseAdjacencyPairing_externalProgram(edges, nodeNumber, "matchGraph.py");
}

/*
 * Greedy matching algorithm.
 */

int chooseMatching_greedyP(const void *a, const void *b) {
    return stIntTuple_get((stIntTuple *)a, 2) - stIntTuple_get((stIntTuple *)b, 2);
}

stList *chooseMatching_greedy(stList *edges, int64_t nodeNumber) {
    /*
     * Greedily picks the edge from the list such that each node has at most one edge.
     */
    //First clone the list..
    edges = stList_copy(edges, NULL);

    stSortedSet *seen = getEmptyNodeOrEdgeSetWithCleanup();
    stList *matching = stList_construct();

    //Sort the adjacency pairs..
    stList_sort(edges, chooseMatching_greedyP);

    double strength = INT64_MAX;
    while (stList_length(edges) > 0) {
        stIntTuple *edge = stList_pop(edges);
        double d = stIntTuple_get(edge, 2);
        assert(d <= strength);
        strength = d;
        if(!nodeInSet(seen, stIntTuple_get(edge, 0)) && !nodeInSet(seen, stIntTuple_get(edge, 1))) {
            addNodeToSet(seen, stIntTuple_get(edge, 0));
            addNodeToSet(seen, stIntTuple_get(edge, 1));
            stList_append(matching,edge);
        }
    }
    assert(stList_length(edges) == 0);
    stList_destruct(edges);
    stSortedSet_destruct(seen);

    return matching;
}
