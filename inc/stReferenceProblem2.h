/*
 * referenceProblem.h
 *
 *  Created on: 15 Jan 2013
 *      Author: benedictpaten
 */

#ifndef REFERENCEPROBLEM2_H_
#define REFERENCEPROBLEM2_H_

#include "sonLib.h"

typedef struct _refEdge refEdge;

typedef struct _refAdjList refAdjList;

typedef struct _refAdjListIt refAdjListIt;

typedef struct _reference refOrdering;

struct _refEdge {
    int64_t to;
    double weight;
};

struct _refAdjListIt {
    stHash *hash;
    stHashIterator *it;
};

/*
 * Adjacency list structure/edge structure
 */
//Negative value indicates the 3' end of segment, positive value indicates 5 prime end.

refEdge refEdge_construct(int64_t to, double weight);

int64_t refEdge_to(refEdge *e);

double refEdge_weight(refEdge *e);

refAdjList *refAdjList_construct(int64_t nodeNumber);

void refAdjList_destruct(refAdjList *aL);

int64_t refAdjList_getNodeNumber(refAdjList *aL);

double refAdjList_getWeight(refAdjList *aL, int64_t n1, int64_t n2);

void refAdjList_setWeight(refAdjList *aL, int64_t n1, int64_t n2, double weight);

void refAdjList_addToWeight(refAdjList *aL, int64_t n1, int64_t n2, double weight);

refAdjListIt adjList_getEdgeIt(refAdjList *aL, int64_t node);

refEdge refAdjListIt_getNext(refAdjListIt *it);

void refAdjListIt_destruct(refAdjListIt *it);

long double refAdjList_getMaxPossibleScore(refAdjList *aL);

int64_t refAdjList_getNumberOfWeights(refAdjList *aL);

//double calculateZScore(int64_t n, int64_t m, int64_t k, double theta);

/*
 * Reference structure
 */

refOrdering *reference_construct(int64_t nodeNumber);

void reference_destruct(refOrdering *ref);

//Need to make one or more intervals before you can insert other nodes into reference.
void reference_makeNewInterval(refOrdering *ref, int64_t leftNode, int64_t rightNode);

void reference_insertNode(refOrdering *ref, int64_t pNode, int64_t node);

bool reference_inGraph(refOrdering *ref, int64_t n);

int64_t reference_getFirstOfInterval(refOrdering *ref, int64_t interval);

int64_t reference_getIntervalNumber(refOrdering *ref);

int64_t reference_getFirst(refOrdering *ref, int64_t n);

int64_t reference_getPrevious(refOrdering *ref, int64_t n);

int64_t reference_getLast(refOrdering *ref, int64_t n);

//Returns nonzero if segment is in reference in same orientation, traversing the reference sequence(s) from 5' to 3'.
bool reference_getOrientation(refOrdering *ref, int64_t n);

//Returns true if the edge (m ,n) is consistent with the reference, where m and n are sides.
bool reference_isConsistent(refOrdering *ref, int64_t m, int64_t n);

//Gets the next position within the reference.
int64_t reference_getNext(refOrdering *ref, int64_t n);

int64_t reference_getRemainingIntervalLength(refOrdering *ref, int64_t n);

//Compares segments position within a reference, ignoring orientation.
int reference_cmp(refOrdering *ref, int64_t n1, int64_t n2);

void reference_log(refOrdering *ref);

//Splits an interval into two, making the existing interval end with stub1, and the new interval (which will be last interval) start with stub2.
void reference_splitInterval(refOrdering *ref, int64_t pNode, int64_t stub1, int64_t stub2);

//Gets the maximum value of a node in the reference. Makes it easy to generate a unique new node.
int64_t reference_getMaximumNode(refOrdering *ref);

//Remove intervals, input is set of first nodes of intervals. The intervals containing these first nodes are removed, including all members of the interval.
void reference_removeIntervals(refOrdering *ref, stSortedSet *firstNodesOfIntervalsToRemove);

/*
 * Reference algorithms
 */

void makeReferenceGreedily2(refAdjList *aL, refAdjList *dAL, refOrdering *ref, double wiggle);

void updateReferenceGreedily(refAdjList *aL, refAdjList *dAL, refOrdering *ref, int64_t permutations);

/*
 * Create a topological sort of each reference interval, trying to place nodes that are connected by direct adjacencies next to one another.
 */
void reorderReferenceToAvoidBreakpoints(refAdjList *aL, refOrdering *ref);

long double getReferenceScore(refAdjList *aL, refOrdering *ref);

/*
 * Nudge the blocks to try to eliminate "bad adjacencies", where to blocks are adjacent but do not have a direct weight between them.
 */
void nudgeGreedily(refAdjList *dAL, refAdjList *aL, refOrdering *ref, int64_t permutations, int64_t maxNudge);

/*
 * Count of adjacent nodes in reference that have no edge connecting them.
 */
int64_t getBadAdjacencyCount(refAdjList *aL, refOrdering *ref);

/*
 * Breaks up chromosomes.
 */

/*
 * Splits the reference up at any adjacency that refSplitFn returns non-zero for.
 */
stList *splitReferenceAtIndicatedLocations(refOrdering *ref, bool (*refSplitFn)(int64_t, refOrdering *, void *), void *extraArgs);

/*
 * Rejoins reference intervals that are specified by referenceIntervalsToPreserve. Returns a modified list of extraStubNodes with
 * unneeded stubs removed.
 */
stList *remakeReferenceIntervals(refOrdering *ref, stList *referenceIntervalsToPreserve, stList *extraStubNodes);

double exponentiallyDecreasingTemperatureFn(double d);

double constantTemperatureFn(double d);

long double calculateZScore(int64_t n, int64_t m, int64_t k, long double theta);

#endif /* REFERENCEPROBLEM2_H_ */
