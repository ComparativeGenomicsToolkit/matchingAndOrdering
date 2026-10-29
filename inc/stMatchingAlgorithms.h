/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * externalAlgorithms.h
 *
 *  Created on: 7 Sep 2010
 *      Author: benedictpaten
 */

#ifndef EXTERNALALGORITHMS_H_
#define EXTERNALALGORITHMS_H_

#include "sonLib.h"

extern const char *MATCHING_EXCEPTION;

stIntTuple *constructEdge(int64_t node1, int64_t node2);

stIntTuple *constructWeightedEdge(int64_t node1, int64_t node2, int64_t weight);

/*
 * Chooses a maximum weight matching between the edges, computed in process by Edmonds' blossom
 * algorithm (see impl/maximumWeightMatching.c). Edges are stIntTuple's of length 3 of the form
 * (node1, node2, weight), where node1 and node2 are distinct indices greater than or equal to zero
 * and less than the total node number, weight is an integer, and there is at most one edge between
 * two nodes. The returned list holds the chosen edges, the input's own tuples, in increasing order
 * of their lesser node. Edges of weight zero or less are never chosen; any perfect matching that
 * pads the result with other pairs is then a maximum weight perfect matching of the clique in
 * which missing edges have weight zero, which is what blossom5 was used to compute. The matching
 * is exact for weights up to INT64_MAX / 6; a connected component (of edges of positive weight)
 * with a heavier edge has its weights divided by up to 7 first, so its matching can be lighter
 * than the maximum by less than 7 per edge.
 */
stList *chooseMatching_edmonds(stList *edges, int64_t nodeNumber);

/*
 * Finds maximal matching (maximum cardinality) matching, using the external matchGraph.py.
 * Same form as chooseMatching_edmonds.
 */
stList *chooseMatching_maximumCardinalityMatching(stList *edges, int64_t nodeNumber);

/*
 * Finds maximum weight matching, using the external matchGraph.py. Same form as
 * chooseMatching_edmonds, and of equal weight, but slower.
 */
stList *chooseMatching_maximumWeightMatching(stList *edges, int64_t nodeNumber);

/*
 * Uses a greedy algorithm to choose the matching, starting with the highest weight pair
 * edges are selected in descending order of weight until no further edges can be added to the matching.
 *
 * Same form as chooseMatching_edmonds.
 */
stList *chooseMatching_greedy(stList *edges, int64_t nodeNumber);

/*
 * Returns number of edges with weight > 0.
 */
int64_t matchingCardinality(stList *matching);

/*
 * Returns sum of weights.
 */
int64_t matchingWeight(stList *matching);

stList *getComponents(stList *edges);

#endif /* EXTERNALALGORITHMS_H_ */
