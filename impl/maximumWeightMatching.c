/*
 * maximumWeightMatching.c
 *
 * A maximum weight matching in a general graph with integer edge weights, by Edmonds' primal-dual
 * blossom algorithm as described by Galil, "Efficient algorithms for finding maximum matching in
 * graphs", ACM Computing Surveys 18 (1986), with lazy updates of the dual variables and reuse of
 * the alternating trees between augmentations, which bring it to O(n m log(n)) time.
 *
 * This is a C translation of the C++ implementation in Joris van Rantwijk's maximum-weight-matching
 * (https://git.jorisvr.nl/joris/maximum-weight-matching, tag v3: cpp/mwmatching.hpp,
 * cpp/concatenable_queue.hpp and cpp/priority_queue.hpp), which is released under the MIT license:
 *
 * Copyright (c) 2023-2024 Joris van Rantwijk
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 * and associated documentation files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next paragraph) shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 * BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <stdlib.h>
#include <string.h>
#include "sonLib.h"
#include "stMatchingAlgorithms.h"

/*
 * The algorithm in brief. Every vertex starts as the root of its own alternating tree, labelled S.
 * Each stage grows the trees along tight edges (edges of zero slack in the dual of the matching
 * LP), shrinking odd cycles into blossoms, until a tight edge joins two different trees, which
 * gives an augmenting path; the matching is augmented, the two trees are dissolved and the others
 * are kept for the next stage. When no tight edge can be used, the dual variables are moved by the
 * largest step (delta) that keeps them feasible. Of the four kinds of step, the first makes the
 * dual of the unmatched vertices zero, which proves the matching optimal and ends the algorithm;
 * the second and third make an edge out of an S-blossom tight; the fourth makes the dual of a
 * T-blossom zero so that it can be expanded. Delta steps are applied lazily, through a running sum
 * and per blossom offsets, and the least slack edges and blossoms are kept in priority queues.
 *
 * As in the C++ code, dual variables and slacks are kept at twice their value so that they stay
 * integers, and the edge weights must be at most a sixth of INT64_MAX to rule out overflow.
 */

#define NO_VERTEX -1
#define NO_EDGE -1
#define NO_BLOSSOM -1
#define INVALID_SLACK INT64_MAX

#define LABEL_NONE 0
#define LABEL_S 1
#define LABEL_T 2

////////////////////////////////////
// Priority queue
////////////////////////////////////

/*
 * A min-priority queue implemented as a binary heap of nodes that live in the caller's own
 * structures, so that an element can be found, re-prioritised or removed in O(log(n)).
 */

typedef struct _heapNode {
    int64_t index; //position in the heap, or -1 if the node is not in a heap
    int64_t prio;
    int64_t data;
} HeapNode;

typedef struct _heap {
    HeapNode **heap;
    int64_t size;
    int64_t capacity;
} Heap;

static void heap_init(Heap *q) {
    q->heap = NULL;
    q->size = 0;
    q->capacity = 0;
}

static void heap_destruct(Heap *q) {
    free(q->heap);
}

static void heapNode_init(HeapNode *node) {
    node->index = -1;
}

static bool heapNode_valid(const HeapNode *node) {
    return node->index != -1;
}

static bool heap_empty(const Heap *q) {
    return q->size == 0;
}

static int64_t heap_minPrio(const Heap *q) {
    assert(q->size > 0);
    return q->heap[0]->prio;
}

static int64_t heap_minElem(const Heap *q) {
    assert(q->size > 0);
    return q->heap[0]->data;
}

static void heap_clear(Heap *q) {
    for (int64_t i = 0; i < q->size; i++) {
        q->heap[i]->index = -1;
    }
    q->size = 0;
}

static void heap_siftUp(Heap *q, int64_t index) {
    HeapNode *node = q->heap[index];
    int64_t prio = node->prio;
    while (index > 0) {
        int64_t nextIndex = (index - 1) / 2;
        HeapNode *nextNode = q->heap[nextIndex];
        if (nextNode->prio <= prio) {
            break;
        }
        q->heap[index] = nextNode;
        nextNode->index = index;
        index = nextIndex;
    }
    node->index = index;
    q->heap[index] = node;
}

static void heap_siftDown(Heap *q, int64_t index) {
    HeapNode *node = q->heap[index];
    int64_t prio = node->prio;
    int64_t elementNumber = q->size;
    while (index < elementNumber / 2) {
        int64_t nextIndex = 2 * index + 1;
        HeapNode *nextNode = q->heap[nextIndex];
        if (nextIndex + 1 < elementNumber) {
            HeapNode *otherNode = q->heap[nextIndex + 1];
            if (otherNode->prio <= nextNode->prio) {
                nextIndex++;
                nextNode = otherNode;
            }
        }
        if (nextNode->prio >= prio) {
            break;
        }
        q->heap[index] = nextNode;
        nextNode->index = index;
        index = nextIndex;
    }
    q->heap[index] = node;
    node->index = index;
}

static void heap_insert(Heap *q, HeapNode *node, int64_t prio, int64_t data) {
    assert(node->index == -1);
    if (q->size == q->capacity) {
        q->capacity = q->capacity == 0 ? 4 : 2 * q->capacity;
        q->heap = st_realloc(q->heap, q->capacity * sizeof(HeapNode *));
    }
    node->index = q->size;
    node->prio = prio;
    node->data = data;
    q->heap[q->size++] = node;
    heap_siftUp(q, node->index);
}

static void heap_setPrio(Heap *q, HeapNode *node, int64_t prio) {
    int64_t index = node->index;
    assert(index != -1);
    assert(q->heap[index] == node);
    int64_t previousPrio = node->prio;
    node->prio = prio;
    if (prio < previousPrio) {
        heap_siftUp(q, index);
    } else if (prio > previousPrio) {
        heap_siftDown(q, index);
    }
}

static void heap_remove(Heap *q, HeapNode *node) {
    int64_t index = node->index;
    assert(index != -1);
    assert(q->heap[index] == node);
    node->index = -1;
    HeapNode *moveNode = q->heap[--q->size];
    if (index < q->size) {
        q->heap[index] = moveNode;
        moveNode->index = index;
        if (moveNode->prio < node->prio) {
            heap_siftUp(q, index);
        } else if (moveNode->prio > node->prio) {
            heap_siftDown(q, index);
        }
    }
}

////////////////////////////////////
// Concatenable queue
////////////////////////////////////

/*
 * A priority queue that supports merging queues and undoing a merge, implemented as an AVL tree
 * with the minimum priority of each subtree kept at its root. Each queue has a name (here a
 * blossom), so that finding the queue holding an element finds the element's top-level blossom.
 * See Blelloch, Ferizovic and Sun, "Parallel ordered sets using join", arXiv:1602.02120, for the
 * join based tree operations.
 */

typedef struct _concatenableQueue ConcatenableQueue;

typedef struct _queueNode {
    int64_t prio;
    int64_t data;
    int64_t minPrio; //the minimum priority in the subtree of the node, and its data
    int64_t minData;
    ConcatenableQueue *owner; //set only at the root of a queue's tree
    struct _queueNode *parent;
    struct _queueNode *left;
    struct _queueNode *right;
    int64_t height; //0 if the node is not in a queue
} QueueNode;

struct _concatenableQueue {
    int64_t name;
    QueueNode *tree; //root node, or NULL if the queue is empty
    QueueNode *firstNode; //left-most node that belongs to this queue
    ConcatenableQueue *firstSubqueue; //first of the queues merged into this one
    ConcatenableQueue *nextSubqueue; //next queue merged into the same queue as this one
};

static void queueNode_init(QueueNode *node) {
    node->owner = NULL;
    node->parent = NULL;
    node->left = NULL;
    node->right = NULL;
    node->height = 0;
}

static void queue_init(ConcatenableQueue *q, int64_t name) {
    q->name = name;
    q->tree = NULL;
    q->firstNode = NULL;
    q->firstSubqueue = NULL;
    q->nextSubqueue = NULL;
}

static int64_t queueNode_find(const QueueNode *node) {
    while (node->parent != NULL) {
        node = node->parent;
    }
    assert(node->owner != NULL);
    return node->owner->name;
}

static int64_t queueNode_prio(const QueueNode *node) {
    assert(node->height != 0);
    return node->prio;
}

static void queueNode_setPrio(QueueNode *node, int64_t prio) {
    assert(node->height != 0);
    node->prio = prio;
    while (node != NULL) {
        int64_t minPrio = node->prio;
        int64_t minData = node->data;
        if (node->left != NULL && node->left->minPrio < minPrio) {
            minPrio = node->left->minPrio;
            minData = node->left->minData;
        }
        if (node->right != NULL && node->right->minPrio < minPrio) {
            minPrio = node->right->minPrio;
            minData = node->right->minData;
        }
        node->minPrio = minPrio;
        node->minData = minData;
        node = node->parent;
    }
}

static void queue_insert(ConcatenableQueue *q, QueueNode *node, int64_t prio, int64_t data) {
    assert(q->tree == NULL && q->firstNode == NULL);
    assert(node->height == 0 && node->parent == NULL && node->left == NULL && node->right == NULL);
    node->prio = prio;
    node->data = data;
    node->minPrio = prio;
    node->minData = data;
    node->owner = q;
    node->height = 1;
    q->tree = node;
    q->firstNode = node;
}

static int64_t queue_minPrio(const ConcatenableQueue *q) {
    assert(q->tree != NULL);
    return q->tree->minPrio;
}

static int64_t queue_minElem(const ConcatenableQueue *q) {
    assert(q->tree != NULL);
    return q->tree->minData;
}

static int64_t queueNode_height(const QueueNode *node) {
    return node != NULL ? node->height : 0;
}

static void queueNode_repair(QueueNode *node) {
    /*
     * Repairs the height and minimum priority of a node after its children changed.
     */
    int64_t leftHeight = queueNode_height(node->left), rightHeight = queueNode_height(node->right);
    node->height = 1 + (leftHeight > rightHeight ? leftHeight : rightHeight);
    int64_t minPrio = node->prio;
    int64_t minData = node->data;
    if (node->left != NULL && node->left->minPrio < minPrio) {
        minPrio = node->left->minPrio;
        minData = node->left->minData;
    }
    if (node->right != NULL && node->right->minPrio < minPrio) {
        minPrio = node->right->minPrio;
        minData = node->right->minData;
    }
    node->minPrio = minPrio;
    node->minData = minData;
}

static void queueNode_replaceChild(QueueNode *parent, QueueNode *oldChild, QueueNode *newChild) {
    if (parent != NULL) {
        if (parent->left == oldChild) {
            parent->left = newChild;
        } else {
            assert(parent->right == oldChild);
            parent->right = newChild;
        }
    }
}

static QueueNode *queueNode_rotateLeft(QueueNode *node) {
    /*
     *     N                C
     *    / \              / \
     *   A   C    --->    N   D
     *      / \          / \
     *     B   D        A   B
     */
    QueueNode *parent = node->parent;
    QueueNode *newTop = node->right;
    assert(newTop != NULL);
    QueueNode *b = newTop->left;
    node->right = b;
    if (b != NULL) {
        b->parent = node;
    }
    newTop->left = node;
    node->parent = newTop;
    newTop->parent = parent;
    queueNode_replaceChild(parent, node, newTop);
    queueNode_repair(node);
    queueNode_repair(newTop);
    return newTop;
}

static QueueNode *queueNode_rotateRight(QueueNode *node) {
    /*
     *       N                B
     *      / \              / \
     *     B   D    --->    A   N
     *    / \                  / \
     *   A   C                C   D
     */
    QueueNode *parent = node->parent;
    QueueNode *newTop = node->left;
    assert(newTop != NULL);
    QueueNode *c = newTop->right;
    node->left = c;
    if (c != NULL) {
        c->parent = node;
    }
    newTop->right = node;
    node->parent = newTop;
    newTop->parent = parent;
    queueNode_replaceChild(parent, node, newTop);
    queueNode_repair(node);
    queueNode_repair(newTop);
    return newTop;
}

static QueueNode *queueNode_joinRight(QueueNode *leftTree, QueueNode *node, QueueNode *rightTree) {
    /*
     * Joins two trees and a middle node where the left tree is the higher: descends the right
     * spine of the left tree to a node of compatible height, puts the middle node and the right
     * tree there, then rebalances on the way back up.
     */
    assert(leftTree != NULL);
    int64_t rightHeight = queueNode_height(rightTree);
    assert(leftTree->height > rightHeight + 1);
    QueueNode *cur = leftTree;
    while (cur->right != NULL && cur->right->height > rightHeight + 1) {
        cur = cur->right;
    }
    node->left = cur->right;
    node->right = rightTree;
    if (node->left != NULL) {
        node->left->parent = node;
    }
    if (rightTree != NULL) {
        rightTree->parent = node;
    }
    cur->right = node;
    node->parent = cur;
    if (cur->left == NULL || cur->left->height <= rightHeight) { //A double rotation is needed
        node = queueNode_rotateRight(node);
        cur = queueNode_rotateLeft(cur);
    } else {
        queueNode_repair(node);
        queueNode_repair(cur);
    }
    while (cur->parent != NULL) {
        cur = cur->parent;
        assert(cur->left != NULL && cur->right != NULL);
        if (cur->left->height + 1 < cur->right->height) {
            cur = queueNode_rotateLeft(cur);
        } else {
            queueNode_repair(cur);
        }
    }
    return cur;
}

static QueueNode *queueNode_joinLeft(QueueNode *leftTree, QueueNode *node, QueueNode *rightTree) {
    /*
     * The mirror image of queueNode_joinRight, for when the right tree is the higher.
     */
    assert(rightTree != NULL);
    int64_t leftHeight = queueNode_height(leftTree);
    assert(leftHeight + 1 < rightTree->height);
    QueueNode *cur = rightTree;
    while (cur->left != NULL && cur->left->height > leftHeight + 1) {
        cur = cur->left;
    }
    node->left = leftTree;
    node->right = cur->left;
    if (leftTree != NULL) {
        leftTree->parent = node;
    }
    if (node->right != NULL) {
        node->right->parent = node;
    }
    cur->left = node;
    node->parent = cur;
    if (cur->right == NULL || cur->right->height <= leftHeight) { //A double rotation is needed
        node = queueNode_rotateLeft(node);
        cur = queueNode_rotateRight(cur);
    } else {
        queueNode_repair(node);
        queueNode_repair(cur);
    }
    while (cur->parent != NULL) {
        cur = cur->parent;
        assert(cur->left != NULL && cur->right != NULL);
        if (cur->left->height > cur->right->height + 1) {
            cur = queueNode_rotateRight(cur);
        } else {
            queueNode_repair(cur);
        }
    }
    return cur;
}

static QueueNode *queueNode_join(QueueNode *leftTree, QueueNode *node, QueueNode *rightTree) {
    /*
     * Joins a left tree, a middle node and a right tree, in that order, and returns the root. The
     * trees may initially be children of the middle node; those links are broken as needed.
     */
    int64_t leftHeight = queueNode_height(leftTree), rightHeight = queueNode_height(rightTree);
    if (leftHeight > rightHeight + 1) {
        leftTree->parent = NULL;
        return queueNode_joinRight(leftTree, node, rightTree);
    }
    if (leftHeight + 1 < rightHeight) {
        rightTree->parent = NULL;
        return queueNode_joinLeft(leftTree, node, rightTree);
    }
    node->parent = NULL;
    node->left = leftTree;
    if (leftTree != NULL) {
        leftTree->parent = node;
    }
    node->right = rightTree;
    if (rightTree != NULL) {
        rightTree->parent = node;
    }
    queueNode_repair(node);
    return node;
}

static QueueNode *queueNode_splitLastNode(QueueNode *tree, QueueNode **lastNode) {
    /*
     * Removes the last node from the tree, returning the rest of the tree rebalanced.
     */
    assert(tree != NULL);
    QueueNode *last = tree;
    while (last->right != NULL) {
        last = last->right;
    }
    //Detach its left subtree and its parent
    QueueNode *leftTree = last->left;
    last->left = NULL;
    if (leftTree != NULL) {
        leftTree->parent = NULL;
    }
    QueueNode *parent = last->parent;
    last->parent = NULL;
    if (parent != NULL) {
        parent->right = NULL;
    }
    //Rebuild along the right spine of the original tree
    while (parent != NULL) {
        QueueNode *cur = parent;
        parent = cur->parent;
        cur->parent = NULL;
        if (parent != NULL) {
            parent->right = NULL;
        }
        leftTree = queueNode_join(cur->left, cur, leftTree);
    }
    *lastNode = last;
    return leftTree;
}

static QueueNode *queueNode_mergeTrees(QueueNode *leftTree, QueueNode *rightTree) {
    QueueNode *node;
    leftTree = queueNode_splitLastNode(leftTree, &node);
    return queueNode_join(leftTree, node, rightTree);
}

static void queueNode_splitTree(QueueNode *splitNode, QueueNode **leftTreeOut, QueueNode **rightTreeOut) {
    /*
     * Splits a tree into the nodes left of splitNode and the rest, splitNode included. Works up from
     * splitNode, joining each node and its remaining children to one side or the other; the joins
     * take O(log(n)) time in total.
     */
    assert(splitNode != NULL);
    QueueNode *leftTree = splitNode->left;
    splitNode->left = NULL;
    if (leftTree != NULL) {
        leftTree->parent = NULL;
    }
    QueueNode *rightTree = NULL;
    QueueNode *node = splitNode;
    bool leftBranch = 1; //The cut runs down the left branch of splitNode, as it goes to the right
    while (node != NULL) {
        QueueNode *parent = node->parent;
        node->parent = NULL;
        bool parentLeftBranch = 1;
        if (parent != NULL) {
            parentLeftBranch = parent->left == node;
            if (parentLeftBranch) {
                parent->left = NULL;
            } else {
                assert(parent->right == node);
                parent->right = NULL;
            }
        }
        if (leftBranch) { //The node goes to the right tree
            assert(node->left == NULL);
            rightTree = queueNode_join(rightTree, node, node->right);
        } else { //The node goes to the left tree
            assert(node->right == NULL);
            leftTree = queueNode_join(node->left, node, leftTree);
        }
        node = parent;
        leftBranch = parentLeftBranch;
    }
    *leftTreeOut = leftTree;
    *rightTreeOut = rightTree;
}

static void queue_merge(ConcatenableQueue *q, ConcatenableQueue **subqueues, int64_t subqueueNumber) {
    /*
     * Moves the elements of the given non-empty queues, in order, into the empty queue q. q keeps
     * a list of the sub-queues so that queue_split can undo the merge.
     */
    assert(q->tree == NULL && q->firstNode == NULL);
    assert(subqueueNumber > 0);
    ConcatenableQueue *sub = subqueues[0];
    assert(sub->tree != NULL && sub->tree->owner == sub);
    QueueNode *mergedTree = sub->tree;
    sub->tree = NULL;
    mergedTree->owner = NULL;
    assert(sub->firstNode != NULL);
    q->firstNode = sub->firstNode;
    assert(sub->nextSubqueue == NULL);
    q->firstSubqueue = sub;
    for (int64_t i = 1; i < subqueueNumber; i++) {
        ConcatenableQueue *previousSub = sub;
        sub = subqueues[i];
        assert(sub->tree != NULL && sub->tree->owner == sub);
        sub->tree->owner = NULL;
        mergedTree = queueNode_mergeTrees(mergedTree, sub->tree);
        sub->tree = NULL;
        assert(sub->nextSubqueue == NULL);
        previousSub->nextSubqueue = sub;
    }
    mergedTree->owner = q;
    q->tree = mergedTree;
}

static void queue_split(ConcatenableQueue *q) {
    /*
     * Undoes the merge that made q, putting the elements back in the queues they came from.
     */
    assert(q->tree != NULL && q->firstSubqueue != NULL);
    assert(q->tree->owner == q);
    q->tree->owner = NULL;
    QueueNode *rightTree = q->tree;
    ConcatenableQueue *sub = q->firstSubqueue;
    while (sub->nextSubqueue != NULL) {
        ConcatenableQueue *nextSub = sub->nextSubqueue;
        sub->nextSubqueue = NULL;
        assert(nextSub->firstNode != NULL);
        QueueNode *leftTree;
        queueNode_splitTree(nextSub->firstNode, &leftTree, &rightTree);
        sub->tree = leftTree;
        leftTree->owner = sub;
        sub = nextSub;
    }
    rightTree->owner = sub;
    sub->tree = rightTree;
    q->tree = NULL;
    q->firstNode = NULL;
    q->firstSubqueue = NULL;
}

////////////////////////////////////
// Blossoms and the matching state
////////////////////////////////////

typedef struct _blossom {
    int64_t parent; //the blossom this is a sub-blossom of, or NO_BLOSSOM if it is top-level
    int64_t baseVertex; //the vertex of the blossom not matched to another vertex in it
    int64_t label; //LABEL_S or LABEL_T for a top-level blossom in an alternating tree
    bool isNonTrivial;
    int64_t treeEdgeX, treeEdgeY; //the edge (x, y), y in this blossom, attaching it to its tree; x is NO_VERTEX at a root
    int64_t treeRoot; //base vertex of the blossom at the root of its alternating tree
    ConcatenableQueue vertexQueue; //the vertices in the blossom
    HeapNode delta2Node; //node in the delta2 queue
    int64_t vertexDualOffset; //pending lazy updates to the duals of the vertices in the blossom
    //The rest is for non-trivial blossoms only
    int64_t subblossomNumber; //0 if the slot is unused
    int64_t *subblossoms; //the sub-blossoms round the alternating cycle, starting with the base's
    int64_t *subblossomEdges; //(x, y) pairs: x in subblossoms[i], y in the next sub-blossom
    int64_t dualVar; //modified dual: z(b) - 2 * deltaSum for S, z(b) + 2 * deltaSum for T, else z(b)
    HeapNode delta4Node; //node in the delta4 queue, for top-level T-blossoms
} Blossom;

typedef struct _matchingContext {
    int64_t vertexNumber;
    int64_t edgeNumber;
    const int64_t *edgeX, *edgeY, *edgeWeight;
    int64_t *adjacencyStarts; //the edges incident to x are adjacentEdges[adjacencyStarts[x] .. adjacencyStarts[x+1]-1]
    int64_t *adjacentEdges;

    int64_t *vertexMate; //the vertex x is matched to, or NO_VERTEX

    //Blossoms 0 .. n-1 are the trivial blossoms of the vertices, n .. 2n-1 slots for non-trivial ones
    Blossom *blossoms;
    int64_t *unusedBlossoms;
    int64_t unusedBlossomNumber;

    /*
     * Modified (doubled) dual of each vertex, invariant under delta steps:
     * u(x) + deltaSum for an S-vertex, u(x) - deltaSum - B(x).vertexDualOffset for a T-vertex and
     * u(x) - B(x).vertexDualOffset otherwise, where B(x) is x's top-level blossom.
     */
    int64_t *vertexDual;
    int64_t initVertexDual; //the initial vertex dual, the maximum edge weight
    int64_t deltaSum; //running sum of the delta steps

    QueueNode *vertexQueueNodes; //for each vertex, its node in the queue of its top-level blossom
    Heap *vertexSEdgeQueues; //for each T- or unlabelled vertex, its edges to S-vertices by modified slack
    HeapNode *vertexSEdgeNodes; //for each edge, its node in a vertexSEdgeQueue
    Heap delta2Queue; //unlabelled top-level blossoms with an edge to an S-blossom, by least slack plus deltaSum
    Heap delta3Queue; //edges between S-vertices in different top-level blossoms, by modified slack
    HeapNode *delta3Nodes; //for each edge, its node in the delta3 queue
    Heap delta4Queue; //top-level non-trivial T-blossoms, by modified dual

    int64_t *scanQueue; //S-vertices whose edges have yet to be scanned
    int64_t scanQueueLength, scanQueueCapacity;

    bool *vertexMarker; //markers placed while tracing alternating paths
    int64_t *markedVertices;
    int64_t markedVertexNumber;

    //Scratch space
    int64_t *path; //a deque of (x, y) pairs, an alternating path, from pathStart to pathEnd (exclusive)
    int64_t pathStart, pathEnd;
    int64_t *blossomVertices; //the vertices of a blossom, as listed by getBlossomVertices
    int64_t *blossomStack;
    int64_t *augmentStack; //(outer, entry) pairs of blossoms still to augment
    int64_t augmentStackLength;
    int64_t *subblossomBuffer;
    ConcatenableQueue **subqueueBuffer;
} MatchingContext;

static int64_t otherVertex(MatchingContext *ctx, int64_t e, int64_t x) {
    return ctx->edgeX[e] != x ? ctx->edgeX[e] : ctx->edgeY[e];
}

static int64_t getBlossomVertices(MatchingContext *ctx, int64_t b) {
    /*
     * Lists the vertices of blossom b in ctx->blossomVertices and returns how many there are.
     * The order, from an explicit stack, is that of the C++ code.
     */
    Blossom *blossom = &ctx->blossoms[b];
    if (!blossom->isNonTrivial) {
        ctx->blossomVertices[0] = blossom->baseVertex;
        return 1;
    }
    int64_t vertexNumber = 0, stackLength = 0;
    ctx->blossomStack[stackLength++] = b;
    while (stackLength > 0) {
        Blossom *b2 = &ctx->blossoms[ctx->blossomStack[--stackLength]];
        for (int64_t i = 0; i < b2->subblossomNumber; i++) {
            Blossom *sub = &ctx->blossoms[b2->subblossoms[i]];
            if (sub->isNonTrivial) {
                ctx->blossomStack[stackLength++] = b2->subblossoms[i];
            } else {
                ctx->blossomVertices[vertexNumber++] = sub->baseVertex;
            }
        }
    }
    return vertexNumber;
}

static int64_t topLevelBlossom(MatchingContext *ctx, int64_t x) {
    return queueNode_find(&ctx->vertexQueueNodes[x]);
}

static void pushScanQueue(MatchingContext *ctx, int64_t x) {
    if (ctx->scanQueueLength == ctx->scanQueueCapacity) {
        ctx->scanQueueCapacity *= 2;
        ctx->scanQueue = st_realloc(ctx->scanQueue, ctx->scanQueueCapacity * sizeof(int64_t));
    }
    ctx->scanQueue[ctx->scanQueueLength++] = x;
}

static void blossom_initTrivial(Blossom *blossom, int64_t x) {
    blossom->parent = NO_BLOSSOM;
    blossom->baseVertex = x;
    blossom->label = LABEL_NONE;
    blossom->isNonTrivial = 0;
    blossom->treeEdgeX = NO_VERTEX;
    blossom->treeEdgeY = NO_VERTEX;
    blossom->treeRoot = NO_VERTEX;
    queue_init(&blossom->vertexQueue, x);
    heapNode_init(&blossom->delta2Node);
    blossom->vertexDualOffset = 0;
    blossom->subblossomNumber = 0;
    blossom->subblossoms = NULL;
    blossom->subblossomEdges = NULL;
    blossom->dualVar = 0;
    heapNode_init(&blossom->delta4Node);
}

static void matchingContext_construct(MatchingContext *ctx, int64_t vertexNumber, int64_t edgeNumber,
        const int64_t *edgeX, const int64_t *edgeY, const int64_t *edgeWeight) {
    int64_t n = vertexNumber, m = edgeNumber;
    ctx->vertexNumber = n;
    ctx->edgeNumber = m;
    ctx->edgeX = edgeX;
    ctx->edgeY = edgeY;
    ctx->edgeWeight = edgeWeight;

    //Adjacency lists, each in the order of the edges
    ctx->adjacencyStarts = st_calloc(n + 1, sizeof(int64_t));
    for (int64_t e = 0; e < m; e++) {
        ctx->adjacencyStarts[edgeX[e] + 1]++;
        ctx->adjacencyStarts[edgeY[e] + 1]++;
    }
    for (int64_t x = 0; x < n; x++) {
        ctx->adjacencyStarts[x + 1] += ctx->adjacencyStarts[x];
    }
    ctx->adjacentEdges = st_malloc(2 * m * sizeof(int64_t));
    int64_t *fill = st_malloc(n * sizeof(int64_t));
    memcpy(fill, ctx->adjacencyStarts, n * sizeof(int64_t));
    for (int64_t e = 0; e < m; e++) {
        ctx->adjacentEdges[fill[edgeX[e]]++] = e;
        ctx->adjacentEdges[fill[edgeY[e]]++] = e;
    }
    free(fill);

    ctx->vertexMate = st_malloc(n * sizeof(int64_t));
    ctx->blossoms = st_malloc(2 * n * sizeof(Blossom));
    ctx->vertexQueueNodes = st_malloc(n * sizeof(QueueNode));
    for (int64_t x = 0; x < n; x++) {
        ctx->vertexMate[x] = NO_VERTEX;
        blossom_initTrivial(&ctx->blossoms[x], x);
        queueNode_init(&ctx->vertexQueueNodes[x]);
        queue_insert(&ctx->blossoms[x].vertexQueue, &ctx->vertexQueueNodes[x], INVALID_SLACK, x);
    }
    ctx->unusedBlossoms = st_malloc(n * sizeof(int64_t));
    ctx->unusedBlossomNumber = 0;
    for (int64_t b = 2 * n - 1; b >= n; b--) { //Slots are handed out lowest first
        blossom_initTrivial(&ctx->blossoms[b], NO_VERTEX);
        ctx->blossoms[b].isNonTrivial = 1;
        ctx->unusedBlossoms[ctx->unusedBlossomNumber++] = b;
    }

    //Vertex duals start at half the maximum edge weight
    int64_t maxWeight = 0;
    for (int64_t e = 0; e < m; e++) {
        if (edgeWeight[e] > maxWeight) {
            maxWeight = edgeWeight[e];
        }
    }
    ctx->initVertexDual = maxWeight;
    ctx->vertexDual = st_malloc(n * sizeof(int64_t));
    for (int64_t x = 0; x < n; x++) {
        ctx->vertexDual[x] = maxWeight;
    }
    ctx->deltaSum = 0;

    ctx->vertexSEdgeQueues = st_malloc(n * sizeof(Heap));
    for (int64_t x = 0; x < n; x++) {
        heap_init(&ctx->vertexSEdgeQueues[x]);
    }
    ctx->vertexSEdgeNodes = st_malloc(m * sizeof(HeapNode));
    ctx->delta3Nodes = st_malloc(m * sizeof(HeapNode));
    for (int64_t e = 0; e < m; e++) {
        heapNode_init(&ctx->vertexSEdgeNodes[e]);
        heapNode_init(&ctx->delta3Nodes[e]);
    }
    heap_init(&ctx->delta2Queue);
    heap_init(&ctx->delta3Queue);
    heap_init(&ctx->delta4Queue);

    ctx->scanQueueCapacity = n;
    ctx->scanQueueLength = 0;
    ctx->scanQueue = st_malloc(ctx->scanQueueCapacity * sizeof(int64_t));
    ctx->vertexMarker = st_calloc(n, sizeof(bool));
    ctx->markedVertices = st_malloc(n * sizeof(int64_t));
    ctx->markedVertexNumber = 0;

    ctx->path = st_malloc(2 * (4 * n + 4) * sizeof(int64_t));
    ctx->blossomVertices = st_malloc(n * sizeof(int64_t));
    ctx->blossomStack = st_malloc(n * sizeof(int64_t));
    ctx->augmentStack = st_malloc(2 * (2 * n + 2) * sizeof(int64_t));
    ctx->augmentStackLength = 0;
    ctx->subblossomBuffer = st_malloc(2 * (n + 1) * sizeof(int64_t));
    ctx->subqueueBuffer = st_malloc((n + 1) * sizeof(ConcatenableQueue *));
}

static void matchingContext_destruct(MatchingContext *ctx) {
    for (int64_t b = 0; b < 2 * ctx->vertexNumber; b++) {
        free(ctx->blossoms[b].subblossoms);
        free(ctx->blossoms[b].subblossomEdges);
    }
    for (int64_t x = 0; x < ctx->vertexNumber; x++) {
        heap_destruct(&ctx->vertexSEdgeQueues[x]);
    }
    free(ctx->adjacencyStarts);
    free(ctx->adjacentEdges);
    free(ctx->vertexMate);
    free(ctx->blossoms);
    free(ctx->unusedBlossoms);
    free(ctx->vertexDual);
    free(ctx->vertexQueueNodes);
    free(ctx->vertexSEdgeQueues);
    free(ctx->vertexSEdgeNodes);
    heap_destruct(&ctx->delta2Queue);
    heap_destruct(&ctx->delta3Queue);
    free(ctx->delta3Nodes);
    heap_destruct(&ctx->delta4Queue);
    free(ctx->scanQueue);
    free(ctx->vertexMarker);
    free(ctx->markedVertices);
    free(ctx->path);
    free(ctx->blossomVertices);
    free(ctx->blossomStack);
    free(ctx->augmentStack);
    free(ctx->subblossomBuffer);
    free(ctx->subqueueBuffer);
}

////////////////////////////////////
// Least slack edge tracking
////////////////////////////////////

static int64_t edgePseudoSlack(MatchingContext *ctx, int64_t e) {
    /*
     * The modified slack of edge e, which is invariant under delta steps.
     */
    return ctx->vertexDual[ctx->edgeX[e]] + ctx->vertexDual[ctx->edgeY[e]] - 2 * ctx->edgeWeight[e];
}

static void delta2AddEdge(MatchingContext *ctx, int64_t e, int64_t y, Blossom *by) {
    /*
     * Tracks edge e, from an S-vertex to the T- or unlabelled vertex y, for delta2.
     */
    int64_t prio = edgePseudoSlack(ctx, e);
    Heap *queue = &ctx->vertexSEdgeQueues[y];
    assert(!heapNode_valid(&ctx->vertexSEdgeNodes[e]));
    bool improved = heap_empty(queue) || prio < heap_minPrio(queue);
    heap_insert(queue, &ctx->vertexSEdgeNodes[e], prio, e);
    if (improved) {
        //Update y's priority in its blossom, and the blossom's in the delta2 queue if it is unlabelled
        queueNode_setPrio(&ctx->vertexQueueNodes[y], prio);
        if (by->label == LABEL_NONE) {
            prio += by->vertexDualOffset;
            if (heapNode_valid(&by->delta2Node)) {
                if (prio < by->delta2Node.prio) {
                    heap_setPrio(&ctx->delta2Queue, &by->delta2Node, prio);
                }
            } else {
                heap_insert(&ctx->delta2Queue, &by->delta2Node, prio, by - ctx->blossoms);
            }
        }
    }
}

static void delta2RemoveEdge(MatchingContext *ctx, int64_t e, int64_t y, Blossom *by) {
    /*
     * Stops tracking edge e for delta2, as its S-vertex has become unlabelled; y is the T- or
     * unlabelled vertex at its other end.
     */
    if (heapNode_valid(&ctx->vertexSEdgeNodes[e])) {
        Heap *queue = &ctx->vertexSEdgeQueues[y];
        heap_remove(queue, &ctx->vertexSEdgeNodes[e]);
        int64_t prio = heap_empty(queue) ? INVALID_SLACK : heap_minPrio(queue);
        if (prio > queueNode_prio(&ctx->vertexQueueNodes[y])) {
            queueNode_setPrio(&ctx->vertexQueueNodes[y], prio);
            if (by->label == LABEL_NONE) {
                //Update or remove the blossom in the delta2 queue
                assert(heapNode_valid(&by->delta2Node));
                prio = queue_minPrio(&by->vertexQueue);
                if (prio < INVALID_SLACK) {
                    prio += by->vertexDualOffset;
                    if (prio > by->delta2Node.prio) {
                        heap_setPrio(&ctx->delta2Queue, &by->delta2Node, prio);
                    }
                } else {
                    heap_remove(&ctx->delta2Queue, &by->delta2Node);
                }
            }
        }
    }
}

static void delta2EnableBlossom(MatchingContext *ctx, Blossom *blossom) {
    /*
     * Puts a blossom that has become an unlabelled top-level blossom in the delta2 queue, if it
     * has an edge to an S-vertex.
     */
    assert(!heapNode_valid(&blossom->delta2Node));
    int64_t prio = queue_minPrio(&blossom->vertexQueue);
    if (prio < INVALID_SLACK) {
        prio += blossom->vertexDualOffset;
        heap_insert(&ctx->delta2Queue, &blossom->delta2Node, prio, blossom - ctx->blossoms);
    }
}

static void delta2DisableBlossom(MatchingContext *ctx, Blossom *blossom) {
    if (heapNode_valid(&blossom->delta2Node)) {
        heap_remove(&ctx->delta2Queue, &blossom->delta2Node);
    }
}

static void delta2ClearVertex(MatchingContext *ctx, int64_t x) {
    /*
     * Forgets the edges tracked for x, which has become an S-vertex.
     */
    heap_clear(&ctx->vertexSEdgeQueues[x]);
    queueNode_setPrio(&ctx->vertexQueueNodes[x], INVALID_SLACK);
}

static int64_t delta2GetMinEdge(MatchingContext *ctx, int64_t *slack) {
    /*
     * Returns the least slack edge from an S-vertex to an unlabelled vertex, and its slack, or
     * NO_EDGE if there is none.
     */
    if (heap_empty(&ctx->delta2Queue)) {
        return NO_EDGE;
    }
    Blossom *blossom = &ctx->blossoms[heap_minElem(&ctx->delta2Queue)];
    assert(blossom->parent == NO_BLOSSOM);
    assert(blossom->label == LABEL_NONE);
    *slack = heap_minPrio(&ctx->delta2Queue) - ctx->deltaSum;
    int64_t x = queue_minElem(&blossom->vertexQueue);
    return heap_minElem(&ctx->vertexSEdgeQueues[x]);
}

static void delta3AddEdge(MatchingContext *ctx, int64_t e) {
    /*
     * Tracks edge e, between S-vertices in different top-level blossoms, for delta3. It may
     * already have been added from its other end.
     */
    if (!heapNode_valid(&ctx->delta3Nodes[e])) {
        heap_insert(&ctx->delta3Queue, &ctx->delta3Nodes[e], edgePseudoSlack(ctx, e), e);
    }
}

static void delta3RemoveEdge(MatchingContext *ctx, int64_t e) {
    if (heapNode_valid(&ctx->delta3Nodes[e])) {
        heap_remove(&ctx->delta3Queue, &ctx->delta3Nodes[e]);
    }
}

static int64_t delta3GetMinEdge(MatchingContext *ctx, int64_t *slack) {
    /*
     * Returns the least slack edge between S-vertices in different top-level blossoms, and its
     * slack, or NO_EDGE if there is none. Edges that merging blossoms has made internal to a
     * blossom are dropped from the queue as they are met.
     */
    while (!heap_empty(&ctx->delta3Queue)) {
        int64_t e = heap_minElem(&ctx->delta3Queue);
        int64_t bx = topLevelBlossom(ctx, ctx->edgeX[e]);
        int64_t by = topLevelBlossom(ctx, ctx->edgeY[e]);
        assert(ctx->blossoms[bx].label == LABEL_S && ctx->blossoms[by].label == LABEL_S);
        if (bx != by) {
            *slack = ctx->delta3Nodes[e].prio - 2 * ctx->deltaSum;
            return e;
        }
        heap_remove(&ctx->delta3Queue, &ctx->delta3Nodes[e]);
    }
    return NO_EDGE;
}

////////////////////////////////////
// Blossom labels
////////////////////////////////////

static void assignBlossomLabelS(MatchingContext *ctx, Blossom *blossom) {
    /*
     * Makes an unlabelled top-level blossom an S-blossom and queues its vertices for scanning.
     */
    assert(blossom->parent == NO_BLOSSOM);
    assert(blossom->label == LABEL_NONE);
    blossom->label = LABEL_S;
    delta2DisableBlossom(ctx, blossom);
    //Unlabelled and S-blossoms modify their duals differently, so adjust to keep the true duals
    if (blossom->isNonTrivial) {
        blossom->dualVar -= 2 * ctx->deltaSum;
    }
    int64_t dualFixup = ctx->deltaSum + blossom->vertexDualOffset;
    blossom->vertexDualOffset = 0;
    int64_t vertexNumber = getBlossomVertices(ctx, blossom - ctx->blossoms);
    for (int64_t i = 0; i < vertexNumber; i++) {
        int64_t x = ctx->blossomVertices[i];
        ctx->vertexDual[x] += dualFixup;
        delta2ClearVertex(ctx, x);
        pushScanQueue(ctx, x);
    }
}

static void assignBlossomLabelT(MatchingContext *ctx, Blossom *blossom) {
    /*
     * Makes an unlabelled top-level blossom a T-blossom.
     */
    assert(blossom->parent == NO_BLOSSOM);
    assert(blossom->label == LABEL_NONE);
    blossom->label = LABEL_T;
    delta2DisableBlossom(ctx, blossom);
    if (blossom->isNonTrivial) {
        blossom->dualVar += 2 * ctx->deltaSum;
        assert(!heapNode_valid(&blossom->delta4Node));
        heap_insert(&ctx->delta4Queue, &blossom->delta4Node, blossom->dualVar, blossom - ctx->blossoms);
    }
    blossom->vertexDualOffset -= ctx->deltaSum;
}

static void removeBlossomLabelS(MatchingContext *ctx, Blossom *blossom) {
    /*
     * Makes a top-level S-blossom unlabelled, moving the edges of its vertices from delta3 to
     * delta2 tracking as appropriate.
     */
    assert(blossom->parent == NO_BLOSSOM);
    assert(blossom->label == LABEL_S);
    blossom->label = LABEL_NONE;
    if (blossom->isNonTrivial) {
        blossom->dualVar += 2 * ctx->deltaSum;
    }
    assert(blossom->vertexDualOffset == 0);
    int64_t dualFixup = -ctx->deltaSum;
    int64_t vertexNumber = getBlossomVertices(ctx, blossom - ctx->blossoms);
    for (int64_t i = 0; i < vertexNumber; i++) {
        int64_t x = ctx->blossomVertices[i];
        ctx->vertexDual[x] += dualFixup;
        for (int64_t j = ctx->adjacencyStarts[x]; j < ctx->adjacencyStarts[x + 1]; j++) {
            int64_t e = ctx->adjacentEdges[j];
            int64_t y = otherVertex(ctx, e, x);
            delta3RemoveEdge(ctx, e);
            Blossom *by = &ctx->blossoms[topLevelBlossom(ctx, y)];
            if (by->label == LABEL_S) { //e now joins the unlabelled x to the S-vertex y
                delta2AddEdge(ctx, e, x, blossom);
            } else { //e was tracked for delta2 from x, which is no longer an S-vertex
                delta2RemoveEdge(ctx, e, y, by);
            }
        }
    }
}

static void removeBlossomLabelT(MatchingContext *ctx, Blossom *blossom) {
    /*
     * Makes a top-level T-blossom unlabelled.
     */
    assert(blossom->parent == NO_BLOSSOM);
    assert(blossom->label == LABEL_T);
    blossom->label = LABEL_NONE;
    if (blossom->isNonTrivial) {
        assert(heapNode_valid(&blossom->delta4Node));
        heap_remove(&ctx->delta4Queue, &blossom->delta4Node);
        blossom->dualVar -= 2 * ctx->deltaSum;
    }
    blossom->vertexDualOffset += ctx->deltaSum;
    delta2EnableBlossom(ctx, blossom);
}

static void changeSBlossomToSubblossom(MatchingContext *ctx, Blossom *blossom) {
    assert(blossom->parent == NO_BLOSSOM);
    assert(blossom->label == LABEL_S);
    blossom->label = LABEL_NONE;
    if (blossom->isNonTrivial) {
        blossom->dualVar += 2 * ctx->deltaSum;
    }
}

static void resetBlossomLabel(MatchingContext *ctx, Blossom *blossom) {
    if (blossom->parent == NO_BLOSSOM) {
        if (blossom->label == LABEL_S) {
            removeBlossomLabelS(ctx, blossom);
        } else if (blossom->label == LABEL_T) {
            removeBlossomLabelT(ctx, blossom);
        }
    }
}

static void removeAlternatingTrees(MatchingContext *ctx, int64_t root1, int64_t root2) {
    /*
     * Unlabels every blossom in the alternating trees with the given roots.
     */
    for (int64_t b = 0; b < 2 * ctx->vertexNumber; b++) {
        Blossom *blossom = &ctx->blossoms[b];
        if ((!blossom->isNonTrivial || blossom->subblossomNumber > 0) && blossom->parent == NO_BLOSSOM
                && blossom->label != LABEL_NONE && (blossom->treeRoot == root1 || blossom->treeRoot == root2)) {
            resetBlossomLabel(ctx, blossom);
        }
    }
}

////////////////////////////////////
// Creating and expanding blossoms
////////////////////////////////////

static void pathPushFront(MatchingContext *ctx, int64_t x, int64_t y) {
    ctx->pathStart--;
    ctx->path[2 * ctx->pathStart] = x;
    ctx->path[2 * ctx->pathStart + 1] = y;
}

static void pathPushBack(MatchingContext *ctx, int64_t x, int64_t y) {
    ctx->path[2 * ctx->pathEnd] = x;
    ctx->path[2 * ctx->pathEnd + 1] = y;
    ctx->pathEnd++;
}

static void traceAlternatingPaths(MatchingContext *ctx, int64_t x, int64_t y) {
    /*
     * Traces back through the alternating trees from the S-vertices x and y, which a tight edge
     * joins, and leaves the path in ctx->path. If they are in the same tree the path is a cycle
     * through a new blossom, from their first common ancestor back to it; if not, it is an
     * augmenting path between the two roots. The two traces alternate, so finding a blossom takes
     * time proportional to its size.
     */
    ctx->pathStart = ctx->pathEnd = 2 * ctx->vertexNumber + 2;
    pathPushBack(ctx, x, y);
    int64_t firstCommon = NO_BLOSSOM;
    while (x != NO_VERTEX || y != NO_VERTEX) {
        if (x != NO_VERTEX) {
            int64_t bx = topLevelBlossom(ctx, x);
            Blossom *blossom = &ctx->blossoms[bx];
            if (ctx->vertexMarker[blossom->baseVertex]) {
                firstCommon = bx;
                break;
            }
            ctx->vertexMarker[blossom->baseVertex] = 1;
            ctx->markedVertices[ctx->markedVertexNumber++] = blossom->baseVertex;
            x = blossom->treeEdgeX;
            if (x != NO_VERTEX) {
                pathPushFront(ctx, blossom->treeEdgeX, blossom->treeEdgeY);
            }
        }
        if (y != NO_VERTEX) {
            int64_t by = topLevelBlossom(ctx, y);
            Blossom *blossom = &ctx->blossoms[by];
            if (ctx->vertexMarker[blossom->baseVertex]) {
                firstCommon = by;
                break;
            }
            ctx->vertexMarker[blossom->baseVertex] = 1;
            ctx->markedVertices[ctx->markedVertexNumber++] = blossom->baseVertex;
            y = blossom->treeEdgeX;
            if (y != NO_VERTEX) {
                pathPushBack(ctx, blossom->treeEdgeY, y);
            }
        }
    }
    for (int64_t i = 0; i < ctx->markedVertexNumber; i++) {
        ctx->vertexMarker[ctx->markedVertices[i]] = 0;
    }
    ctx->markedVertexNumber = 0;
    //Trim the paths so that both end at the common ancestor
    if (firstCommon != NO_BLOSSOM) {
        assert(ctx->blossoms[firstCommon].label == LABEL_S);
        while (topLevelBlossom(ctx, ctx->path[2 * ctx->pathStart]) != firstCommon) {
            ctx->pathStart++;
        }
        while (topLevelBlossom(ctx, ctx->path[2 * (ctx->pathEnd - 1) + 1]) != firstCommon) {
            ctx->pathEnd--;
        }
    }
    assert((ctx->pathEnd - ctx->pathStart) % 2 == 1); //A path between S-blossoms has odd length
}

static void makeBlossom(MatchingContext *ctx) {
    /*
     * Makes a new S-blossom from the alternating cycle in ctx->path, relabelling its T-sub-blossoms
     * as S and queueing their vertices for scanning.
     */
    int64_t subblossomNumber = ctx->pathEnd - ctx->pathStart;
    assert(subblossomNumber % 2 == 1 && subblossomNumber >= 3);
    int64_t *subblossoms = st_malloc(subblossomNumber * sizeof(int64_t));
    int64_t *subblossomEdges = st_malloc(2 * subblossomNumber * sizeof(int64_t));
    for (int64_t i = 0; i < subblossomNumber; i++) {
        subblossomEdges[2 * i] = ctx->path[2 * (ctx->pathStart + i)];
        subblossomEdges[2 * i + 1] = ctx->path[2 * (ctx->pathStart + i) + 1];
        subblossoms[i] = topLevelBlossom(ctx, subblossomEdges[2 * i]);
    }
#ifndef NDEBUG
    for (int64_t i = 0; i < subblossomNumber; i++) { //Check the path is a cycle
        assert(topLevelBlossom(ctx, subblossomEdges[2 * i + 1]) == subblossoms[(i + 1) % subblossomNumber]);
    }
#endif
    assert(ctx->blossoms[subblossoms[0]].label == LABEL_S);
    for (int64_t i = 0; i < subblossomNumber; i++) {
        Blossom *sub = &ctx->blossoms[subblossoms[i]];
        if (sub->label == LABEL_T) { //The vertices of former T-blossoms become S-vertices
            removeBlossomLabelT(ctx, sub);
            assignBlossomLabelS(ctx, sub);
        }
        changeSBlossomToSubblossom(ctx, sub);
    }

    assert(ctx->unusedBlossomNumber > 0);
    int64_t b = ctx->unusedBlossoms[--ctx->unusedBlossomNumber];
    Blossom *blossom = &ctx->blossoms[b];
    assert(blossom->isNonTrivial && blossom->subblossomNumber == 0);
    blossom->parent = NO_BLOSSOM;
    blossom->baseVertex = ctx->blossoms[subblossoms[0]].baseVertex;
    blossom->label = LABEL_NONE;
    blossom->treeRoot = NO_VERTEX;
    queue_init(&blossom->vertexQueue, b);
    heapNode_init(&blossom->delta2Node);
    blossom->vertexDualOffset = 0;
    blossom->subblossomNumber = subblossomNumber;
    blossom->subblossoms = subblossoms;
    blossom->subblossomEdges = subblossomEdges;
    blossom->dualVar = -2 * ctx->deltaSum; //So that the true dual starts at zero
    heapNode_init(&blossom->delta4Node);

    for (int64_t i = 0; i < subblossomNumber; i++) {
        ctx->blossoms[subblossoms[i]].parent = b;
        ctx->subqueueBuffer[i] = &ctx->blossoms[subblossoms[i]].vertexQueue;
    }
    queue_merge(&blossom->vertexQueue, ctx->subqueueBuffer, subblossomNumber);

    blossom->label = LABEL_S;
    blossom->treeEdgeX = ctx->blossoms[subblossoms[0]].treeEdgeX;
    blossom->treeEdgeY = ctx->blossoms[subblossoms[0]].treeEdgeY;
    blossom->treeRoot = ctx->blossoms[subblossoms[0]].treeRoot;
}

static int64_t findSubblossom(Blossom *blossom, int64_t sub) {
    int64_t i = 0;
    while (blossom->subblossoms[i] != sub) {
        i++;
        assert(i < blossom->subblossomNumber);
    }
    return i;
}

static void extendTreeSToT(MatchingContext *ctx, int64_t x, int64_t y);

static void expandTBlossom(MatchingContext *ctx, Blossom *blossom) {
    /*
     * Expands a top-level T-blossom whose dual has reached zero, relinks into its alternating tree
     * the sub-blossoms on the even length path from where it was attached to its base, and frees
     * it.
     */
    assert(blossom->parent == NO_BLOSSOM);
    assert(blossom->label == LABEL_T);
    removeBlossomLabelT(ctx, blossom);
    delta2DisableBlossom(ctx, blossom);
    queue_split(&blossom->vertexQueue); //This restores the queues of the sub-blossoms
    int64_t vertexDualOffset = blossom->vertexDualOffset; //To be pushed down to the sub-blossoms
    blossom->vertexDualOffset = 0;
    for (int64_t i = 0; i < blossom->subblossomNumber; i++) {
        Blossom *sub = &ctx->blossoms[blossom->subblossoms[i]];
        assert(sub->parent == blossom - ctx->blossoms);
        assert(sub->label == LABEL_NONE);
        sub->parent = NO_BLOSSOM;
        assert(sub->vertexDualOffset == 0);
        sub->vertexDualOffset = vertexDualOffset;
        delta2EnableBlossom(ctx, sub);
    }

    //Find the sub-blossom by which the blossom was attached to its tree, and walk from it to the
    //base in an even number of steps, labelling the sub-blossoms on the way
    int64_t entry = topLevelBlossom(ctx, blossom->treeEdgeY);
    int64_t entryPosition = findSubblossom(blossom, entry);
    int64_t x = blossom->treeEdgeX, y = blossom->treeEdgeY;
    int64_t k = blossom->subblossomNumber;
    int64_t position = entryPosition;
    while (position != 0 && position != k) {
        extendTreeSToT(ctx, x, y); //Labels the current sub-blossom T and the next one S
        if (entryPosition % 2 == 0) { //Walk two steps backward to the base
            position -= 2;
            assert(position >= 0);
            y = blossom->subblossomEdges[2 * position];
            x = blossom->subblossomEdges[2 * position + 1];
        } else { //Walk two steps forward to the base
            position++;
            assert(position < k);
            x = blossom->subblossomEdges[2 * position];
            y = blossom->subblossomEdges[2 * position + 1];
            position++;
        }
    }
    //Finally label the base sub-blossom T
    Blossom *base = &ctx->blossoms[blossom->subblossoms[0]];
    assignBlossomLabelT(ctx, base);
    base->treeEdgeX = x;
    base->treeEdgeY = y;
    base->treeRoot = blossom->treeRoot;

    //Free the blossom
    free(blossom->subblossoms);
    free(blossom->subblossomEdges);
    blossom->subblossoms = NULL;
    blossom->subblossomEdges = NULL;
    blossom->subblossomNumber = 0;
    blossom->parent = NO_BLOSSOM;
    blossom->label = LABEL_NONE;
    ctx->unusedBlossoms[ctx->unusedBlossomNumber++] = blossom - ctx->blossoms;
}

////////////////////////////////////
// Augmenting
////////////////////////////////////

static void pushAugmentStack(MatchingContext *ctx, int64_t outer, int64_t entry) {
    ctx->augmentStack[2 * ctx->augmentStackLength] = outer;
    ctx->augmentStack[2 * ctx->augmentStackLength + 1] = entry;
    ctx->augmentStackLength++;
}

static void augmentBlossomRec(MatchingContext *ctx, int64_t b, int64_t entry) {
    /*
     * Augments along the alternating path through blossom b from its sub-blossom entry to its
     * base, then rotates the blossom so that entry holds the new base. The non-trivial
     * sub-blossoms the path touches are queued on the augment stack.
     */
    Blossom *blossom = &ctx->blossoms[b];
    int64_t k = blossom->subblossomNumber;
    int64_t entryPosition = findSubblossom(blossom, entry);
    int64_t position = entryPosition;
    while (position != 0 && position != k) {
        int64_t x, y, bx, by;
        if (entryPosition % 2 == 0) { //Walk backward to the base
            by = blossom->subblossoms[position - 1];
            position -= 2;
            assert(position >= 0);
            bx = blossom->subblossoms[position];
            x = blossom->subblossomEdges[2 * position];
            y = blossom->subblossomEdges[2 * position + 1];
        } else { //Walk forward to the base, which may wrap past the end of the list
            position++;
            assert(position < k);
            x = blossom->subblossomEdges[2 * position];
            y = blossom->subblossomEdges[2 * position + 1];
            bx = blossom->subblossoms[position];
            position++;
            by = position == k ? blossom->subblossoms[0] : blossom->subblossoms[position];
        }
        //Pull this edge into the matching
        ctx->vertexMate[x] = y;
        ctx->vertexMate[y] = x;
        //Augment through any non-trivial sub-blossoms touching the edge
        if (ctx->blossoms[bx].isNonTrivial) {
            pushAugmentStack(ctx, bx, x);
        }
        if (ctx->blossoms[by].isNonTrivial) {
            pushAugmentStack(ctx, by, y);
        }
    }
    //Rotate the sub-blossoms to put entry first
    if (entryPosition != 0) {
        int64_t *buffer = ctx->subblossomBuffer;
        memcpy(buffer, blossom->subblossoms, entryPosition * sizeof(int64_t));
        memmove(blossom->subblossoms, blossom->subblossoms + entryPosition, (k - entryPosition) * sizeof(int64_t));
        memcpy(blossom->subblossoms + k - entryPosition, buffer, entryPosition * sizeof(int64_t));
        int64_t *edges = blossom->subblossomEdges;
        memcpy(buffer, edges, 2 * entryPosition * sizeof(int64_t));
        memmove(edges, edges + 2 * entryPosition, 2 * (k - entryPosition) * sizeof(int64_t));
        memcpy(edges + 2 * (k - entryPosition), buffer, 2 * entryPosition * sizeof(int64_t));
    }
    //The entry sub-blossom has already been augmented, so its base is the new base
    blossom->baseVertex = ctx->blossoms[entry].baseVertex;
}

static void augmentBlossom(MatchingContext *ctx, int64_t b, int64_t entry) {
    /*
     * Augments along the alternating path through blossom b from the sub-blossom entry (at any
     * depth) to its base, handling nested blossoms with an explicit stack.
     */
    ctx->augmentStackLength = 0;
    pushAugmentStack(ctx, b, entry);
    while (ctx->augmentStackLength > 0) {
        int64_t top = ctx->augmentStackLength - 1;
        int64_t outer = ctx->augmentStack[2 * top];
        int64_t innerEntry = ctx->augmentStack[2 * top + 1];
        int64_t inner = ctx->blossoms[innerEntry].parent;
        assert(inner != NO_BLOSSOM);
        if (inner != outer) { //After inner, continue with its parent
            ctx->augmentStack[2 * top + 1] = inner;
        } else { //After inner, outer is done
            ctx->augmentStackLength--;
        }
        augmentBlossomRec(ctx, inner, innerEntry);
    }
}

static void augmentMatching(MatchingContext *ctx) {
    /*
     * Augments the matching along the augmenting path in ctx->path, whose even numbered edges are
     * the unmatched ones.
     */
    assert((ctx->pathEnd - ctx->pathStart) % 2 == 1);
    for (int64_t i = ctx->pathStart; i < ctx->pathEnd; i += 2) {
        int64_t x = ctx->path[2 * i], y = ctx->path[2 * i + 1];
        int64_t bx = topLevelBlossom(ctx, x);
        if (ctx->blossoms[bx].isNonTrivial) {
            augmentBlossom(ctx, bx, x);
        }
        int64_t by = topLevelBlossom(ctx, y);
        if (ctx->blossoms[by].isNonTrivial) {
            augmentBlossom(ctx, by, y);
        }
        ctx->vertexMate[x] = y;
        ctx->vertexMate[y] = x;
    }
}

////////////////////////////////////
// Alternating trees
////////////////////////////////////

static void extendTreeSToT(MatchingContext *ctx, int64_t x, int64_t y) {
    /*
     * Labels T the unlabelled blossom holding y, which a tight edge joins to the S-vertex x, and
     * labels S the blossom its base is matched into.
     */
    Blossom *bx = &ctx->blossoms[topLevelBlossom(ctx, x)];
    Blossom *by = &ctx->blossoms[topLevelBlossom(ctx, y)];
    assert(bx->label == LABEL_S);
    assignBlossomLabelT(ctx, by);
    by->treeEdgeX = x;
    by->treeEdgeY = y;
    by->treeRoot = bx->treeRoot;
    int64_t y2 = by->baseVertex;
    int64_t z = ctx->vertexMate[y2];
    assert(z != NO_VERTEX);
    Blossom *bz = &ctx->blossoms[topLevelBlossom(ctx, z)];
    assignBlossomLabelS(ctx, bz);
    bz->treeEdgeX = y2;
    bz->treeEdgeY = z;
    bz->treeRoot = by->treeRoot;
}

static bool addSToSEdge(MatchingContext *ctx, int64_t x, int64_t y) {
    /*
     * Uses a tight edge between S-vertices in different top-level blossoms: within one tree it
     * makes a blossom, between two trees it augments the matching, and returns true if it did.
     */
    Blossom *bx = &ctx->blossoms[topLevelBlossom(ctx, x)];
    Blossom *by = &ctx->blossoms[topLevelBlossom(ctx, y)];
    assert(bx->label == LABEL_S && by->label == LABEL_S);
    assert(bx != by);
    traceAlternatingPaths(ctx, x, y);
    if (bx->treeRoot == by->treeRoot) {
        makeBlossom(ctx);
        return 0;
    }
    removeAlternatingTrees(ctx, bx->treeRoot, by->treeRoot);
    augmentMatching(ctx);
    return 1;
}

static void scanNewSVertices(MatchingContext *ctx) {
    /*
     * Tracks the edges of the newly labelled S-vertices, for delta2 or delta3 depending on the
     * label at the other end. Tight edges are taken up later, by zero delta steps.
     */
    for (int64_t i = 0; i < ctx->scanQueueLength; i++) {
        int64_t x = ctx->scanQueue[i];
        int64_t bx = topLevelBlossom(ctx, x);
        assert(ctx->blossoms[bx].label == LABEL_S);
        for (int64_t j = ctx->adjacencyStarts[x]; j < ctx->adjacencyStarts[x + 1]; j++) {
            int64_t e = ctx->adjacentEdges[j];
            int64_t y = otherVertex(ctx, e, x);
            int64_t by = topLevelBlossom(ctx, y);
            if (bx == by) { //Internal to a blossom
                continue;
            }
            if (ctx->blossoms[by].label == LABEL_S) {
                delta3AddEdge(ctx, e);
            } else {
                delta2AddEdge(ctx, e, y, &ctx->blossoms[by]);
            }
        }
    }
    ctx->scanQueueLength = 0;
}

////////////////////////////////////
// Main algorithm
////////////////////////////////////

static bool runStage(MatchingContext *ctx) {
    /*
     * Runs one stage: delta steps until an augmenting path is found and used, returning true, or
     * until the first kind of step proves the matching optimal, returning false.
     */
    while (1) {
        scanNewSVertices(ctx);

        //Compute the delta step: the least of the four kinds, preferring the later kinds on ties
        int64_t deltaType = 1;
        int64_t delta = ctx->initVertexDual - ctx->deltaSum; //The dual of the unmatched vertices
        int64_t deltaEdge = NO_EDGE;
        Blossom *deltaBlossom = NULL;
        int64_t slack = 0;
        int64_t e = delta2GetMinEdge(ctx, &slack); //Least slack from an S-vertex to an unlabelled one
        if (e != NO_EDGE && slack <= delta) {
            deltaType = 2;
            delta = slack;
            deltaEdge = e;
        }
        e = delta3GetMinEdge(ctx, &slack); //Half the least slack between S-blossoms
        if (e != NO_EDGE && slack / 2 <= delta) {
            deltaType = 3;
            delta = slack / 2;
            deltaEdge = e;
        }
        if (!heap_empty(&ctx->delta4Queue)) { //Half the least dual of a T-blossom
            Blossom *blossom = &ctx->blossoms[heap_minElem(&ctx->delta4Queue)];
            assert(blossom->parent == NO_BLOSSOM);
            assert(blossom->label == LABEL_T);
            int64_t blossomDual = blossom->dualVar - 2 * ctx->deltaSum;
            if (blossomDual / 2 <= delta) {
                deltaType = 4;
                delta = blossomDual / 2;
                deltaBlossom = blossom;
            }
        }

        //Apply it, which implicitly updates the duals
        ctx->deltaSum += delta;

        if (deltaType == 2) { //Use the edge from an S-vertex to an unlabelled vertex it made tight
            int64_t x = ctx->edgeX[deltaEdge], y = ctx->edgeY[deltaEdge];
            if (ctx->blossoms[topLevelBlossom(ctx, x)].label != LABEL_S) {
                int64_t z = x;
                x = y;
                y = z;
            }
            extendTreeSToT(ctx, x, y);
        } else if (deltaType == 3) { //Use the edge between S-blossoms it made tight
            if (addSToSEdge(ctx, ctx->edgeX[deltaEdge], ctx->edgeY[deltaEdge])) {
                return 1;
            }
        } else if (deltaType == 4) { //Expand the T-blossom whose dual it made zero
            expandTBlossom(ctx, deltaBlossom);
        } else { //No further improvement is possible
            assert(deltaType == 1);
            return 0;
        }
    }
}

static void cleanup(MatchingContext *ctx) {
    /*
     * Unlabels the remaining alternating trees and applies the pending updates to the vertex
     * duals.
     */
    assert(ctx->scanQueueLength == 0);
    for (int64_t b = 0; b < 2 * ctx->vertexNumber; b++) {
        Blossom *blossom = &ctx->blossoms[b];
        if (blossom->isNonTrivial && blossom->subblossomNumber == 0) { //Unused slot
            continue;
        }
        if (blossom->parent == NO_BLOSSOM && blossom->label != LABEL_NONE) {
            resetBlossomLabel(ctx, blossom);
        }
        if (blossom->vertexDualOffset != 0) {
            int64_t dualFixup = blossom->vertexDualOffset;
            blossom->vertexDualOffset = 0;
            int64_t vertexNumber = getBlossomVertices(ctx, b);
            for (int64_t i = 0; i < vertexNumber; i++) {
                ctx->vertexDual[ctx->blossomVertices[i]] += dualFixup;
            }
        }
    }
    assert(heap_empty(&ctx->delta2Queue));
    assert(heap_empty(&ctx->delta3Queue));
    assert(heap_empty(&ctx->delta4Queue));
}

static void run(MatchingContext *ctx) {
    //Every vertex starts as an S-labelled root of its own alternating tree
    for (int64_t x = 0; x < ctx->vertexNumber; x++) {
        assert(ctx->vertexMate[x] == NO_VERTEX);
        Blossom *bx = &ctx->blossoms[topLevelBlossom(ctx, x)];
        assert(bx->baseVertex == x);
        assignBlossomLabelS(ctx, bx);
        bx->treeEdgeX = NO_VERTEX;
        bx->treeEdgeY = x;
        bx->treeRoot = x;
    }
    //Each stage but the last adds an edge to the matching
    while (runStage(ctx)) {
        ;
    }
    cleanup(ctx);
}

#ifndef NDEBUG

////////////////////////////////////
// Verification
////////////////////////////////////

/*
 * Checks that the matching is optimal, by checking that the duals are feasible and meet the
 * complementary slackness conditions with it: every dual is non-negative, every unmatched vertex
 * has zero dual, every blossom with a positive dual is full, every edge has non-negative slack and
 * every matched edge has zero slack. This takes time O(n^2), so it is only done in debug builds.
 */

static bool checkedAdd(int64_t *result, int64_t a, int64_t b) {
    //Returns true on overflow; only used on non-negative values
    if (a > INT64_MAX - b) {
        return 1;
    }
    *result = a + b;
    return 0;
}

static bool verifyBlossom(MatchingContext *ctx, int64_t topBlossom, int64_t *edgeDuals) {
    /*
     * Descends the blossom tree under a top-level blossom, adding to each edge the duals of the
     * blossoms holding it, and checks on the way back up that the blossoms are full.
     */
    int64_t n = ctx->vertexNumber;
    int64_t *vertexDepth = st_calloc(n, sizeof(int64_t)); //depth of the smallest blossom on the current path holding x
    int64_t *pathSumDual = st_calloc(n + 2, sizeof(int64_t)); //sum of the blossom duals down to each depth
    int64_t *pathMatchedNumber = st_calloc(n + 2, sizeof(int64_t)); //matched edges found at each depth
    int64_t *stackBlossoms = st_malloc((n + 1) * sizeof(int64_t));
    int64_t *stackPositions = st_malloc((n + 1) * sizeof(int64_t));
    int64_t *vertices = st_malloc(n * sizeof(int64_t));
    int64_t stackLength = 0;
    bool ok = 1;
    stackBlossoms[stackLength] = topBlossom;
    stackPositions[stackLength++] = 0;
    while (stackLength > 0 && ok) {
        int64_t depth = stackLength;
        Blossom *blossom = &ctx->blossoms[stackBlossoms[stackLength - 1]];
        int64_t position = stackPositions[stackLength - 1];
        if (position == 0) { //Entering the blossom
            int64_t vertexNumber = getBlossomVertices(ctx, stackBlossoms[stackLength - 1]);
            memcpy(vertices, ctx->blossomVertices, vertexNumber * sizeof(int64_t));
            for (int64_t i = 0; i < vertexNumber; i++) {
                vertexDepth[vertices[i]] = depth;
            }
            pathSumDual[depth] = pathSumDual[depth - 1];
            if (blossom->subblossomNumber < 3 || checkedAdd(&pathSumDual[depth], pathSumDual[depth], blossom->dualVar)) {
                ok = 0;
                break;
            }
            pathMatchedNumber[depth] = 0;
        }
        if (position < blossom->subblossomNumber) {
            stackPositions[stackLength - 1]++;
            Blossom *sub = &ctx->blossoms[blossom->subblossoms[position]];
            if (sub->isNonTrivial) {
                stackBlossoms[stackLength] = blossom->subblossoms[position];
                stackPositions[stackLength++] = 0;
            } else { //Find the smallest blossom holding each edge out of the vertex
                int64_t x = sub->baseVertex;
                for (int64_t j = ctx->adjacencyStarts[x]; j < ctx->adjacencyStarts[x + 1]; j++) {
                    int64_t e = ctx->adjacentEdges[j];
                    if (ctx->edgeX[e] == x) {
                        int64_t y = ctx->edgeY[e];
                        int64_t edgeDepth = vertexDepth[y];
                        if (edgeDepth > 0) {
                            if (checkedAdd(&edgeDuals[e], edgeDuals[e], pathSumDual[edgeDepth])) {
                                ok = 0;
                                break;
                            }
                            if (ctx->vertexMate[x] == y) {
                                pathMatchedNumber[edgeDepth]++;
                            }
                        }
                    }
                }
            }
        } else { //Leaving the blossom: check it is full
            int64_t vertexNumber = getBlossomVertices(ctx, stackBlossoms[stackLength - 1]);
            if (vertexNumber != 2 * pathMatchedNumber[depth] + 1) {
                ok = 0;
                break;
            }
            pathMatchedNumber[depth - 1] += pathMatchedNumber[depth];
            for (int64_t i = 0; i < vertexNumber; i++) {
                vertexDepth[ctx->blossomVertices[i]] = depth - 1;
            }
            stackLength--;
        }
    }
    free(vertexDepth);
    free(pathSumDual);
    free(pathMatchedNumber);
    free(stackBlossoms);
    free(stackPositions);
    free(vertices);
    return ok;
}

static bool verifyOptimum(MatchingContext *ctx) {
    int64_t n = ctx->vertexNumber, m = ctx->edgeNumber;
    //The mates are symmetric and each matched pair is an edge
    int64_t matchedVertexNumber = 0, matchedEdgeNumber = 0;
    for (int64_t x = 0; x < n; x++) {
        int64_t y = ctx->vertexMate[x];
        if (y != NO_VERTEX) {
            matchedVertexNumber++;
            if (ctx->vertexMate[y] != x) {
                return 0;
            }
        }
    }
    for (int64_t e = 0; e < m; e++) {
        if (ctx->vertexMate[ctx->edgeX[e]] == ctx->edgeY[e]) {
            matchedEdgeNumber++;
        }
    }
    if (matchedVertexNumber != 2 * matchedEdgeNumber) {
        return 0;
    }
    //Vertex duals are non-negative, and zero if unmatched; blossom duals are non-negative
    for (int64_t x = 0; x < n; x++) {
        if (ctx->vertexDual[x] < 0 || (ctx->vertexMate[x] == NO_VERTEX && ctx->vertexDual[x] != 0)) {
            return 0;
        }
    }
    for (int64_t b = n; b < 2 * n; b++) {
        if (ctx->blossoms[b].subblossomNumber > 0 && ctx->blossoms[b].dualVar < 0) {
            return 0;
        }
    }
    //Blossoms are full, and each edge's duals sum to at least twice its weight, exactly if matched
    int64_t *edgeDuals = st_calloc(m > 0 ? m : 1, sizeof(int64_t));
    bool ok = 1;
    for (int64_t e = 0; e < m && ok; e++) {
        ok = !checkedAdd(&edgeDuals[e], ctx->vertexDual[ctx->edgeX[e]], ctx->vertexDual[ctx->edgeY[e]]);
    }
    for (int64_t b = n; b < 2 * n && ok; b++) {
        if (ctx->blossoms[b].subblossomNumber > 0 && ctx->blossoms[b].parent == NO_BLOSSOM) {
            ok = verifyBlossom(ctx, b, edgeDuals);
        }
    }
    for (int64_t e = 0; e < m && ok; e++) {
        int64_t weight = 2 * ctx->edgeWeight[e];
        if (weight > edgeDuals[e] || (ctx->vertexMate[ctx->edgeX[e]] == ctx->edgeY[e] && edgeDuals[e] != weight)) {
            ok = 0;
        }
    }
    free(edgeDuals);
    return ok;
}

#endif

////////////////////////////////////
// Interface
////////////////////////////////////

static void computeMaximumWeightMatching(int64_t vertexNumber, int64_t edgeNumber,
        const int64_t *edgeX, const int64_t *edgeY, const int64_t *edgeWeight, int64_t *mates) {
    /*
     * Fills mates with the vertex each vertex is matched to, or -1, in a maximum weight matching of
     * the graph of edges (edgeX[e], edgeY[e]) with weights edgeWeight[e] >= 0.
     */
    MatchingContext ctx;
    matchingContext_construct(&ctx, vertexNumber, edgeNumber, edgeX, edgeY, edgeWeight);
    run(&ctx);
#ifndef NDEBUG
    if (!verifyOptimum(&ctx)) {
        st_errAbort("The maximum weight matching failed to verify as optimal");
    }
#endif
    memcpy(mates, ctx.vertexMate, vertexNumber * sizeof(int64_t));
    matchingContext_destruct(&ctx);
}

static int64_t findComponent(int64_t *parents, int64_t x) {
    while (parents[x] != x) {
        parents[x] = parents[parents[x]];
        x = parents[x];
    }
    return x;
}

stList *chooseMatching_edmonds(stList *edges, int64_t nodeNumber) {
    /*
     * Edges with weight <= 0 never add to the weight of a matching, so they are left out, and
     * each connected component of what remains is matched separately: the union of maximum
     * weight matchings of the components is one of the whole graph. A component that is a single
     * edge, the common case, is simply matched.
     */
    int64_t edgeNumber = 0;
    for (int64_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        int64_t x = stIntTuple_get(edge, 0), y = stIntTuple_get(edge, 1);
        if (x < 0 || x >= nodeNumber || y < 0 || y >= nodeNumber || x == y) {
            st_errAbort("Matching edge (%" PRIi64 ", %" PRIi64 ") is not between two distinct nodes in [0, %" PRIi64 ")",
                    x, y, nodeNumber);
        }
        if (stIntTuple_get(edge, 2) > 0) {
            edgeNumber++;
        }
    }
    stList *matching = stList_construct();
    if (edgeNumber == 0) {
        return matching;
    }

    //Gather the edges with positive weight, and join their ends into components
    stIntTuple **positiveEdges = st_malloc(edgeNumber * sizeof(stIntTuple *));
    int64_t *parents = st_malloc(nodeNumber * sizeof(int64_t));
    for (int64_t x = 0; x < nodeNumber; x++) {
        parents[x] = x;
    }
    int64_t j = 0;
    for (int64_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        if (stIntTuple_get(edge, 2) > 0) {
            positiveEdges[j++] = edge;
            int64_t cx = findComponent(parents, stIntTuple_get(edge, 0));
            int64_t cy = findComponent(parents, stIntTuple_get(edge, 1));
            if (cx != cy) {
                parents[cx > cy ? cx : cy] = cx < cy ? cx : cy;
            }
        }
    }

    //Bucket the edges by component, keeping their order, and number the nodes of each component
    //from 0 in increasing order
    int64_t *componentEdgeStarts = st_calloc(nodeNumber + 1, sizeof(int64_t));
    for (int64_t e = 0; e < edgeNumber; e++) {
        componentEdgeStarts[findComponent(parents, stIntTuple_get(positiveEdges[e], 0)) + 1]++;
    }
    for (int64_t x = 0; x < nodeNumber; x++) {
        componentEdgeStarts[x + 1] += componentEdgeStarts[x];
    }
    int64_t *componentEdges = st_malloc(edgeNumber * sizeof(int64_t));
    int64_t *fill = st_malloc(nodeNumber * sizeof(int64_t));
    memcpy(fill, componentEdgeStarts, nodeNumber * sizeof(int64_t));
    for (int64_t e = 0; e < edgeNumber; e++) {
        componentEdges[fill[findComponent(parents, stIntTuple_get(positiveEdges[e], 0))]++] = e;
    }
    int64_t *localNodes = st_malloc(nodeNumber * sizeof(int64_t));
    int64_t *componentSizes = st_calloc(nodeNumber, sizeof(int64_t));
    for (int64_t x = 0; x < nodeNumber; x++) {
        int64_t c = findComponent(parents, x);
        localNodes[x] = componentSizes[c]++;
    }

    //Match each component, recording for each node its matched edge
    stIntTuple **matchedEdges = st_calloc(nodeNumber, sizeof(stIntTuple *));
    int64_t *localX = st_malloc(edgeNumber * sizeof(int64_t));
    int64_t *localY = st_malloc(edgeNumber * sizeof(int64_t));
    int64_t *localWeight = st_malloc(edgeNumber * sizeof(int64_t));
    int64_t *localToGlobal = st_malloc(nodeNumber * sizeof(int64_t));
    int64_t *mates = st_malloc(nodeNumber * sizeof(int64_t));
    for (int64_t c = 0; c < nodeNumber; c++) {
        int64_t start = componentEdgeStarts[c], componentEdgeNumber = componentEdgeStarts[c + 1] - start;
        if (componentEdgeNumber == 0) {
            continue;
        }
        if (componentEdgeNumber == 1) {
            stIntTuple *edge = positiveEdges[componentEdges[start]];
            matchedEdges[stIntTuple_get(edge, 0)] = edge;
            matchedEdges[stIntTuple_get(edge, 1)] = edge;
            continue;
        }
        int64_t maxWeight = 0;
        for (int64_t i = 0; i < componentEdgeNumber; i++) {
            stIntTuple *edge = positiveEdges[componentEdges[start + i]];
            int64_t x = stIntTuple_get(edge, 0), y = stIntTuple_get(edge, 1);
            localX[i] = localNodes[x];
            localY[i] = localNodes[y];
            localWeight[i] = stIntTuple_get(edge, 2);
            localToGlobal[localX[i]] = x;
            localToGlobal[localY[i]] = y;
            if (localWeight[i] > maxWeight) {
                maxWeight = localWeight[i];
            }
        }
        if (maxWeight > INT64_MAX / 6) {
            //Too heavy to match exactly, so divide the weights by the least factor that brings them
            //under the limit, rounding up to keep them positive. The factor is at most 7, so this
            //only loses the last three bits of weights that are above 2^60.
            int64_t factor = maxWeight / (INT64_MAX / 6) + 1;
            st_logInfo("Dividing the matching weights of a component of %" PRIi64 " edges by %" PRIi64
                    ", as its heaviest edge weighs %" PRIi64 "\n", componentEdgeNumber, factor, maxWeight);
            for (int64_t i = 0; i < componentEdgeNumber; i++) {
                localWeight[i] = localWeight[i] / factor + (localWeight[i] % factor != 0);
            }
        }
        int64_t componentNodeNumber = componentSizes[c];
        computeMaximumWeightMatching(componentNodeNumber, componentEdgeNumber, localX, localY, localWeight, mates);
        //Each matched pair is joined by one edge, which there must be only one of
        for (int64_t i = 0; i < componentEdgeNumber; i++) {
            if (mates[localX[i]] == localY[i]) {
                int64_t x = localToGlobal[localX[i]], y = localToGlobal[localY[i]];
                if (matchedEdges[x] != NULL) {
                    st_errAbort("Matching edges must be unique, but (%" PRIi64 ", %" PRIi64 ") is given twice", x, y);
                }
                matchedEdges[x] = matchedEdges[y] = positiveEdges[componentEdges[start + i]];
            }
        }
    }

    //Report the matched edges in order of their lesser node, as blossom5 did
    for (int64_t x = 0; x < nodeNumber; x++) {
        stIntTuple *edge = matchedEdges[x];
        if (edge != NULL && x < stIntTuple_get(edge, 0) + stIntTuple_get(edge, 1) - x) { //x is the lesser of its two nodes
            stList_append(matching, edge);
        }
    }
    st_logDebug("The maximum weight matching for %" PRIi64 " nodes with %" PRIi64 " edges of positive weight has %" PRIi64 " edges\n",
            nodeNumber, edgeNumber, stList_length(matching));

    free(positiveEdges);
    free(parents);
    free(componentEdgeStarts);
    free(componentEdges);
    free(fill);
    free(localNodes);
    free(componentSizes);
    free(matchedEdges);
    free(localX);
    free(localY);
    free(localWeight);
    free(localToGlobal);
    free(mates);
    return matching;
}
