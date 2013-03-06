/**
 * dijkstra.cpp
 *
 * [D, P] = DIJKSTRA(G, S)
 *
 *   This function computes the shortest distance tree(s) from one (or more)
 *   source node(s) in a graph.
 *
 *   G is a sparse matrix where element (i,j) contains the cost of going
 *   from node i to node j in the graph. That is, G contains the edge
 *   weights or distances between nodes.
 *
 *   Note: Matlab's convention is that the absence of an element in a sparse
 *   matrix means edge weight = 0. However, in this function, the lack of an
 *   entry in the sparse matrix is understood as a lack of edge (or edge
 *   weight = Inf).
 *
 *   S is a vector with a list of source node indices. A shortest path
 *   tree will be computed for each one of the nodes in S.
 *
 *   D, P are matrices where each row corresponds to a source node in S. D
 *   is the shortest distance from each node to the source node.
 *
 *   P is the predecessor of each node in the shortest path tree.
 *
 *   Because this function is implemented as a MEX function in C++,
 *   takes a sparse matrix at the input, and uses a Fibonacci heap
 *   implementation [1], it is well suited for huge sparse graphs.
 *
 *   [1]
 *   http://www.ahhf45.com/info/Data_Structures_and_Algorithms/resources/technical_artile/fibonacci_heap/fibonacci.htm
 *
 * Author: Mark Steyvers, Stanford University, 19 Dec 2000.
 * Version: 0.2.0
 * $Rev$
 * $Date$
 *
 * Modified by Ramon Casero <rcasero@gmail.com>, University of Oxford,
 * 23 Mar 2010 to also provide the predecessor list at the output.
 * 6  Feb 2013 to accept a list of targets. The algorithm will stop after hitting all targets
 *
 * This file is distributed as a derivative work of a third-party function
 * with project Gerardus.
 *
 * http://code.google.com/p/gerardus/
 *
 *
 **/

#include <math.h>
#include "mex.h"
#include <list>

extern void _main();

#include <stdlib.h>
// #include <iostream.h>
#include <iostream>
#include <stdio.h>
// #include <conio.h>
#include <ctype.h>
#include <memory.h>
#include <time.h>

#include "fibheap.h"

#define _FIBHEAP_CPP

//***************************************************************************
// This Fibonacci heap implementation is Copyright (c) 1996 by John Boyer.
// See the header file for free usage information.
//***************************************************************************

//***************************************************************************
// The classes in this package are designed to allow the package user
// to quickly and easily develop applications that require a heap data
// structure.  Using amortized analysis, the asymptotically fastest heap
// data structure is the Fibonacci heap.  The constants are a little
// high so the real speed gain will not be seen until larger data sets
// are required, but in most cases, if the data set is small, then the
// run-time will be neglible anyway.
//
// To use this heap class you need do only two things.  First, subclass
// the FibHeapNode class to create the class of objects that you'd
// like to store in a heap.  Second, create an instance of the FibHeap
// class, which can then be used to Insert(), ExtractMin(), etc.,
// instances of your FibHeapNode subclass.  Notice that you don't need
// to create a subclass of the FibHeap class.
//
// The application-specific data object that you'd like to store in a heap
// will have a key value.  In the class that you derive from FibHeapNode,
// you will need to define the key structure then provide assignment (=),
// equality (==) and less-than operators and a destructor.  These functions
// are declared virtual so that the code in the FibHeap class can compare,
// assign and destroy your objects by calling on your code.
//
// The overloaded operators in your defined class MUST call functions in
// the Fibonacci heap node class first.  For assignment, the function
// FHN_Assign() should be called before code that deals with the copy of
// the key value.  For comparison operators, the function FHN_Cmp() should
// appear first.  If it returns 0, then keys can be compared as normal.
// The following indicates what the three most common operators must do
// based on the return value of FHN_Cmp() 
//
// For ==, if zero returned, then compare keys
//	   if non-zero X returned, then return 0
// For <,  if zero returned, then compare keys
//         if non-zero X returned, then return X<0?1:0
// For >,  if zero returned, then compare keys
//         if non-zero X returned, then return X>0?1:0   
//***************************************************************************

#include <stdlib.h>
#include <iostream>
#include <stdio.h>

#include "fibheap.h"

//***************************************************************************
//=========================================================
// FibHeapNode Constructor
//=========================================================
//***************************************************************************

FibHeapNode::FibHeapNode()
{
     Left = Right = Parent = Child = NULL;
     Degree = Mark = NegInfinityFlag = 0;
}

//=========================================================
// FibHeapNode Destructor
//
// Body is empty, but declaration is required in order to
// force virtual.  This will ensure that FibHeap class
// calls derived class destructors.
//=========================================================

FibHeapNode::~FibHeapNode()
{
}

//=========================================================
// FHN_Assign()
//
// To be used as first step of an assignment operator in a
// derived class.  The derived class will handle assignment
// of key value, and this function handles copy of the
// NegInfinityFlag (which overrides the key value if it is
// set).
//=========================================================

void FibHeapNode::FHN_Assign(FibHeapNode& RHS)
{
     NegInfinityFlag = RHS.NegInfinityFlag;
}

//=========================================================
// FHN_Cmp()
//
// To be used as the first step of ALL comparators in a
// derived class.
//
// Compares the relative state of the two neg. infinity
// flags.  Note that 'this' is the left hand side.  If
// LHS neg. infinity is set, then it will be less than (-1)
// the RHS unless RHS neg. infinity flag is also set.
// Only if function returns 0 should the key comparison
// defined in the derived class be performed, e.g.
//
// For ==, if zero returned, then compare keys
//	   if non-zero X returned, then return 0
// For <,  if zero returned, then compare keys
//         if non-zero X returned, then return X<0?1:0
// For >,  if zero returned, then compare keys
//         if non-zero X returned, then return X>0?1:0    
//=========================================================

int  FibHeapNode::FHN_Cmp(FibHeapNode& RHS)
{
     if (NegInfinityFlag)
	 return RHS.NegInfinityFlag ? 0 : -1;
     return RHS.NegInfinityFlag ? 1 : 0; 
}

//========================================================================
// We do, on occasion, compare and assign objects of type FibHeapNode, but
// only when the NegInfinityFlag is set.  See for example FibHeap::Delete().
//
// Also, these functions exemplify what a derived class should do.
//========================================================================

void FibHeapNode::operator =(FibHeapNode& RHS)
{
     FHN_Assign(RHS);
     // Key assignment goes here in derived classes
}

int  FibHeapNode::operator ==(FibHeapNode& RHS)
{
     if (FHN_Cmp(RHS)) return 0;
     // Key compare goes here in derived classes
     return 1;
}

int  FibHeapNode::operator <(FibHeapNode& RHS)
{
int X;

     if ((X=FHN_Cmp(RHS)) != 0)
	  return X < 0 ? 1 : 0;
     // Key compare goes here in derived classes
     return 0;
}

//=========================================================
// Print()
//=========================================================

void FibHeapNode::Print()
{
     if (NegInfinityFlag)
       std::cout << "-inf.";
}

//***************************************************************************
//===========================================================================
// FibHeap Constructor
//===========================================================================
//***************************************************************************

FibHeap::FibHeap()
{
     MinRoot = NULL;
     NumNodes = NumTrees = NumMarkedNodes = 0;
     ClearHeapOwnership();
}

//===========================================================================
// FibHeap Destructor
//===========================================================================

FibHeap::~FibHeap()
{
FibHeapNode *Temp;

     if (GetHeapOwnership())
     {
         while (MinRoot != NULL)
         {
             Temp = ExtractMin();
             delete Temp;
	 }
     }
}

//===========================================================================
// Insert() - O(1) actual; O(2) amortized
//
// I am using O(2) here to indicate that although Insert() is
// constant time, its amortized rating is more costly because some
// of the work of inserting is done by other operations such as
// ExtractMin(), which is where tree-balancing occurs.
//
// The child pointer is deliberately not set to NULL because Insert()
// is also used internally to help put whole trees onto the root list.
//===========================================================================

void FibHeap::Insert(FibHeapNode *NewNode)
{
     if (NewNode == NULL) return;

// If the heap is currently empty, then new node becomes singleton
// circular root list
 
     if (MinRoot == NULL)
	 MinRoot = NewNode->Left = NewNode->Right = NewNode;

     else
     {
// Pointers from NewNode set to insert between MinRoot and MinRoot->Right

         NewNode->Right = MinRoot->Right;
	 NewNode->Left = MinRoot;

// Set Pointers to NewNode  

	 NewNode->Left->Right = NewNode;
         NewNode->Right->Left = NewNode;

// The new node becomes new MinRoot if it is less than current MinRoot

         if (*NewNode < *MinRoot)
	     MinRoot = NewNode;
     }

// We have one more node in the heap, and it is a tree on the root list

     NumNodes++;

     NumTrees++;
     NewNode->Parent = NULL;
}

//===========================================================================
// Union() - O(1) actual; O(1) amortized
//===========================================================================

void FibHeap::Union(FibHeap *OtherHeap)
{
FibHeapNode *Min1, *Min2, *Next1, *Next2;

     if (OtherHeap == NULL || OtherHeap->MinRoot == NULL) return;

// We join the two circular lists by cutting each list between its
// min node and the node after the min.  This code just pulls those
// nodes into temporary variables so we don't get lost as changes
// are made.

     Min1 = MinRoot;
     Min2 = OtherHeap->MinRoot;
     Next1 = Min1->Right;
     Next2 = Min2->Right;

// To join the two circles, we join the minimum nodes to the next
// nodes on the opposite chains.  Conceptually, it looks like the way
// two bubbles join to form one larger bubble.  They meet at one point
// of contact, then expand out to make the bigger circle.
 
     Min1->Right = Next2;
     Next2->Left = Min1;
     Min2->Right = Next1;
     Next1->Left = Min2;

// Choose the new minimum for the heap
 
     if (*Min2 < *Min1)
         MinRoot = Min2;

// Set the amortized analysis statistics and size of the new heap
                   
     NumNodes += OtherHeap->NumNodes;
     NumMarkedNodes += OtherHeap->NumMarkedNodes;
     NumTrees += OtherHeap->NumTrees;

// Complete the union by setting the other heap to emptiness
// then destroying it

     OtherHeap->MinRoot  = NULL;
     OtherHeap->NumNodes =
     OtherHeap->NumTrees =
     OtherHeap->NumMarkedNodes = 0;

     delete OtherHeap;
}

//===========================================================================
// Minimum - O(1) actual; O(1) amortized
//===========================================================================

FibHeapNode *FibHeap::Minimum()
{
     return MinRoot;
}

//===========================================================================
// ExtractMin() - O(n) worst-case actual; O(lg n) amortized
//===========================================================================

FibHeapNode *FibHeap::ExtractMin()
{
FibHeapNode *Result;
FibHeap *ChildHeap = NULL;

// Remove minimum node and set MinRoot to next node

     if ((Result = Minimum()) == NULL)
          return NULL;

     MinRoot = Result->Right;
     Result->Right->Left = Result->Left;
     Result->Left->Right = Result->Right;
     Result->Left = Result->Right = NULL;

     NumNodes --;
     if (Result->Mark)
     {
	 NumMarkedNodes --;
         Result->Mark = 0;
     }
     Result->Degree = 0;

// Attach child list of Minimum node to the root list of the heap
// If there is no child list, then do no work

     if (Result->Child == NULL)
     {
	 if (MinRoot == Result)
	     MinRoot = NULL;
     }

// If MinRoot==Result then there was only one root tree, so the
// root list is simply the child list of that node (which is
// NULL if this is the last node in the list)

     else if (MinRoot == Result)
         MinRoot = Result->Child;

// If MinRoot is different, then the child list is pushed into a
// new temporary heap, which is then merged by Union() onto the
// root list of this heap.

     else 
     {
         ChildHeap = new FibHeap();
         ChildHeap->MinRoot = Result->Child;
     }

// Complete the disassociation of the Result node from the heap

     if (Result->Child != NULL)
	 Result->Child->Parent = NULL;
     Result->Child = Result->Parent = NULL;

// If there was a child list, then we now merge it with the
//	rest of the root list

     if (ChildHeap)
         Union(ChildHeap);

// Consolidate heap to find new minimum and do reorganize work

     if (MinRoot != NULL)
         _Consolidate();

// Return the minimum node, which is now disassociated with the heap
// It has Left, Right, Parent, Child, Mark and Degree cleared.

     return Result;
}

//===========================================================================
// DecreaseKey() - O(lg n) actual; O(1) amortized
//
// The O(lg n) actual cost stems from the fact that the depth, and
// therefore the number of ancestor parents, is bounded by O(lg n).
//===========================================================================

int  FibHeap::DecreaseKey(FibHeapNode *theNode, FibHeapNode& NewKey)
{
FibHeapNode *theParent;

     if (theNode==NULL || *theNode < NewKey)
         return NOTOK;

     *theNode = NewKey;

     theParent = theNode->Parent;
     if (theParent != NULL && *theNode < *theParent)
     {
         _Cut(theNode, theParent);
         _CascadingCut(theParent);
     }

     if (*theNode < *MinRoot)
         MinRoot = theNode;

     return OK;
}

//===========================================================================
// Delete() - O(lg n) amortized; ExtractMin() dominates
//
// Notice that if we don't own the heap nodes, then we clear the
// NegInfinityFlag on the deleted node.  Presumably, the programmer
// will be reusing the node.
//===========================================================================

int  FibHeap::Delete(FibHeapNode *theNode)
{
FibHeapNode Temp;
int Result;

     if (theNode == NULL) return NOTOK;

     Temp.NegInfinityFlag = 1;
     Result = DecreaseKey(theNode, Temp);

     if (Result == OK)
         if (ExtractMin() == NULL)
             Result = NOTOK;

     if (Result == OK)
     {
         if (GetHeapOwnership())
	      delete theNode;
	 else theNode->NegInfinityFlag = 0;
     }
         
     return Result;
}

//===========================================================================

void FibHeap::_Exchange(FibHeapNode*& N1, FibHeapNode*& N2)
{
FibHeapNode *Temp;

     Temp = N1;
     N1 = N2;
     N2 = Temp;
}

//===========================================================================
// Consolidate()
//
// Internal function that reorganizes heap as part of an ExtractMin().
// We must find new minimum in heap, which could be anywhere along the
// root list.  The search could be O(n) since all nodes could be on
// the root list.  So, we reorganize the tree into more of a binomial forest
// structure, and then find the new minimum on the consolidated O(lg n) sized
// root list, and in the process set each Parent pointer to NULL, and count
// the number of resulting subtrees.
//
// Note that after a list of n inserts, there will be n nodes on the root
// list, so the first ExtractMin() will be O(n) regardless of whether or
// not we consolidate.  However, the consolidation causes subsequent
// ExtractMin() operations to be O(lg n).  Furthermore, the extra cost of
// the first ExtractMin() is covered by the higher amortized cost of
// Insert(), which is the real governing factor in how costly the first
// ExtractMin() will be.
//===========================================================================

void FibHeap::_Consolidate()
{
FibHeapNode *x, *y, *w;
FibHeapNode *A[1+8*sizeof(long)]; // 1+lg(n)
int  I=0, Dn = 1+8*sizeof(long);
short d;

// Initialize the consolidation detection array

     for (I=0; I < Dn; I++)
          A[I] = NULL;

// We need to loop through all elements on root list.
// When a collision of degree is found, the two trees
// are consolidated in favor of the one with the lesser
// element key value.  We first need to break the circle
// so that we can have a stopping condition (we can't go
// around until we reach the tree we started with
// because all root trees are subject to becoming a
// child during the consolidation).

     MinRoot->Left->Right = NULL;
     MinRoot->Left = NULL;
     w = MinRoot;

     do {
//cout << "Top of Consolidate's loop\n";
//Print(w);

	x = w;
        d = x->Degree;
        w = w->Right;

// We need another loop here because the consolidated result
// may collide with another large tree on the root list.

        while (A[d] != NULL)
        {
             y = A[d];
	     if (*y < *x)
		 _Exchange(x, y);
             if (w == y) w = y->Right;
             _Link(y, x);
             A[d] = NULL;
             d++;
//cout << "After a round of Linking\n";
//Print(x);
	}
	A[d] = x;

     } while (w != NULL);

// Now we rebuild the root list, find the new minimum,
// set all root list nodes' parent pointers to NULL and
// count the number of subtrees.

     MinRoot = NULL;
     NumTrees = 0;
     for (I = 0; I < Dn; I++)
          if (A[I] != NULL)
              _AddToRootList(A[I]);
}

//===========================================================================
// The node y is removed from the root list and becomes a subtree of node x.
//===========================================================================

void FibHeap::_Link(FibHeapNode *y, FibHeapNode *x)
{
// Remove node y from root list

     if (y->Right != NULL)
	 y->Right->Left = y->Left;
     if (y->Left != NULL)
         y->Left->Right = y->Right;
     NumTrees--;

// Make node y a singleton circular list with a parent of x

     y->Left = y->Right = y;
     y->Parent = x;

// If node x has no children, then list y is its new child list

     if (x->Child == NULL)
	 x->Child = y;

// Otherwise, node y must be added to node x's child list

     else
     {
         y->Left = x->Child;
         y->Right = x->Child->Right;
         x->Child->Right = y;
         y->Right->Left = y;
     }

// Increase the degree of node x because it's now a bigger tree

     x->Degree ++;

// Node y has just been made a child, so clear its mark

     if (y->Mark) NumMarkedNodes--;
     y->Mark = 0;
}

//===========================================================================
//===========================================================================

void FibHeap::_AddToRootList(FibHeapNode *x)
{
     if (x->Mark) NumMarkedNodes --;
     x->Mark = 0;

     NumNodes--;
     Insert(x);
}

//===========================================================================
// Remove node x from the child list of its parent node y
//===========================================================================

void FibHeap::_Cut(FibHeapNode *x, FibHeapNode *y)
{
     if (y->Child == x)
         y->Child = x->Right;
     if (y->Child == x)
	 y->Child = NULL;

     y->Degree --;

     x->Left->Right = x->Right;
     x->Right->Left = x->Left;

     _AddToRootList(x);
}

//===========================================================================
// Cuts each node in parent list, putting successive ancestor nodes on the
// root list until we either arrive at the root list or until we find an
// ancestor that is unmarked.  When a mark is set (which only happens during
// a cascading cut), it means that one child subtree has been lost; if a
// second subtree is lost later during another cascading cut, then we move
// the node to the root list so that it can be re-balanced on the next
// consolidate. 
//===========================================================================

void FibHeap::_CascadingCut(FibHeapNode *y)
{
FibHeapNode *z = y->Parent;

     while (z != NULL)
     {
         if (y->Mark == 0)
         {
             y->Mark = 1;
             NumMarkedNodes++;
             z = NULL;
         }
         else
         {
             _Cut(y, z);
             y = z;
	     z = y->Parent;
         }
     }
}


class HeapNode : public FibHeapNode
{
      double   N;
      long int IndexV;

public:

      HeapNode() : FibHeapNode() { N = 0; };   

      virtual void operator =(FibHeapNode& RHS);
      virtual int  operator ==(FibHeapNode& RHS);
      virtual int  operator <(FibHeapNode& RHS);

      virtual void operator =(double NewKeyVal);
      virtual void Print();

      double GetKeyValue() { return N; };       /* !!!! */
      void SetKeyValue(double n) { N = n; };

      long int GetIndexValue() { return IndexV; };
      void SetIndexValue(long int v) { IndexV = v; };

};

void HeapNode::Print()
{
     FibHeapNode::Print();
     mexPrintf("%f (%d)" , N , IndexV);
}

void HeapNode::operator =(double NewKeyVal)
{
     HeapNode Temp;
     Temp.N = N = NewKeyVal;
     FHN_Assign(Temp);
}

void HeapNode::operator =(FibHeapNode& RHS)
{
     FHN_Assign(RHS);
     N = ((HeapNode&) RHS).N;
}

int  HeapNode::operator ==(FibHeapNode& RHS)
{
     if (FHN_Cmp(RHS)) return 0;
     return N == ((HeapNode&) RHS).N ? 1 : 0;
}

int  HeapNode::operator <(FibHeapNode& RHS)
{
int X;

     if ((X=FHN_Cmp(RHS)) != 0)
	  return X < 0 ? 1 : 0;

     return N < ((HeapNode&) RHS).N ? 1 : 0;
};

int IntCmp(const void *pA, const void *pB)
{
int A, B;

    A = *((const int *) pA);
    B = *((const int *) pB);
    if (A < B) return -1;
    if (A == B) return 0;
    return 1; 
}

// Note: "target": we call by value because we want to make a local
// copy of the list of targets that we can depopulate as we hit
// targets
void dodijk_sparse(long int M,
		   long int N,
		   long int S,
		   double *P, // parents
		   double   *D, // distances
		   double   *sr,
		   mwSize      *irs,
		   mwSize      *jcs,
		   HeapNode *A,
		   FibHeap  *theHeap,
		   std::list<mwIndex> target) {

  bool     finished;
  long int i,startind,endind,whichneighbor,ndone,closest;
  double   closestD,arclength; 
  double   INF,SMALL,olddist;
  HeapNode *Min;
  HeapNode Temp;
  
  INF   = mxGetInf();
  SMALL = mxGetEps();
  
  /* initialize */
  for (i=0; i<M; i++) {
    if (i!=S) A[i] = (double) INF; else A[i] = (double) SMALL;
    if (i!=S) D[i] = (double) INF; else D[i] = (double) SMALL;
    theHeap->Insert(&A[i]);
    A[i].SetIndexValue((long int) i);
    P[i] = mxGetNaN();
  }
  
  // Insert 0 then extract it. This will cause the
  // Fibonacci heap to get balanced.
  
  theHeap->Insert(&Temp);
  theHeap->ExtractMin();
  
  // theHeap->Print();
  // for (i=0; i<M; i++)
  // {
  //    closest = A[i].GetIndexValue();
  //    closestD = A[i].GetKeyValue();
  //    mexPrintf("Index at i=%d =%d  value=%f\n" , i , closest , closestD);
  // }

  // has the user provided a list of targets, so that the algorithm
  // will stop when all targets are hit?
  bool withTargets = target.size() > 0;

  /* loop over nonreached nodes */
  finished = false;
  ndone    = 0;
  while ((finished==false) && (ndone < M)) {
    
    // if ((ndone % 100) == 0) mexPrintf("Done with node %d\n" , ndone);
    
    Min = (HeapNode *) theHeap->ExtractMin();
    closest  = Min->GetIndexValue();
    closestD = Min->GetKeyValue();
    
    if ((closest<0) || (closest>=M)) mexErrMsgTxt("Minimum Index out of bounds...");
    
    // theHeap->Print();
    // mexPrintf("EXTRACTED MINIMUM  NDone=%d S=%d closest=%d closestD=%f\n" , ndone , S , closest , closestD);//TT
    // mexErrMsgTxt("Exiting...");
     
    D[closest] = closestD;
    
    // the rest of the graph is not reachable from here, or we have
    // hit all the targets, we have finished
    if (closestD == INF) { 
      
      finished = true; 
      
    } else {

      /* add the closest to the determined list */
      ndone++;

      // if the user provided a list the targets
      if (withTargets) {

	// drop this node from the list, if it's in it (note, the user
	// provides nodes as 1,...,N, while C++ stores them as
	// 0,...,N-1
	target.remove(closest+1);

	// if we have hit all the targets, this is the last iteration
	if (target.size() == 0) {
	  finished = true;
	}

      }
      

      /* relax all nodes adjacent to closest */
      startind = jcs[closest];
      endind   = jcs[closest+1] - 1;
      
      if (startind!=endind+1)
	for (i=startind; i<=endind; i++) {
	  whichneighbor = irs[i];
	  arclength = sr[i];
	  olddist   = D[whichneighbor];
	   
	  // mexPrintf("INSPECT NEIGHBOR #%d  olddist=%f newdist=%f\n" , whichneighbor , olddist , closestD+arclength);//TT
	  
	  if (olddist > (closestD + arclength)) {
	    D[whichneighbor] = closestD + arclength;
	    P[whichneighbor] = closest;
	     
	    Temp = A[whichneighbor];
	    Temp.SetKeyValue(closestD + arclength);
	    theHeap->DecreaseKey(&A[whichneighbor], Temp);
	    
	    // mexPrintf("UPDATING NODE #%d  olddist=%f newdist=%f newpred=%d\n", 
	    // 		  whichneighbor, olddist, closestD+arclength, closest);//TT
	    
	  } // if (olddist > (closestD + arclength))
	} // for (i=startind; i<=endind; i++)
      
    } // if (closestD == INF)
    
  } // while ((finished==false) && (ndone < M))

  // source node has no parent
  P[S] = -1;
  
}


//========================================================================
// Print()
//
// Used internally for debugging purposes.  The function prints the key
// value for each node along the root list, then it calls itself on each
// child list.   
//========================================================================

void FibHeap::Print(FibHeapNode *Tree, FibHeapNode *theParent)
{
FibHeapNode* Temp = NULL;


 if (Tree == NULL) {
   Tree = MinRoot;
 }

     Temp = Tree;
     do {
	if (Temp->Left == NULL)
	    mexPrintf("(Left is NULL)");
	Temp->Print();
	if (Temp->Parent != theParent)
	    mexPrintf("(Parent is incorrect)");
	if (Temp->Right == NULL)
	     mexPrintf("(Right is NULL)");
	else if (Temp->Right->Left != Temp)
	     mexPrintf("(Error in left link left) ->");
	else mexPrintf(" <-> ");

	Temp = Temp->Right;

// 	if (kbhit() && getch() == 27)
// 	{
// 	  std::cout << "Hit a key to resume or ESC to break\n";
// 	    if (getch() == 27)
//                 break;
//         }
     } while (Temp != NULL && Temp != Tree);
     mexPrintf("\n");

     Temp = Tree;
     do {
	mexPrintf("Children of ");
	Temp->Print();
	mexPrintf(": ");
	if (Temp->Child == NULL)
	     mexPrintf("NONE\n");
        else Print(Temp->Child, Temp);
	Temp = Temp->Right;
     } while (Temp!=NULL && Temp != Tree);

     if (theParent == NULL)
     {
     char ch;

	 mexPrintf("\n\n\n");
	 std::cin >> ch;
     }
}

//===========================================================================

void mexFunction(int          nlhs,
		 mxArray      *plhs[],
		 int          nrhs,
		 const mxArray *prhs[]) {

  double    *sr,*D,*P,*SS,*Dsmall;
  double *Psmall;
  mwIndex       *irs,*jcs;
  mwSize  M,N,S,MS,NS,i,j;
  
  HeapNode *A = NULL;
  FibHeap  *theHeap = NULL;

  // check input arguments
  if ((nrhs < 2) || (nrhs > 3)) {
    mexErrMsgTxt("Incorrect number of input arguments");
  } else if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments");
  }

  if (!mxIsSparse(prhs[0])) {
    mexErrMsgTxt("Function only implemented for sparse arrays");
  } 

  // input distance matrix dimensions
  M = mxGetM(prhs[0]);
  N = mxGetN(prhs[0]);
  
  if (M != N) mexErrMsgTxt("Input matrix needs to be square");
  
  // list of source nodes
  SS = mxGetPr(prhs[1]);
  MS = mxGetM(prhs[1]);
  NS = mxGetN(prhs[1]);
  
  if ((MS==0) || (NS==0) || ((MS>1) && (NS>1))) {
    mexErrMsgTxt("Source nodes are specified in one dimensional matrix only");
  }
  if (NS>MS) {
    MS=NS;
  }
  
  // distance values output
  plhs[0] = mxCreateDoubleMatrix(MS,M, mxREAL);
  D = mxGetPr(plhs[0]);
  
  // predecessors output
  plhs[1] = mxCreateDoubleMatrix(MS,M, mxREAL);
  P = mxGetPr(plhs[1]);

  // to simplify things, we expect the list of targets to be of type
  // double, the default type in Matlab, and a row vector
  if ((nrhs > 2) && !mxIsDouble(prhs[2])) {
    mexErrMsgTxt("List of target nodes must be of type double");
  }
  if ((nrhs > 2) && mxGetM(prhs[2])!=1) {
    mexErrMsgTxt("List of target nodes must be a row vector");
  }

  // temporal storage for the output of Dijkstra's algorithm applied
  // to any one source node
  Dsmall = (double *) mxCalloc(M , sizeof(double));
  Psmall = (double *) mxCalloc(M , sizeof(double));
  
  if (Dsmall == NULL || Psmall == NULL) {
    mexErrMsgTxt("Memory allocation failed");
  }
  
  // dealing with sparse array
  sr      = mxGetPr(prhs[0]);
  irs     = mxGetIr(prhs[0]);
  jcs     = mxGetJc(prhs[0]);

  // list to keep the array of target nodes. A list is a good data
  // structure for this, because it allows efficient removing of
  // elements, so we know which targets are still missing
  std::list<mwIndex> target;
  
  // pointer to the input array with the targets
  double *ptarget = NULL;
  if (nrhs > 2) {
    ptarget = mxGetPr(prhs[2]);
    if (ptarget == NULL) {
      mexErrMsgTxt("Cannot get pointer to the input list of targets");
    }
  }
  
  // populate the list of targets with the array of targets, if one
  // has been provided (note that even if the user doesn't pass an
  // argument for the list of targets, Matlab is going to consider
  // that an argument of size (1,1) with value 1 was passed. The only
  // way to know wether the list was given is to check the number of
  // input argument
  if ((nrhs > 2) && (prhs[2] != NULL) && !mxIsEmpty(prhs[2])) {
    target.assign(mxGetPr(prhs[2]), mxGetPr(prhs[2]) + mxGetN(prhs[2]));
  }

  // NB: at this point in the program, the list target is empty if the
  // user didn't provide a list of targets, or has some elements if
  // she did

  // setup for the Fibonacci heap
  for (i=0; i<MS; i++) {

    if ((theHeap = new FibHeap) == NULL || (A = new HeapNode[M+1]) == NULL) {
      mexErrMsgTxt("Memory allocation failed-- ABORTING.\n");
    }
    
    theHeap->ClearHeapOwnership();
    
    S = (long int) *(SS + i);
    S--;
    
    if ((S < 0) || (S > M-1)) mexErrMsgTxt("Source node out of bounds");

    /* -------------------------------------------------------------------------------------------------
       run the dijkstra code 
       ------------------------------------------------------------------------------------------------- */
    
    dodijk_sparse(M,N,S,Psmall,Dsmall,sr,irs,jcs,A,theHeap,target);
    
    for (j=0; j<M; j++) {
      // copy distance values and predecessor indices to output 
      *(D + j*MS + i) = *(Dsmall + j);
      *(P + j*MS + i) = *(Psmall + j) + 1;
      //   P[j] = (double)(Psmall[j] + 1);
      
    }
    
    /* -------------------------------------------------------------------------------------------------
       end of the dijkstra code 
       ------------------------------------------------------------------------------------------------- */
    
    delete theHeap;
    delete[] A;
  } 
    
}
