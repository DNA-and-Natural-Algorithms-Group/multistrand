/*
 Copyright (c) 2007-2010 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 */

#ifndef __STRANDORDERING_H__
#define __STRANDORDERING_H__

#include "loop.h"
#include "scomplex.h"
#include "optionlists.h"
#include <string>

// needed for the openloop components of a strand ordering

class OpenLoop;

class orderinglist {
public:
	orderinglist(int insize, int in_id, char *inTag, char *inSeq,
			char *inCodeSeq, char* inStruct);
	~orderinglist(void);
	orderinglist *next, *prev;
	char *thisTag, *thisSeq, *thisCodeSeq, *thisStruct;
	OpenLoop *thisLoop; // corresponds to the OpenLoop to the 'left' of this strand
	int size;
	int uid;
};

class StrandOrdering {
public:
	StrandOrdering(void);
	StrandOrdering(char *in_seq, char *in_struc, char *in_cseq);
	StrandOrdering(orderinglist *beginning, orderinglist *ending, int numitems);
	StrandOrdering(char *in_seq, char *in_structure, char *in_cseq,
			class identList *strandids);
	~StrandOrdering(void);
	void cleanup(void);
	static StrandOrdering * joinOrdering(StrandOrdering *first,
			StrandOrdering *second);
	StrandOrdering *breakOrdering(Loop *firstOldBreak, Loop *secondOldBreak,
			Loop *firstNewBreak, Loop *secondNewBreak); // maybe id or openloop pointer
	void reorder(OpenLoop *index); // reorder so that open loop passed is the available openloop
	void addBasepair(char *first_bp, char *second_bp);
	void breakBasepair(char *first_bp, char *second_bp);

	OpenLoop *checkIDList(class identList *stoplist, int count);
	int checkIDBound(char *id);

	// following three functions are used by SComplex::generateLoops
	// to generate the loop structure of a given complex, using a flat representation of the starting sequence and structure.
	// Note that the first function, generateFlatSequence, is given pointers to appopriate char * markers to hold the flat representation. Structure is very difficult to have based on the ordering, though, so perhaps it needs to be handled differently.
	void generateFlatSequence(char **sequence, char **structure,
			char **code_sequence);

	// this function converts an index into a previously given sequence from generateFlatSequence into a char * pointer into the appropriate strand's sequence at the given location.
	char *convertIndex(int index);

	// addOpenLoop links up the appropriate strand with the open loop involving the nick immediately before that strand in the ordering.
	void addOpenLoop(OpenLoop *newLoop, int index);

	// end functions for SComplex::generateLoops()

	// getIndex finds which open loop is associated with a particular index into
	// one exterior base type for the complex, and returns that Openloop, as well
	// as the updated index into that open loop only, and a char * pointing at
	// the particular base in the open loop.
	OpenLoop *getIndex(char type, int *index, char **location);

	// i/o routines and accessors for strandcomplex
	char *getSequence(void);
	char *getStructure(void);
	char *getStrandNames(void);
	Loop *getLoop(void);
	int getStrandCount(void);

	// updates and returns the current exterior base count.
	exterior_bases *getExteriorBases(void);
	void updateHalfContexts();
	string toString(void);

	// replaces the first open loop in the ordering with the second.
	void replaceOpenLoop(Loop *oldLoop, Loop *newLoop);

private:
	char *seq, *struc, *strandnames;
	orderinglist *first, *last;
	int count;
	exterior_bases total_exterior_bases;

};

#endif
