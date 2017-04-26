/*
 Copyright (c) 2007-2008 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 Frits Dannenberg (fdann@caltech.edu)
 */

// Implementation of the StrandComplex object found in scomplex.h
#include <string.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <sstream>
#include "scomplex.h"

#include <utility.h>

using std::cout;

BaseCount emptyBaseCount;

// JS: i'd like to optimize this lookup. It really should be just a bitwise
//  or, and an array lookup, the extra function call annoys me.
extern int baseLookup(char base);

StrandComplex::StrandComplex(char *seq, char *struc) {

	char *tempseq = (char *) new char[strlen(seq) + 1];
	char *tempstruct = (char *) new char[strlen(struc) + 1];
	char * tempcseq = (char *) new char[strlen(seq) + 1];
	strcpy(tempseq, seq);
	strcpy(tempcseq, seq);
	strcpy(tempstruct, struc);
	for (int loop = 0; loop < strlen(tempcseq); loop++)
		tempcseq[loop] = baseLookup(tempcseq[loop]);

	beginLoop = NULL;
	kineticMoves = NULL;
	ordering = new StrandOrdering(tempseq, tempstruct, tempcseq);
	delete[] tempseq;
	delete[] tempstruct;
	delete[] tempcseq;

}

StrandComplex::StrandComplex(char *seq, char *struc, class identList *id_list) {

	char *tempseq = (char *) new char[strlen(seq) + 1];
	char *tempstruct = (char *) new char[strlen(struc) + 1];
	char * tempcseq = (char *) new char[strlen(seq) + 1];
	strcpy(tempseq, seq);
	strcpy(tempcseq, seq);
	strcpy(tempstruct, struc);
	for (int loop = 0; loop < strlen(tempcseq); loop++)
		tempcseq[loop] = baseLookup(tempcseq[loop]);

	beginLoop = NULL;
	kineticMoves = NULL;
	ordering = new StrandOrdering(tempseq, tempstruct, tempcseq, id_list);
	delete[] tempseq;
	delete[] tempstruct;
	delete[] tempcseq;

}

StrandComplex::StrandComplex(StrandOrdering *newOrdering) {
	ordering = newOrdering;
	beginLoop = ordering->getLoop();
	kineticMoves = NULL;

}

StrandComplex::~StrandComplex(void) {
	// we cannot delete this here now, as they could be associated with a strandordering that will live on when the complex dies.
	if (ordering != NULL)
		delete ordering;
}

typedef std::vector<Loop*> LoopVector;

void StrandComplex::cleanup(void) {
	LoopVector loops;
	loops.push_back(beginLoop);
	while (loops.size() > 0) {
		Loop* current = loops.back();
		loops.pop_back();
		for (int i = 0; i < current->getCurAdjacent(); i++)
			if (current->getAdjacent(i) != NULL)
				loops.push_back(current->getAdjacent(i));
		delete current;
	}
	beginLoop = NULL;
	ordering->cleanup();
}

/* 
 int StrandComplex::checkIDList( class identlist *stoplist, int id_count )

 checks a given ID list to see if it could be consistent with this complex's strands, and if so it reorders the current complex (if necessary) and returns 1, otherwise returns 0;

 */

int StrandComplex::checkIDList(class identList *stoplist, int id_count) {
	OpenLoop *temp;
	temp = ordering->checkIDList(stoplist, id_count);
	if (temp == NULL)
		return 0;

	ordering->reorder(temp);

	return 1;
}

int StrandComplex::checkIDBound(char *id) {
	return ordering->checkIDBound(id);
}

StrandComplex *StrandComplex::performComplexJoin(JoinCriteria crit, bool useArr) {

// FD 2016 Nov 14: Adjusting this to ignore the exterior nucleotides if useArr= TRUE;
// FD 2016 Dec 15: This comment is no longer applicable (full model now implemented)

	if (utility::debugTraces) {
		cout << "Perform Complex Join 1/3 ***************" << std::endl;
	}

	StrandComplex** complexes = crit.complexes;
	char* types = crit.types;
	int* index = crit.index;

	OpenLoop *loops[2];
	OpenLoop *new_loops[2] = { NULL, NULL };
	StrandOrdering *new_ordering = NULL;
	char *locations[2] = { NULL, NULL };

	if (utility::debugTraces) {
		cout << "Perform Complex Join 2/3 ***************" << std::endl;
	}

	// find the affected loops, and update indexes to be into those loops.
	loops[0] = complexes[0]->ordering->getIndex(crit, 0, &locations[0], useArr);
	loops[1] = complexes[1]->ordering->getIndex(crit, 1, &locations[1], useArr);

	// Strand Orderings are now ready to be joined.
	complexes[0]->ordering->reorder(loops[0]);
	complexes[1]->ordering->reorder(loops[1]);

	if (utility::debugTraces) {
		cout << "Perform Complex Join 3/3 ***************" << std::endl;
	}

	// Join the strand orderings.
	new_ordering = StrandOrdering::joinOrdering(complexes[0]->ordering, complexes[1]->ordering);

	// Join the open loops
	OpenLoop::performComplexJoin(loops, new_loops, types, index, crit.half, useArr);

	// add the base pair into the output structure.
	new_ordering->addBasepair(locations[0], locations[1]);

	// replace the old open loops with the new ones in the ordering
	new_ordering->replaceOpenLoop(loops[0], new_loops[0]);
	new_ordering->replaceOpenLoop(loops[1], new_loops[1]);

	complexes[0]->beginLoop = new_ordering->getLoop();

	complexes[0]->beginLoop->verifyLoop(NULL,  NULL);

	loops[0]->cleanupAdjacent();
	delete loops[0];
	loops[1]->cleanupAdjacent();
	delete loops[1];

	delete complexes[1]->ordering;
	complexes[1]->ordering = NULL;
	return complexes[1];
}

StrandComplex * StrandComplex::doChoice(Move * move) {
	// TODO: fix for two affected loops being deleted, must get a 'good' starting loop for the complex still.
	Loop *temp = NULL, *temp2 = NULL, *temp3 = NULL;
	char id2, id3;

	temp2 = move->affected[0];
	temp3 = move->affected[1];

	assert(temp2 != NULL);

	id2 = temp2->getType();
	if (temp3 != NULL)
		id3 = temp3->getType();
	else
		id3 = 0;

	if (id2 == 'O' && id3 == 'O') { // Break the complex.

		Loop *newLoop[2] = { NULL, NULL };
		StrandOrdering *newOrdering = NULL;

//		cout << "Going to break the complex!! 1/3 ********************** " << std::endl;

		ordering->breakBasepair(move->getAffected(0)->getLocation(move, 0), move->getAffected(1)->getLocation(move, 1));
		Loop::performComplexSplit(move, &newLoop[0], &newLoop[1]);

//		cout << "Going to break the complex!! 2/3 ********************** " << std::endl;

		// We now have open loop pointers to the two resulting open loops.
		// Now need to link up the new open loops correctly in the strand ordering
		// and then split the ordering, using one piece to build a new complex, which we then need to return to the calling function (presumably a system which can then add the new complex into the environment.

		newOrdering = ordering->breakOrdering(temp2, temp3, newLoop[0], newLoop[1]);
		beginLoop = ordering->getLoop();

		if (utility::debugTraces) {
			cout << "Going to break the complex!! 3/3 ********************** " << std::endl;
		}

		return (new StrandComplex(newOrdering)); // newComplex

	} else {
		if (move->getType() & MOVE_CREATE)	 // FD: test if we have a create-basepair move
			ordering->addBasepair(move->getAffected(0)->getLocation(move, 0), move->getAffected(0)->getLocation(move, 1));
		else if (move->getType() & MOVE_DELETE) // FD: test if we have a delete-basepair move
			ordering->breakBasepair(move->getAffected(0)->getLocation(move, 0), move->getAffected(1)->getLocation(move, 1));

		temp = move->doChoice();

		if (id2 == 'O')
			ordering->replaceOpenLoop(temp2, temp);
		if (temp3 != NULL && id3 == 'O')
			ordering->replaceOpenLoop(temp3, temp);

		if (beginLoop == temp2 || beginLoop == temp3) {
			if (temp != NULL)
				beginLoop = temp;
			else
				assert(0);
		}
		beginLoop->verifyLoop( NULL, NULL);
	}
	return NULL;
}

// used by generateLoops to handle loop traversals well.
struct intlist {
	int data;
	int seqlen;
	int pairtype;
	Loop *predec;
	struct intlist *next;
};

// ZIFNAB NEEDS CHANGE FOR SEQUENCE?STRUCTURE
int StrandComplex::generateLoops(void) {
	int *pairlist;
	char *newstruc;
	char *newseq;
	int loop, loop2, startpos, traverse, listlength, seqlen;
	int olflag = -1; // set to non -1 for internal open loops.
	int olseqlen = 0;
	int openloopcount; // used to track the offset for setting up adjacencies in the open loop.
	int depth = 0;
	struct intlist *stacklist, *stacklisttail, *templist = NULL, *templisttail =
	NULL, *temp_intlist;
	Loop *newLoop;
	char *sequence, *structure, *charsequence;

	// ZIFNAB: begin work here 8/2.
	// ZIFNAB: done, completed change to convertIndex notation.
	// ZIFNAB: more work 8/22: sequence, charsequence are used oddly, which one is actually the character sequence? do loops get the character sequence or the code sequence pointer?
	// ZIFNAB: completed: sequence is the code sequence, which has translated A/G/C/T but non translated special characters. get index should be returning into the code sequence.
	ordering->generateFlatSequence(&charsequence, &structure, &sequence);

	pairlist = (int *) new int[strlen(sequence) + 1];
	newstruc = (char *) new char[strlen(sequence) + 1];
	newseq = (char *) new char[strlen(sequence) + 1];

	stacklist = new struct intlist;
	stacklist->data = -1;
	stacklist->seqlen = 0;
	stacklist->predec = NULL;
	stacklist->next = NULL;
	stacklisttail = stacklist;

	strcpy(newstruc, structure);

	for (loop = 0; loop < strlen(sequence); loop++) {
		newseq[loop] = baseLookup(charsequence[loop]);
		pairlist[loop] = -1;
		if (structure[loop] == '(')
			depth++;
		else if (structure[loop] == ')') {
			for (loop2 = loop; loop2 >= 0; loop2--) {
				if (newstruc[loop2] == '(') {
					newstruc[loop2] = '.';
					newstruc[loop] = '.';
					pairlist[loop] = loop2;
					pairlist[loop2] = loop;
					loop2 = -1;
				}
			}
			depth--;
		}
	}

	if (pairlist[0] != -1)
		stacklist->pairtype = pairtypes[newseq[0]][newseq[pairlist[0]]];

	if (depth != 0) {
		printf("Mismatched Parens in Start Structure.");
		return -1;
	}

	/* Algorithm which the following while loop implements:

	 Queue q, gq  (implemented by templist, stacklist respectively)
	 the queues contain data in the form of side lengths, the type of base pairing for the branch (or loop), the predecessor loop.

	 While ( gq is not empty )
	 {
	 Pop a loop l off gq (implemented as a character position, startpos)
	 Traverse l and count branches and side lengths.
	 Add branches to q (without info on predecessor loop)
	 Classify l and generate the Loop structure L
	 Modify q to contain L as predecessor.
	 Modify l's predecessor with L (use addAdjacent)
	 add q to gq (maintain order)
	 }
	 */

	while (stacklist != NULL) // as long as we have unexplored base pairs on the stack, keep going.
	{
		templist = NULL;
		templisttail = NULL;
		listlength = 0;
		seqlen = 0;
		startpos = stacklist->data;

		if (startpos != -1) {
			if (pairlist[startpos] != -1) {

				if (pairlist[startpos + 1] != -1) // there was an immediate connection after the starting position. Stack or bulge, typically.
						{
					traverse = pairlist[startpos + 1] + 1; // add one otherwise we take the same link backwards when the while loops starts.
					startpos = startpos + 1;
					templist = (struct intlist *) new struct intlist;
					templisttail = templist;
					templist->data = startpos;
					templist->seqlen = 0;
					templist->predec = NULL;
					templist->next = NULL;
					templist->pairtype = pairtypes[newseq[startpos]][newseq[pairlist[startpos]]];
					// CHECK to make sure startpos+1 is the right index. FIXME 5/26
					listlength++;
				} else // we have unpaired bases after the initiating branch
				{
					traverse = startpos + 2;
					startpos++;
					seqlen++;
				}
			} else {
				traverse = startpos + 1;
				seqlen++;
			}
		} else // startpos == -1
		{
			traverse = startpos + 1;
		}
		if (startpos >= 0)
			if (sequence[startpos] == '_' || sequence[startpos] == '+') {
				//	    printf("Open Loop at olflag = %d\n",startpos);
				if (olflag != -1) // error, we shouldn't have more than one open loop specifier in a loop.
					printf("Multiple open loop specifiers in one loop!\n");

				olflag = startpos;
				seqlen--; // does not count towards sequence length
				olseqlen = seqlen;
			}

		// Current problem: last item generated will be the initial loop (the one which started this computation. Identify and eliminate addition/creation.
		// 2/11/04. START HERE - Resolved, see comment below
		while (traverse != startpos && traverse < strlen(sequence)) {
			if (sequence[traverse] == '_' || sequence[traverse] == '+') {
				//printf("Open Loop at olflag = %d\n",traverse);
				if (olflag != -1)			// error, we shouldn't have more than one open loop specifier in a loop.
					printf("Multiple open loop specifiers in one loop!\n");

				olflag = traverse;
				olseqlen = seqlen;
				seqlen--;			// does not count towards sequence length
			}
			if (pairlist[traverse] != -1) {
				if (pairlist[traverse] + 1 != startpos) // make sure this is not the initial pairing
						{
					if (templisttail == NULL) {
						templist = (struct intlist *) new struct intlist;
						templisttail = templist;
						templist->data = traverse;
						templist->seqlen = seqlen;
						templist->predec = NULL;
						templist->next = NULL;
						templist->pairtype = pairtypes[newseq[traverse]][newseq[pairlist[traverse]]];
						seqlen = 0;
					} else {
						templisttail->next = (struct intlist *) new struct intlist;
						templisttail->next->data = traverse;
						templisttail->next->seqlen = seqlen;
						templisttail->next->predec = NULL;
						templisttail->next->next = NULL;
						templisttail->next->pairtype = pairtypes[newseq[traverse]][newseq[pairlist[traverse]]];
						seqlen = 0;
						templisttail = templisttail->next;
					}
				}
				traverse = pairlist[traverse] + 1;
				listlength++;
			} else {
				traverse++;
				seqlen++;
			}
		}

		// classification of loop type time.
		// classification should end up with a pointer to the new loop, newLoop.
		if (olflag != -1)			// 'internal' open loop
				{
			//printf("listlength: %d\n",listlength);
			//printf("seqlen#0,1,olsq: %d,%d, %d\n", seqlen, templist?templist->seqlen:-1,olseqlen);
			int *OL_pairtypes;
			int *OL_sidelengths;
			char **OL_sequences;

			openloopcount = 0;
			// listlength is at least one.
			OL_pairtypes = (int *) new int[listlength];
			OL_sidelengths = (int *) new int[listlength + 1];
			OL_sequences = (char **) new char *[listlength + 1];
			// deletion for these is handled in the OpenLoop destructor.
			temp_intlist = templist;

			if (listlength == 1) {
				OL_sequences[0] = ordering->convertIndex(olflag); // CHANGED 01/06
				// removed +1 in index to hopefully fix the offset problems with open loops. This may require olflag to always be the last _ before the open loop, but that seems acceptable.
				OL_sidelengths[0] = seqlen - (olflag - stacklist->data - 1);
				OL_pairtypes[0] = stacklist->pairtype;
				OL_sequences[1] = ordering->convertIndex(stacklist->data);
				OL_sidelengths[1] = olflag - stacklist->data - 1;
				openloopcount = -1;
			} else // Algorithm follows:
				   // We need to find the circular rotation such that we always
				   // have the sequences in the correct 5'->3' ordering.
				   // 1. step through the list of adjacent helices till we find
				   //    the one immediately after the nick.
				   // 2. start adding adjacent helices, wrap around when we reach
				   //    the end of the list for the first time, until we've added
				   //    all listlength helices.
				   // 3. The final sequence and sidelength is the 5' dangle
				   //    adjacent to the nick.
			{
				temp_intlist = templist;
				for (loop = 0; loop < listlength - 1; loop++, temp_intlist = temp_intlist->next) {
					if (temp_intlist->data > olflag) // this data item is after the nick.
						break; // cause this loop to end.
					// temp_intlist will then be the first pairing after the nick.
				}

				OL_sequences[0] = ordering->convertIndex(olflag);
				if (temp_intlist == NULL)
					OL_sidelengths[0] = seqlen - olseqlen;
				else
					OL_sidelengths[0] = temp_intlist->seqlen - olseqlen;
				for (loop = 0; loop < listlength; loop++) {
					if (temp_intlist == NULL) {
						temp_intlist = templist;
						openloopcount = -openloopcount - 1;

						OL_pairtypes[loop] = stacklist->pairtype;
						if (loop != 0)
							OL_sidelengths[loop] = seqlen;
						OL_sequences[loop + 1] = ordering->convertIndex(stacklist->data);
					} else {
						if (openloopcount >= 0)
							openloopcount++;
						OL_pairtypes[loop] = temp_intlist->pairtype;
						if (loop != 0)
							OL_sidelengths[loop] = temp_intlist->seqlen;
						OL_sequences[loop + 1] = ordering->convertIndex(pairlist[temp_intlist->data]);
						temp_intlist = temp_intlist->next;
					}
				}
				OL_sidelengths[listlength] = olseqlen;
			}
			//for( loop = 0; loop <= listlength ; loop ++ )
			//printf("Seq %d: %s\n Length %d: %d\n",loop,OL_sequences[loop],loop,OL_sidelengths[loop]);
			newLoop = new OpenLoop(listlength, OL_pairtypes, OL_sidelengths, OL_sequences);
			newLoop->initAdjacency(-(openloopcount + 1));
			ordering->addOpenLoop((OpenLoop *) newLoop, olflag);
			olflag = -1;
		} else if (traverse > strlen(sequence) - 1) // Open Loop
				// Will need another classifier here. (for non initiating open loops) (CHECK: This should now be covered by the above case.)
				{
			int *OL_pairtypes;
			int *OL_sidelengths;
			char **OL_sequences;
			if (listlength != 0) {
				OL_pairtypes = (int *) new int[listlength];
				OL_sidelengths = (int *) new int[listlength + 1];
				OL_sequences = (char **) new char *[listlength + 1];
				// deletion for these is handled in the OpenLoop destructor.
				temp_intlist = templist;
				/*	      OL_pairtypes[0] = stacklist->pairtype;
				 OL_sidelengths[0] = seqlen;
				 OL_sequences[0] = &sequence[stacklist->data]; */
				/* Hmmm, this doesn't work. I'm currently commenting it out
				 ... the problem appears to be that we need this for non open
				 loops, but for open loops it doesn't make any sense. */
				// Possibly a problem here, need to make sure sequences get paired correctly with lengths. FIXME
				OL_sequences[0] = ordering->convertIndex(stacklist->data);
				OL_sidelengths[listlength] = seqlen;
				for (loop = 0; loop < listlength; loop++, temp_intlist = temp_intlist->next) {
					OL_pairtypes[loop] = temp_intlist->pairtype;
					OL_sidelengths[loop] = temp_intlist->seqlen;
					OL_sequences[loop + 1] = ordering->convertIndex(pairlist[temp_intlist->data]);
				}

				newLoop = new OpenLoop(listlength, OL_pairtypes, OL_sidelengths, OL_sequences);
				ordering->addOpenLoop((OpenLoop *) newLoop, stacklist->data);
			} else {
				OL_sidelengths = (int *) new int[listlength + 1];
				OL_sequences = (char **) new char *[listlength + 1];
				OL_sidelengths[0] = seqlen;
				OL_sequences[0] = ordering->convertIndex(-1);
				newLoop = new OpenLoop(0, NULL, OL_sidelengths, OL_sequences); // open chain
				ordering->addOpenLoop((OpenLoop *) newLoop, -1);
			}
		} else if (listlength > 2) // MultiLoop
				{
			int *ML_pairtypes;
			int *ML_sidelengths;
			char **ML_sequences;
			ML_pairtypes = (int *) new int[listlength];
			ML_sidelengths = (int *) new int[listlength];
			ML_sequences = (char **) new char *[listlength];
			// deletion for these is handled in the OpenLoop destructor.
			temp_intlist = templist;
			// Possibly a problem here, need to make sure sequences get paired correctly with lengths. FIXME
			/*
			 // old code for pairtypes and sidelengths for multiloop. uncomment the above and below to restore old sequencing.
			 ML_pairtypes[0] = stacklist->pairtype;
			 ML_sidelengths[0] = seqlen;
			 ML_sequences[1] = &sequence[stacklist->data];
			 for( loop = 1; loop < listlength; loop++, temp_intlist = temp_intlist->next)
			 {
			 ML_pairtypes[loop] = temp_intlist->pairtype;
			 ML_sidelengths[loop] = temp_intlist->seqlen;
			 if( loop == listlength - 1 )
			 ML_sequences[0] = &sequence[pairlist[temp_intlist->data]];
			 else
			 ML_sequences[loop+1] = &sequence[pairlist[temp_intlist->data]];
			 }
			 */

			// new code for pairtypes, sidelengths, seqs for multiloop, matching sequencing correctly.
			ML_pairtypes[0] = stacklist->pairtype;
			ML_sidelengths[0] = temp_intlist->seqlen;
			ML_sequences[1] = ordering->convertIndex(pairlist[temp_intlist->data]);
			for (loop = 1; loop < listlength; loop++) {
				ML_pairtypes[loop] = temp_intlist->pairtype;
				temp_intlist = temp_intlist->next;
				if (loop == listlength - 1) {
					ML_sidelengths[loop] = seqlen;
					ML_sequences[0] = ordering->convertIndex(stacklist->data);
				} else {
					ML_sidelengths[loop] = temp_intlist->seqlen;
					ML_sequences[loop + 1] = ordering->convertIndex(pairlist[temp_intlist->data]);
				}
			}
			// end new code.

			newLoop = new MultiLoop(listlength, ML_pairtypes, ML_sidelengths, ML_sequences);
		} else if (listlength == 1 && seqlen >= 3) // Hairpin Loop
				{
			newLoop = new HairpinLoop(stacklist->pairtype, seqlen, ordering->convertIndex(stacklist->data));
		} else if (listlength == 2 && (seqlen > 0 && templist->seqlen > 0)) // Interior Loop
				{
			newLoop = new InteriorLoop(stacklist->pairtype, templist->pairtype, templist->seqlen, seqlen, ordering->convertIndex(startpos - 1),
					ordering->convertIndex(pairlist[templist->data]), NULL,
					NULL);
		} else if (listlength == 2 && (seqlen == 0 && templist->seqlen == 0)) // Stack Loop
				{
			newLoop = new StackLoop(stacklist->pairtype, templist->pairtype, ordering->convertIndex(startpos - 1),
					ordering->convertIndex(pairlist[templist->data]));
		} else if (listlength == 2 && (seqlen == 0 || templist->seqlen == 0)) // Bulge Loop
				// this must be after stackloop, otherwise it may catch stackloop's conditions.
				{
			newLoop = new BulgeLoop(stacklist->pairtype, templist->pairtype, templist->seqlen, seqlen, ordering->convertIndex(startpos - 1),
					ordering->convertIndex(pairlist[templist->data]), NULL,
					NULL);
		} else if (1) // Error-generating Loop
		{
			// This should never happen.
		}

		// add correct predecessor to all adjacent loops.
		templisttail = templist;
		if (templisttail != NULL) {
			while (templisttail->next != NULL) {
				templisttail->predec = newLoop;
				templisttail = templisttail->next;
			}
			templisttail->predec = newLoop;
		}

		// classification is done, we should add the predecessor...
		if (stacklist->predec == NULL) //  we are either at the start, or an error occurred.
		{
			// we'll assume no error, and so the complex's beginning loop should be this one.
			beginLoop = newLoop;
		} else {
			newLoop->addAdjacent(stacklist->predec);
			stacklist->predec->addAdjacent(newLoop);
		}

		// add q to gq
		stacklisttail->next = templist;
		if (templisttail != NULL)
			stacklisttail = templisttail;
		templist = NULL;
		templisttail = NULL;

		// remove the read item from stacklist.
		templist = stacklist;
		stacklist = stacklist->next;
		templist->next = NULL;
		delete templist;

		// newLoop = NULL;   // uncomment this when all forks are implemented.
	}
	delete[] pairlist;
	delete[] newstruc;
	delete[] newseq;
	if (sequence != NULL)
		delete[] sequence;
	if (structure != NULL)
		delete[] structure;
	if (charsequence != NULL)
		delete[] charsequence;
}

void StrandComplex::printAllMoves(void) {

	beginLoop->printAllMoves(NULL);

}

string StrandComplex::toString() {

//	std::stringstream ss;
//
//	// doesn't do anything right now
////	printMyLoops(output, beginLoop);
//
//	// printing the ordering
//	ss << ordering->toString();

//	if (ordering != NULL) {
//
//		return ordering->toString();
//
//	} else {
//
//		return "ordering is NULL \n";
//	}

	return "";

//	return ss.str();

}

OpenInfo & StrandComplex::getOpenInfo() {

	return ordering->getOpenInfo();

}

StrandOrdering * StrandComplex::getOrdering() {

	return ordering;

}

double StrandComplex::getTotalFlux(void) {
	return totalFlux = beginLoop->returnFlux(NULL);
}

char *StrandComplex::getSequence(void) {
// ZIFNAB - use the ordering to return a valid character sequence representation for this complex.
	return ordering->getSequence();
}

char *StrandComplex::getStructure(void) {

	return ordering->getStructure();

}

char *StrandComplex::getStrandNames(void) {
	return ordering->getStrandNames();
}

BaseCount& StrandComplex::getExteriorBases( HalfContext* lowerHalf) {

	if (lowerHalf == NULL) {

		return ordering->getExteriorBases();

	} else {

		OpenInfo& info = ordering->getOpenInfo();

		if (info.tally.count(*lowerHalf)) {

			return info.tally.find(*lowerHalf)->second;

		} else {

//			cout << "Could not find the exteriorBases for " << lowerHalf << "\n";

			return emptyBaseCount;

		}

	}

}

double StrandComplex::getEnergy(void) {

	return beginLoop->returnEnergies( NULL);

}

void StrandComplex::generateMoves(void) {
	beginLoop->firstGen( NULL);
}

Move *StrandComplex::getChoice(double *rand_choice) {
	return beginLoop->getChoice(rand_choice, NULL);
}

int StrandComplex::getStrandCount(void) {
	return ordering->getStrandCount();
}
