/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

// Implementation of the StrandComplex object found in scomplex.h
#include <string.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <sstream>
#include "scomplex.h"
#include "simtimer.h"
#include "simoptions.h"
#include "basetype.h"

#include <utility.h>

using std::cout;

BaseCount emptyBaseCount;

// JS: i'd like to optimize this lookup. It really should be just a bitwise
//  or, and an array lookup, the extra function call annoys me.

extern BaseType baseLookup(char base);

StrandComplex::StrandComplex(char *seq, char *struc)
{

    int len = static_cast<int>(strlen(seq));
	char *tempseq = (char *)new char[len + 1];
	char *tempstruct = (char *)new char[len + 1];
	BaseType *temp_bseq = (BaseType *)new BaseType[len];
	strcpy(tempseq, seq);
	strcpy(tempstruct, struc);
	for (int loop = 0; loop < len; loop++)
		temp_bseq[loop] = baseLookup(tempseq[loop]);

	beginLoop = NULL;
	ordering = new StrandOrdering(tempseq, tempstruct, temp_bseq);
	// These pointers are copied in the orderingList constructor
	// so we can delete them here.
	delete[] tempseq;
	delete[] tempstruct;
	delete[] temp_bseq;
}

StrandComplex::StrandComplex(char *seq, char *struc, class identList *id_list, bool debug)
{
    int len = static_cast<int>(strlen(seq));
	char *tempseq = (char *)new char[len + 1];
	char *tempstruct = (char *)new char[len + 1];
	BaseType *temp_bseq = (BaseType *)new BaseType[len];
	strcpy(tempseq, seq);
	strcpy(tempstruct, struc);
	for (int loop = 0; loop < len; loop++)
	{
		temp_bseq[loop] = baseLookup(tempseq[loop]);
	}

	if (debug) {
		cout << "StrandComplex:" << endl;
		printf(" seq='%s', cseq='", tempseq);
		for (int i = 0; i < len; i++)
			cout << (int) temp_bseq[i] << ",";
		printf(" struct='%s'\n", tempstruct);
		cout << flush;
	}

	beginLoop = NULL;
	ordering = new StrandOrdering(tempseq, tempstruct, temp_bseq, id_list);
	delete[] tempseq;
	delete[] tempstruct;
	delete[] temp_bseq;
}

StrandComplex::StrandComplex(StrandOrdering *newOrdering)
{
	ordering = newOrdering;
	beginLoop = ordering->getLoop();
}

StrandComplex::~StrandComplex(void)
{
	// we cannot delete this here now, as they could be associated with a strandordering that will live on when the complex dies.
	if (ordering != NULL)
		delete ordering;
}

typedef std::vector<Loop *> LoopVector;

void StrandComplex::cleanup(void)
{
	LoopVector loops;
	loops.push_back(beginLoop);
	while (loops.size() > 0)
	{
		Loop *current = loops.back();
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

int StrandComplex::checkIDList(class identList *stoplist, int id_count)
{
	OpenLoop *temp;
	temp = ordering->checkIDList(stoplist, id_count);
	if (temp == NULL)
		return 0;

	ordering->reorder(temp);

	return 1;
}

int StrandComplex::checkIDBound(char *id)
{
	return ordering->checkIDBound(id);
}

StrandComplex *StrandComplex::performComplexJoin(JoinCriteria crit, EnergyModel *energyModel)
{
	bool useArr = energyModel->useArrhenius(),
		 debug = energyModel->simOptions->debug;

	if (debug)
		cout << "Perform Complex Join 1/3 ***************" << std::endl;

	StrandComplex **complexes = crit.complexes;
	BaseType *types = crit.types;
	int *index = crit.index;

	OpenLoop *loops[2];
	OpenLoop *new_loops[2] = {NULL, NULL};
	StrandOrdering *new_ordering = NULL;
	BaseType *locations[2] = {NULL, NULL};

	if (debug)
		cout << "Perform Complex Join 2/3 ***************" << std::endl;

	// find the affected loops, and update indexes to be into those loops.
	loops[0] = complexes[0]->ordering->getIndex(crit, 0, &locations[0], useArr);
	loops[1] = complexes[1]->ordering->getIndex(crit, 1, &locations[1], useArr);

	// Strand Orderings are now ready to be joined.
	complexes[0]->ordering->reorder(loops[0]);
	complexes[1]->ordering->reorder(loops[1]);

	if (debug)
		cout << "Perform Complex Join 3/3 ***************" << std::endl;

	// Join the strand orderings.
	new_ordering = StrandOrdering::joinOrdering(complexes[0]->ordering, complexes[1]->ordering);

	// Join the open loops
	OpenLoop::performComplexJoin(
		loops, new_loops, types, index, crit.half, energyModel);

	// add the base pair into the output structure.
	new_ordering->addBasepair(locations[0], locations[1]);

	// replace the old open loops with the new ones in the ordering
	new_ordering->replaceOpenLoop(loops[0], new_loops[0]);
	new_ordering->replaceOpenLoop(loops[1], new_loops[1]);

	complexes[0]->beginLoop = new_ordering->getLoop();

	complexes[0]->beginLoop->verifyLoop(NULL, NULL);

	loops[0]->cleanupAdjacent();
	delete loops[0];
	loops[1]->cleanupAdjacent();
	delete loops[1];

	delete complexes[1]->ordering;
	complexes[1]->ordering = NULL;
	return complexes[1];
}

StrandComplex *StrandComplex::doChoice(Move *move, SimTimer &timer, EnergyModel *energyModel)
{
	// TODO: fix for two affected loops being deleted, must get a 'good' starting loop for the complex still.
	Loop *temp = NULL, *temp2 = NULL, *temp3 = NULL;
	char id2, id3;
	bool debug = energyModel->simOptions->debug;

	temp2 = move->affected[0];
	temp3 = move->affected[1];

	if (debug)
	{

		cout << "Triggering move: " << endl;
		cout << move->toString(false);
		cout << "Affected:" << endl;
		cout << temp2->toString() << endl;

		if (temp3 != NULL)
		{
			cout << temp3->toString() << endl;
		}
	}

	assert(temp2 != NULL);

	id2 = temp2->getType();
	if (temp3 != NULL)
		id3 = temp3->getType();
	else
		id3 = 0;

	if (id2 == 'O' && id3 == 'O')
	{ // Break the complex.
		// for co-transcriptional work, this is never the case.

		Loop *newLoop[2] = {NULL, NULL};
		StrandOrdering *newOrdering = NULL;

		ordering->breakBasepair(move->getAffected(0)->getLocation(move, 0), move->getAffected(1)->getLocation(move, 1));
		Loop::performComplexSplit(move, &newLoop[0], &newLoop[1], energyModel);

		// We now have open loop pointers to the two resulting open loops.
		// Now need to link up the new open loops correctly in the strand ordering
		// and then split the ordering, using one piece to build a new complex, which we then need to return to the calling function (presumably a system which can then add the new complex into the environment.

		newOrdering = ordering->breakOrdering(temp2, temp3, newLoop[0], newLoop[1]);
		beginLoop = ordering->getLoop();

		if (debug)
			cout << "Going to break the complex!! 3/3 ********************** " << std::endl;

		return (new StrandComplex(newOrdering)); // newComplex
	}
	else
	{

		if (timer.simOptions->cotranscriptional)
		{
		// 	 if cotranscriptional is active, print the distance from the origin.
		// 	 this has to be done with pointer arithmatic, since we have not refactored this (yet).

			uint64_t dist_left_bp = move->getAffected(0)->getLocation(move, 0) - ordering->first->thisBaseSeq;
			uint64_t dist_right_bp = 0;

			if (move->getType() & MOVE_CREATE){ // FD: test if we have a create-basepair move
				dist_right_bp  =  move->getAffected(0)->getLocation(move, 1) -  ordering->first->thisBaseSeq;
			} else if (move->getType() & MOVE_DELETE) { // FD: test if we have a delete-basepair move
				dist_right_bp = move->getAffected(1)->getLocation(move, 1) -  ordering->first->thisBaseSeq;
			}

			// exit if one of the transitions are touching a frozen base
			// right now, frozen if more than 100 nt from the most active one.
			uint64_t actives= timer.nuclAdded + timer.simOptions->initialActiveNT;

			constexpr uint64_t bound = 200;
			if( dist_left_bp + bound < actives || dist_right_bp + bound < actives ){
				return NULL;
			}
		}

		if (move->getType() & MOVE_CREATE) // FD: test if we have a create-basepair move
			ordering->addBasepair(move->getAffected(0)->getLocation(move, 0), move->getAffected(0)->getLocation(move, 1));
		else if (move->getType() & MOVE_DELETE) // FD: test if we have a delete-basepair move
			ordering->breakBasepair(move->getAffected(0)->getLocation(move, 0), move->getAffected(1)->getLocation(move, 1));

		temp = move->doChoice(energyModel);

		if (id2 == 'O')
			ordering->replaceOpenLoop(temp2, temp);
		if (temp3 != NULL && id3 == 'O')
			ordering->replaceOpenLoop(temp3, temp);

		if (beginLoop == temp2 || beginLoop == temp3)
		{
			if (temp != NULL)
				beginLoop = temp;
			else
				assert(0);
		}
		beginLoop->verifyLoop(NULL, NULL);
	}
	return NULL;
}

// used by generateLoops to handle loop traversals well.
struct generateLoopsData {
	int data;
	int seqlen;
	Loop *predec = NULL;
	struct generateLoopsData *next = NULL;
};

// ZIFNAB NEEDS CHANGE FOR SEQUENCE?STRUCTURE
int StrandComplex::generateLoops(bool debug)
{

	int startpos, traverse, listlength, seqlen;
	int olflag = -1; // set to non -1 for internal open loops.
	int olseqlen = 0;
	int openloopcount; // used to track the offset for setting up adjacencies in the open loop.
	int depth = 0;
	generateLoopsData *stacklist, *stacklisttail, *temp_intlist;
	generateLoopsData *templist = NULL, *templisttail = NULL;
	Loop *newLoop;
	char *structure, *charsequence;
	BaseType *sequence;

	if (debug)
		cout << "Generating flat sequences." << endl << flush;

	// ZIFNAB: begin work here 8/2.
	// ZIFNAB: done, completed change to convertIndex notation.
	// ZIFNAB: more work 8/22: sequence, charsequence are used oddly, which one is actually the character sequence? do loops get the character sequence or the code sequence pointer?
	// ZIFNAB: completed: sequence is the code sequence, which has translated A/G/C/T but non translated special characters. get index should be returning into the code sequence.
	ordering->generateFlatSequence(&charsequence, &structure, &sequence, debug);

    int len = static_cast<int>(strlen(charsequence));
	if (debug) {
		printf(" Generating loops for: cseq='%s', seq='", charsequence);
		for (int loop = 0; loop < len; loop++)
			cout << (int) sequence[loop] << ",";
		printf(" struct='%s'\n", structure);
		cout << " strlen(charsequence) = " << strlen(charsequence) << endl << flush;
	}

	int *pairlist = (int *) new int[len + 1];
	char *newstruc = (char *) new char[len + 1];
	BaseType *newseq = (BaseType *) new BaseType[len];

	stacklist = new generateLoopsData;
	stacklist->data = -1;
	stacklist->seqlen = 0;
	stacklist->predec = NULL;
	stacklist->next = NULL;
	stacklisttail = stacklist;

	strcpy(newstruc, structure);

	for (int loop = 0; loop < len; loop++) {
		newseq[loop] = baseLookup(charsequence[loop]);
		pairlist[loop] = -1;
		if (structure[loop] == '(')
			depth++;
		else if (structure[loop] == ')') {
			for (int loop2 = loop; loop2 >= 0; loop2--) {
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

	if (depth != 0) {
		printf("Mismatched Parens in Start Structure.");
		return -1;
	}

	if (debug) {
		cout << " pairlist='";
		for (int loop = 0; loop < len; loop++)
			cout << pairlist[loop] << ",";
		cout << "'" << endl << " newstruc='";
		for (int loop = 0; loop < len; loop++)
			cout << newstruc[loop];
		cout << "'" << endl << " newseq=";
		for (int loop = 0; loop < len; loop++)
			cout << (int) newseq[loop] << ",";
		cout << "'" << endl << flush;
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

				if (pairlist[startpos + 1] != -1) { // there was an immediate connection after the starting position. Stack or bulge, typically.

					traverse = pairlist[startpos + 1] + 1; // add one otherwise we take the same link backwards when the while loops starts.
					startpos = startpos + 1;
					templist = new struct generateLoopsData;
					templisttail = templist;
					templist->data = startpos;
					templist->seqlen = 0;
					templist->predec = NULL;
					templist->next = NULL;
					// CHECK to make sure startpos+1 is the right index. FIXME 5/26
					listlength++;

				} else { // we have unpaired bases after the initiating branch

					traverse = startpos + 2;
					startpos++;
					seqlen++;
				}

			} else {

				traverse = startpos + 1;
				seqlen++;
			}

		} else { // startpos == -1

			traverse = startpos + 1;

		}

		if (startpos >= 0)
			if (sequence[startpos] == baseInvalid || charsequence[startpos] == '+')
			{
				//	    printf("Open Loop at olflag = %d\n",startpos);
				if (olflag != -1) // error, we shouldn't have more than one open loop specifier in a loop.
					printf("Multiple open loop specifiers in one loop!\n");

				olflag = startpos;
				seqlen--; // does not count towards sequence length
				olseqlen = seqlen;
			}

		// Current problem: last item generated will be the initial loop (the one which started this computation. Identify and eliminate addition/creation.
		// 2/11/04. START HERE - Resolved, see comment below
		while (traverse != startpos && traverse < len)
		{
			if (sequence[traverse] == baseInvalid || charsequence[traverse] == '+')
			{
				//printf("Open Loop at olflag = %d\n",traverse);
				if (olflag != -1) // error, we shouldn't have more than one open loop specifier in a loop.
					printf("Multiple open loop specifiers in one loop!\n");

				olflag = traverse;
				olseqlen = seqlen;
				seqlen--; // does not count towards sequence length
			}
			if (pairlist[traverse] != -1)
			{
				if (pairlist[traverse] + 1 != startpos) // make sure this is not the initial pairing
                {
					if (templisttail == NULL) {
						templist = new struct generateLoopsData;
						templisttail = templist;
						templist->data = traverse;
						templist->seqlen = seqlen;
						seqlen = 0;
					} else {
						templisttail->next = new struct generateLoopsData;
						templisttail->next->data = traverse;
						templisttail->next->seqlen = seqlen;
						seqlen = 0;
						templisttail = templisttail->next;
					}
				}
				traverse = pairlist[traverse] + 1;
				listlength++;
			}
			else
			{
				traverse++;
				seqlen++;
			}
		}

		// JS: classification of loop type time.
		// classification should end up with a pointer to the new loop, newLoop.

		if (olflag != -1) {			// 'internal' open loop

			int *OL_sidelengths;
			BaseType **OL_sequences;

			openloopcount = 0;
			// listlength is at least one.
			OL_sidelengths = (int *)new int[listlength + 1];
			OL_sequences = (BaseType **)new BaseType *[listlength + 1];

			// deletion for these is handled in the OpenLoop destructor.
			temp_intlist = templist;

			if (listlength == 1)
			{
				OL_sequences[0] = ordering->convertIndex(olflag); // CHANGED 01/06
				// removed +1 in index to hopefully fix the offset problems with open loops. This may require olflag to always be the last _ before the open loop, but that seems acceptable.
				OL_sidelengths[0] = seqlen - (olflag - stacklist->data - 1);
				OL_sequences[1] = ordering->convertIndex(stacklist->data);
				OL_sidelengths[1] = olflag - stacklist->data - 1;
				openloopcount = -1;
			}
			else // JS: Algorithm follows:
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
				for (int loop = 0; loop < listlength - 1; loop++, temp_intlist = temp_intlist->next) {
					if (temp_intlist->data > olflag) { // this data item is after the nick.
						break; // cause this loop to end.
					}
					// temp_intlist will then be the first pairing after the nick.
				}

				OL_sequences[0] = ordering->convertIndex(olflag); // label a
				if (temp_intlist == NULL)
					OL_sidelengths[0] = seqlen - olseqlen;
				else
					OL_sidelengths[0] = temp_intlist->seqlen - olseqlen;
				for (int loop = 0; loop < listlength; loop++) {
					if (temp_intlist == NULL) {
						temp_intlist = templist;
						openloopcount = -openloopcount - 1;

						if (loop != 0)
						{
							OL_sidelengths[loop] = seqlen;
						}
						OL_sequences[loop + 1] = ordering->convertIndex(stacklist->data);
					}
					else
					{
						if (openloopcount >= 0)
						{
							openloopcount++;
						}
						if (loop != 0)
						{
							OL_sidelengths[loop] = temp_intlist->seqlen;
						}
						OL_sequences[loop + 1] = ordering->convertIndex(pairlist[temp_intlist->data]);
						temp_intlist = temp_intlist->next;
					}
				}
				OL_sidelengths[listlength] = olseqlen;
			}

			newLoop = new OpenLoop(listlength, OL_sidelengths, OL_sequences);
			newLoop->initAdjacency(-(openloopcount + 1));
			ordering->addOpenLoop((OpenLoop *)newLoop, olflag);
			olflag = -1;
		}
		else if (traverse > len - 1) // Open Loop
		// Will need another classifier here. (for non initiating open loops) (CHECK: This should now be covered by the above case.)
		{
			int *OL_sidelengths;
			BaseType **OL_sequences;
			if (listlength != 0)
			{

				OL_sidelengths = (int *)new int[listlength + 1];
				OL_sequences = (BaseType **)new BaseType *[listlength + 1];
				// deletion for these is handled in the OpenLoop destructor.
				temp_intlist = templist;
				/*	      OL_pairtypes[0] = stacklist->pairtype;
				 OL_sidelengths[0] = seqlen;
				 OL_sequences[0] = &sequence[stacklist->data]; */
				/* JS: Hmmm, this doesn't work. I'm currently commenting it out
				 ... the problem appears to be that we need this for non open
				 loops, but for open loops it doesn't make any sense. */
				// Possibly a problem here, need to make sure sequences get paired correctly with lengths. FIXME

				// Jake: Starting point on an open loop is -1. We don't actually read the first value of the
				// sequence but it's still so dangerous to be doing this. I added a "wrapper" array so that if that
				// -1 is ever indexed it won't cause memory issues.

				OL_sequences[0] = ordering->convertIndex(stacklist->data);
				OL_sidelengths[listlength] = seqlen;
				for (int loop = 0; loop < listlength; loop++, temp_intlist = temp_intlist->next) {
					OL_sidelengths[loop] = temp_intlist->seqlen;
					OL_sequences[loop + 1] = ordering->convertIndex(pairlist[temp_intlist->data]);
				}
				
				newLoop = new OpenLoop(listlength, OL_sidelengths, OL_sequences);
				//exit(1);
				ordering->addOpenLoop((OpenLoop *)newLoop, stacklist->data);
			}
			else
			{

				OL_sidelengths = (int *)new int[listlength + 1];
				OL_sequences = (BaseType **)new BaseType *[listlength + 1];
				OL_sidelengths[0] = seqlen;
				OL_sequences[0] = ordering->convertIndex(-1); // this feels very dangerous (need to think about it more)
				newLoop = new OpenLoop(0, OL_sidelengths, OL_sequences); // open chain
				ordering->addOpenLoop((OpenLoop *)newLoop, -1);
			}
		}
		else if (listlength > 2) // MultiLoop
		{
			int *ML_sidelengths;
			BaseType **ML_sequences;
			ML_sidelengths = (int *)new int[listlength];
			ML_sequences = (BaseType **)new BaseType *[listlength];
			// deletion for these is handled in the OpenLoop destructor.
			temp_intlist = templist;
			// JS: Possibly a problem here, need to make sure sequences get paired correctly with lengths. FIXME

			// JS: new code for pairtypes, sidelengths, seqs for multiloop, matching sequencing correctly.
			ML_sidelengths[0] = temp_intlist->seqlen;
			ML_sequences[1] = ordering->convertIndex(pairlist[temp_intlist->data]);
			for (int loop = 1; loop < listlength; loop++) {
				temp_intlist = temp_intlist->next;
				if (loop == listlength - 1)
				{
					ML_sidelengths[loop] = seqlen;
					ML_sequences[0] = ordering->convertIndex(stacklist->data);
				}
				else
				{
					ML_sidelengths[loop] = temp_intlist->seqlen;
					ML_sequences[loop + 1] = ordering->convertIndex(pairlist[temp_intlist->data]);
				}
			}
			// end new code.

			newLoop = new MultiLoop(listlength, ML_sidelengths, ML_sequences);
		}
		else if (listlength == 1 && seqlen >= 3) // Hairpin Loop
		{
			newLoop = new HairpinLoop(seqlen, ordering->convertIndex(stacklist->data));
		}
		else if (listlength == 2 && (seqlen > 0 && templist->seqlen > 0)) // Interior Loop
		{
			newLoop = new InteriorLoop(templist->seqlen, seqlen, ordering->convertIndex(startpos - 1), ordering->convertIndex(pairlist[templist->data]), NULL,
									   NULL);
		}
		else if (listlength == 2 && (seqlen == 0 && templist->seqlen == 0)) // Stack Loop
		{
			newLoop = new StackLoop(ordering->convertIndex(startpos - 1), ordering->convertIndex(pairlist[templist->data]));
		}
		else if (listlength == 2 && (seqlen == 0 || templist->seqlen == 0)) // Bulge Loop
		// this must be after stackloop, otherwise it may catch stackloop's conditions.
		{
			newLoop = new BulgeLoop(templist->seqlen, seqlen, ordering->convertIndex(startpos - 1), ordering->convertIndex(pairlist[templist->data]), NULL,
                                    NULL);
		} else {
			// This should never happen.
			printf("StrandComplex::generateLoops() - Invalid loop!");
			throw std::bad_alloc();
		}

		// add correct predecessor to all adjacent loops.
		templisttail = templist;
		if (templisttail != NULL)
		{
			while (templisttail->next != NULL)
			{
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
		}
		else
		{
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

		newLoop = NULL;   // uncomment this when all forks are implemented.
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

	return 0;
}

void StrandComplex::printAllMoves(bool useArrhenius)
{

	beginLoop->printAllMoves(NULL, useArrhenius);
}

string StrandComplex::toString()
{

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

OpenInfo &StrandComplex::getOpenInfo()
{

	return ordering->getOpenInfo();
}

StrandOrdering *StrandComplex::getOrdering()
{

	return ordering;
}

double StrandComplex::getTotalFlux(void)
{
	return beginLoop->returnFlux(NULL);
}

int StrandComplex::getMoveCount(void)
{
	return beginLoop->getMoveCount(NULL);
}

string &StrandComplex::getSequence(void)
{
	// ZIFNAB - use the ordering to return a valid character sequence representation for this complex.
	return ordering->getSequence();
}

string &StrandComplex::getStructure(void)
{

	return ordering->getStructure();
}

char *StrandComplex::getStrandNames(void) {
	return ordering->getStrandNames();
}

BaseCount &StrandComplex::getExteriorBases(HalfContext *lowerHalf)
{

	if (lowerHalf == NULL)
	{

		return ordering->getExteriorBases();
	}
	else
	{

		OpenInfo &info = ordering->getOpenInfo();

		if (info.tally.count(*lowerHalf))
		{

			return info.tally.find(*lowerHalf)->second;
		}
		else
		{

			return emptyBaseCount;
		}
	}
}

double StrandComplex::getEnergy(EnergyModel *energyModel)
{
	return beginLoop->returnEnergies(NULL, energyModel);
}

double StrandComplex::getEnthalpy(EnergyModel *energyModel)
{
	return beginLoop->returnEnthalpies(NULL, energyModel);
}

void StrandComplex::generateMoves(EnergyModel *energyModel) {
	beginLoop->firstGen(NULL, energyModel);
}

Move *StrandComplex::getChoice(SimTimer &timer)
{
	return beginLoop->getChoice(timer, NULL);
}

int StrandComplex::getStrandCount(void)
{
	return ordering->getStrandCount();
}
