/*
 Copyright (c) 2007-2016 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 Frits Dannenberg (fdann@caltech.edu)
 */

// StrandOrdering object
// used to track sequences and structures within a complex for easy printing, etc.
#include "scomplex.h" // implicitly includes strandordering.h, and is necessary for proper ordering. TODO: decorrelate these headers, they should be independent.
#include <string>
#include <sstream>
#include <assert.h>
#include <iostream>
#include <utility.h>

using std::cout;

orderingList::orderingList(int insize, int n_id, char *inTag, char *inSeq, char *inCodeSeq, char* inStruct) {
	size = insize;
	uid = n_id;

	// ???
	thisTag = new char[strlen(inTag) + 1];
	thisSeq = new char[size + 1];
	thisCodeSeq = new char[size + 1];
	thisStruct = new char[size + 1];

	assert(thisTag != NULL);
	assert(thisSeq != NULL);
	assert(thisCodeSeq != NULL);
	assert(thisStruct != NULL);

	strncpy(thisTag, inTag, strlen(inTag) + 1);
	strncpy(thisSeq, inSeq, size);
	strncpy(thisCodeSeq, inCodeSeq, size);
	strncpy(thisStruct, inStruct, size);
	thisSeq[size] = '\0';
	thisCodeSeq[size] = '\0';
	thisStruct[size] = '\0';

	next = prev = NULL;
	thisLoop = NULL;

}

orderingList::~orderingList(void) {
	if (thisTag != NULL)
		delete[] thisTag;
	if (thisSeq != NULL)
		delete[] thisSeq;
	if (thisCodeSeq != NULL)
		delete[] thisCodeSeq;
	if (thisStruct != NULL)
		delete[] thisStruct;
	next = prev = NULL;
	if (thisLoop != NULL)
		delete thisLoop;
	thisLoop = NULL;
}

StrandOrdering::~StrandOrdering(void) {
	orderingList *temp = first;
	orderingList *temp2;
	while (temp != NULL) {
		temp2 = temp->next;
		delete temp;
		temp = temp2;
	}

	first = last = NULL;
	count = 0;
	if (seq != NULL)
		delete[] seq;
	if (struc != NULL)
		delete[] struc;
	if (strandnames != NULL)
		delete[] strandnames;
	return;
}

void StrandOrdering::cleanup(void) {
	orderingList *temp = first;
	orderingList *temp2;
	while (temp != NULL) {
		temp2 = temp->next;
		temp->thisLoop = NULL;
		temp = temp2;
	}
}

StrandOrdering::StrandOrdering(void) {

}

StrandOrdering::StrandOrdering(orderingList *beginning, orderingList *ending, int numitems) {

	first = beginning;
	last = ending;
	count = numitems;

}

// Note that in_cseq is the code sequence (ie, not printable) and in_seq is the printable version.
StrandOrdering::StrandOrdering(char *in_seq, char *in_structure, char *in_cseq) {
	char def_tag[] = "default";

	// count the number of strands, verify balanced parentheses and connectedness.
	int total_counter = 0, strand_counter = 0, strand_size = 0, sflag = 0;
	unsigned int index = 0;
	orderingList *new_elem = NULL;

	for (index = 0; index < strlen(in_cseq); index++) {
		switch (in_structure[index]) {
		case '(':
			total_counter++;
			strand_counter++;
			strand_size++;
			break;
		case ')':
			total_counter--;
			strand_counter--;
			strand_size++;
			if (strand_counter < 0)
				sflag = 1;
			break;
		case '.':
			strand_size++;
			break;
		case '+':
			count++;
			if (strand_counter == 0 && sflag == 0)
				printf("Unconnected strand in initialized complex. Strandordering.cc\n");
			new_elem = new orderingList(strand_size, -1, def_tag, &in_seq[index - strand_size], &in_cseq[index - strand_size],
					&in_structure[index - strand_size]);
			// default tag is passed and copied, so this stack alloc should be fine.

			// TODO: do i want orderinglist to be circular? does it help anything?
			if (first == NULL)
				first = last = new_elem;
			else {
				last->next = new_elem;
				new_elem->prev = last;
				last = new_elem;
			}
			new_elem = NULL;
			strand_counter = sflag = 0;
			strand_size = 0;
			while (index < strlen(in_seq) - 1 && in_seq[index + 1] == '+')
				index++;
			break;
		}
		if (total_counter < 0)
			printf("Mismatched ) at index %d in input.\n", index);
	}
	if (total_counter > 0)
		printf("Mismatched ( in input.\n");

	if (in_seq[index - 1] != '+') {
		count++;
		new_elem = new orderingList(strand_size, -1, def_tag, &in_seq[index - strand_size], &in_cseq[index - strand_size], &in_structure[index - strand_size]);
		// TODO: do i want orderinglist to be circular? does it help anything?
		if (first == NULL)
			first = last = new_elem;
		else {
			last->next = new_elem;
			new_elem->prev = last;
			last = new_elem;
		}
		new_elem = NULL;
	}

}

StrandOrdering::StrandOrdering(char *in_seq, char *in_structure, char *in_cseq, class identList *strandids) {

	class identList *traverse = strandids;

	// count the number of strands, verify balanced parentheses and connectedness.
	int total_counter = 0, strand_counter = 0, strand_size = 0, sflag = 0;
	unsigned int index = 0;
	orderingList *new_elem = NULL;

	for (index = 0; index < strlen(in_cseq); index++) {
		switch (in_structure[index]) {
		case '(':
			total_counter++;
			strand_counter++;
			strand_size++;
			break;
		case ')':
			total_counter--;
			strand_counter--;
			strand_size++;
			if (strand_counter < 0)
				sflag = 1;
			break;
		case '.':
			strand_size++;
			break;
		case '+':
			count++;
			if (strand_counter == 0 && sflag == 0)
				printf("Unconnected strand in initialized complex. Strandordering.cc\n");
			assert(traverse != NULL);
			new_elem = new orderingList(strand_size, traverse->uid, traverse->id, &in_seq[index - strand_size], &in_cseq[index - strand_size],
					&in_structure[index - strand_size]);
			traverse = traverse->next;
			// TODO: do i want orderinglist to be circular? does it help anything?
			if (first == NULL)
				first = last = new_elem;
			else {
				last->next = new_elem;
				new_elem->prev = last;
				last = new_elem;
			}
			new_elem = NULL;
			strand_counter = sflag = 0;
			strand_size = 0;

			while (index < (strlen(in_seq) - 1) && in_seq[index + 1] == '+')
				index++;
			break;
		}
		if (total_counter < 0)
			printf("Mismatched ) at index %d in input.\n", index);
	}
	if (total_counter > 0)
		printf("Mismatched ( in input.\n");

	if (in_seq[index - 1] != '+') {
		count++;
		assert(traverse != NULL);
		new_elem = new orderingList(strand_size, traverse->uid, traverse->id, &in_seq[index - strand_size], &in_cseq[index - strand_size],
				&in_structure[index - strand_size]);
		traverse = traverse->next;

		// TODO: do i want orderinglist to be circular? does it help anything?
		if (first == NULL)
			first = last = new_elem;
		else {
			last->next = new_elem;
			new_elem->prev = last;
			last = new_elem;
		}
		new_elem = NULL;
	}
	delete strandids;
}

StrandOrdering * StrandOrdering::joinOrdering(StrandOrdering *first, StrandOrdering *second) {

	first->last->next = second->first;
	second->first->prev = first->last;
	first->last = second->last;

	first->count += second->count;

	if (first->seq != NULL) {
		delete[] first->seq;
		first->seq = NULL;
	}
	if (first->struc != NULL) {
		delete[] first->struc;
		first->struc = NULL;
	}
	if (first->strandnames != NULL) {
		delete[] first->strandnames;
		first->strandnames = NULL;
	}
	second->first = NULL;
	second->last = NULL;

	first->openInfo.upToDate = false;

	return first;
}

// StrandOrdering::reorder( OpenLoop *index )

void StrandOrdering::reorder(OpenLoop *index) {
	orderingList *traverse = first, *traverse_second = NULL;
	int loop, count;

	for (traverse = first; traverse != NULL && traverse->thisLoop != index; traverse = traverse->next)
		;

	// traverse is now either NULL (bad!) or pointing atthe entry for hte open loop we wish to pivot around. Shouldcheck to makesure this isn't the original one, as no reordering is needed in that case.
	assert(traverse != NULL);

	if (traverse == first)
		return; // no reordering needed.

	count = 0;
	for (traverse_second = traverse; traverse_second != NULL; traverse_second = traverse_second->next) {
		for (loop = 0; loop < traverse_second->size; loop++) {
			if (traverse_second->thisStruct[loop] == '(')
				count++;
			if (traverse_second->thisStruct[loop] == ')') {
				if (count == 0)
					traverse_second->thisStruct[loop] = '(';
				else
					count--;
			}
		}
	}

	count = 0;
	for (traverse_second = traverse->prev; traverse_second != NULL; traverse_second = traverse_second->prev) {
		for (loop = traverse_second->size - 1; loop >= 0; loop--) {
			if (traverse_second->thisStruct[loop] == ')')
				count++;
			if (traverse_second->thisStruct[loop] == '(') {
				if (count == 0)
					traverse_second->thisStruct[loop] = ')';
				else
					count--;
			}
		}
	}

	traverse_second = traverse->prev;
	last->next = first;
	first->prev = last;
	traverse->prev = NULL;
	traverse_second->next = NULL;
	first = traverse;
	last = traverse_second;

	if (seq != NULL) {
		delete[] seq;
		seq = NULL;
	}
	if (struc != NULL) {
		delete[] struc;
		struc = NULL;
	}
	if (strandnames != NULL) {
		delete[] strandnames;
		strandnames = NULL;
	}
}

/*


 OpenLoop *StrandOrdering::checkIDList( class identlist * stoplist, int id_count )

 */

OpenLoop *StrandOrdering::checkIDList(class identList * stoplist, int id_count) {
	orderingList *traverse = first;
	class identList *id_traverse = stoplist;
	class OpenLoop* thingtoreturn = NULL;
	int num_matched = 0;
	if (id_count != count)
		return NULL;

	while (num_matched < id_count && id_traverse != NULL) {
		if (strcmp(traverse->thisTag, id_traverse->id) == 0) {
			if (num_matched == 0)
				thingtoreturn = traverse->thisLoop;
			num_matched++;
			id_traverse = id_traverse->next;
		} else {
			if (num_matched > 0)
				return NULL;
		}

		if (traverse == last) {
			if (num_matched == 0)
				return NULL;
			traverse = first;
		} else
			traverse = traverse->next;
	}

	if (num_matched == id_count)
		return thingtoreturn;
	else
		return NULL;

}

/*

 int StrandOrdering::checkIDBound( char *id )

 */

int StrandOrdering::checkIDBound(char *id) {
	orderingList *traverse = first;
	//class identlist *id_traverse = stoplist;
	//class OpenLoop *thingtoreturn;
	//int num_matched = 0;
	//if( id_count != count ) return NULL;

	unsigned int loop;
	int flag;
	while (traverse != NULL) {
		if (strcmp(traverse->thisTag, id) == 0) {
			flag = 0;
			for (loop = 0; loop < strlen(traverse->thisStruct) && (flag == 0); loop++)
				if (traverse->thisStruct[loop] == '.')
					flag = 1;

			if (flag == 0)
				return 1;
		}
		traverse = traverse->next;
	}
	return 0;
}

// void generateFlatSequence( char **sequence, char **structure, char **code_sequence
// -- Returns a flat representation of the strand ordering's sequence, structure and coded sequence. Used by SComplex::generateLoops() to re-use the old generate loops code.
// Note that the returned arrays are allocated here, but expected to be deallocated by the calling function.
// Sequence seperation is indicated by a single '_' character.
void StrandOrdering::generateFlatSequence(char **sequence, char **structure, char **code_sequence) {

	int totallength = 0, index = 0, cpos = 0;
	orderingList *traverse = first;

	openInfo.upToDate = false;

	for (index = 0; index < count; index++, traverse = traverse->next)
		totallength += traverse->size;
	totallength += count - 1;
	//  printf("Total sequence length w/breaks: %d\n",totallength);

	*sequence = new char[totallength + 1];
	*structure = new char[totallength + 1];
	*code_sequence = new char[totallength + 1];

	for (index = 0, cpos = 0, traverse = first; index < count; index++, traverse = traverse->next) {
		strncpy(&((*sequence)[cpos]), traverse->thisSeq, traverse->size);
		strncpy(&((*structure)[cpos]), traverse->thisStruct, traverse->size);
		strncpy(&((*code_sequence)[cpos]), traverse->thisCodeSeq, traverse->size);
		cpos += traverse->size;

		if (index != count - 1) {
			(*sequence)[cpos] = '_';
			(*structure)[cpos] = '.';
			(*code_sequence)[cpos] = '_';

			cpos++;
		} else {
			(*sequence)[cpos] = '\0';
			(*structure)[cpos] = '\0';
			(*code_sequence)[cpos] = '\0';
		}
	}
}

// char *convertIndex( int index)
// -- converts an index into a flat char sequence returned by generateFlatSequence into an appropriate pointer into the particular strand's code sequence.
char *StrandOrdering::convertIndex(int index) {
	int cpos, cstrand;
	orderingList *traverse;
	for (cpos = 0, cstrand = 0, traverse = first; cstrand < count; cstrand++, traverse = traverse->next) {
		if (index < cpos + traverse->size) // index is into the current strand
			return &traverse->thisCodeSeq[index - cpos];
		else if (index == cpos + traverse->size)
			; //fprintf(stderr, "Strandordering.cc, convertIndex: indexed into a strand break.\n");

		cpos += traverse->size + 1;
	}
	fprintf(stderr, "strandordering.cc, convertIndex: index out of bounds.\n");
	return NULL;
}

// Used for delete moves to get the actual Open loop and location within which is to be joined.
OpenLoop* StrandOrdering::getIndex(JoinCriteria& crit, int site, char **location, bool useArr) {

	// FD 2016 11-14: Adjusting this to work with useArr
	// Type: this is the type of base we are looking for?
	// location -- this is an OUTPUT variable.
	// location is a pointer to a char in the char** seq array

//	traverse;

	int* index = &crit.index[site];
	char type = crit.types[site];

	if (!useArr) {

		for (orderingList* traverse = first; traverse != NULL; traverse = traverse->next) {

			assert(traverse->thisLoop != NULL);

			BaseCount& baseCount = traverse->thisLoop->getFreeBases(useArr);

			if (*index < baseCount.count[type]) {

				*location = traverse->thisLoop->getBase(type, *index, useArr);

				return traverse->thisLoop;

			} else {

				*index = *index - baseCount.count[type];

			}

		}

	} else {

//		cout << crit << std::endl;

		// FD: nearly the same, but we have to call functions that consider the local context.
		// It is possible to overlap the two code paths, if we just asign the same local context
		// to each nucleotide, depending on a static global useArr trigger.

		for (orderingList* traverse = first; traverse != NULL; traverse = traverse->next) {

			assert(traverse->thisLoop != NULL);

			map<HalfContext, BaseCount>& myTally = traverse->thisLoop->getOpenInfo().tally;

//			if (!myTally.count(crit.half[site])) {
//
//				OpenInfo& openInfo = traverse->thisLoop->getOpenInfo();
//
//				cout << "printing myTally ************************* \n";
////				cout << baseCount << "\n";
//
//				cout << "site= " << site << "\n";
//				cout << "crit-half[site]= " << crit.half[site] << "\n";
//
//				cout << " *index       = " << *index << "\n";
//				cout << " OpenInfo     = " << openInfo << " " << std::endl;
//				cout << "(int) type    = " << (int) type << std::endl;
//
//			}

//			assert(openInfo.tally.count(crit.half[site]));

//			assert(myTally.count(crit.half[site]));

//			assert(traverse->thisLoop->getOpenInfo().tally.count(crit.half[site]));

//			BaseCount& baseCount = openInfo.tally.find(crit.half[site])->second;

			//BaseCount& baseCount = traverse->thisLoop->getOpenInfo().tally.find(crit.half[site])->second;

//			if (utility::debugTraces) {
//
//				cout << "printing baseCount ************************* \n";
//				cout << baseCount << "\n";
//
//				cout << "site= " << site << "\n";
//				cout << "crit-half[site]= " << crit.half[site] << "\n";
//
//				cout << " *index               = " << *index << "\n";
//				cout << "baseCount.count[type] = " << baseCount.count[type] << "\n";
//				cout << "(int) type            = " << (int) type << std::endl;
//
//			}

			if (myTally.count(crit.half[site])) { // FD: there is at least one of the correct type in this openLoop

				BaseCount& baseCount = myTally.find(crit.half[site])->second;

				if (*index < baseCount.count[type]) {

					if (utility::debugTraces) {

						cout << traverse->thisLoop->toString() << endl;

					}

					*location = traverse->thisLoop->getBase(type, *index, crit.half[site]);

//					cout << "Executing merging! ****" << std::endl;

					return traverse->thisLoop;

				} else {

					*index = *index - baseCount.count[type];

				}

			}

		}

	}

	assert(0);
	return NULL;
}

// In this case, the thisLoop data members are not initialized when using the standard constructor, and need to be associated during the scomplex's initialization (as only at that point will the open loops be created. 
// void addOpenLoop( OpenLoop *newLoop, int index)
// 

void StrandOrdering::addOpenLoop(OpenLoop *newLoop, int index) {

	openInfo.upToDate = false;

	int cpos, cstrand;
	orderingList *traverse;
	for (cpos = 0, cstrand = 0, traverse = first; cstrand < count; cstrand++, traverse = traverse->next) {
		if (index < cpos + traverse->size) // index is into the current strand
				{
			assert(traverse->thisLoop == NULL);
			traverse->thisLoop = newLoop;
			return; // loop terminates once association occurs, function complete.
		} else if (index == cpos + traverse->size)
			; //fprintf(stderr, "Strandordering.cc, convertIndex: indexed into a strand break.\n");
			  // message removed for spam purposes.

		cpos += traverse->size + 1;
	}
	fprintf(stderr, "strandordering.cc, addOpenLoop: index out of bounds, no openloop associated.\n");
	return;
}

StrandOrdering *StrandOrdering::breakOrdering(Loop *firstOldBreak, Loop *secondOldBreak, Loop *firstNewBreak, Loop *secondNewBreak) {
	orderingList *temp = NULL, *temp2 = NULL, *traverse, *extra = NULL;
	StrandOrdering *newOrdering;

	openInfo.upToDate = false;

	int numitems = 0;
	for (traverse = first; traverse != NULL; traverse = traverse->next) {
		if (traverse->thisLoop == (OpenLoop *) firstOldBreak) {
			if (temp == NULL)
				temp = traverse;
			else
				temp2 = traverse;
			traverse->thisLoop = (OpenLoop *) firstNewBreak;
		}
		if (traverse->thisLoop == (OpenLoop *) secondOldBreak) {
			if (temp == NULL)
				temp = traverse;
			else
				temp2 = traverse;
			traverse->thisLoop = (OpenLoop *) secondNewBreak;
		}
		if (temp != NULL && temp2 == NULL)
			numitems++;
	}
	if (temp == first) // one was initial open loop
			{
		extra = temp2->prev;
		extra->next = NULL;
		temp2->prev = NULL;
		newOrdering = new StrandOrdering(temp2, last, count - numitems);
		last = extra;
		count = numitems;
	} else {
		extra = temp2->prev;
		extra->next = NULL;
		temp2->prev = temp->prev;
		temp->prev->next = temp2;
		temp->prev = NULL;

		newOrdering = new StrandOrdering(temp, extra, numitems);
		count = count - numitems;
	}

	if (seq != NULL) {
		delete[] seq;
		seq = NULL;
	}
	if (struc != NULL) {
		delete[] struc;
		struc = NULL;
	}
	if (strandnames != NULL) {
		delete[] strandnames;
		strandnames = NULL;
	}
	firstOldBreak->cleanupAdjacent();
	secondOldBreak->cleanupAdjacent();
	delete firstOldBreak;
	delete secondOldBreak;

	return newOrdering;
}

char *StrandOrdering::getSequence(void) {
	if (seq != NULL)
		return seq;

	int totallength = 0, index = 0, cpos = 0;
	orderingList *traverse = first;
	for (index = 0; index < count; index++, traverse = traverse->next)
		totallength += traverse->size;
	totallength += count - 1;
	//  printf("Total sequence length w/breaks: %d\n",totallength);

	seq = new char[totallength + 1];
	struc = new char[totallength + 1];

	for (index = 0, cpos = 0, traverse = first; index < count; index++, traverse = traverse->next) {
		strncpy(&(seq[cpos]), traverse->thisSeq, traverse->size);
		strncpy(&(struc[cpos]), traverse->thisStruct, traverse->size);

		cpos += traverse->size;

		if (index != count - 1) {
			seq[cpos] = '+';
			struc[cpos] = '+';
			cpos++;
		} else {
			seq[cpos] = '\0';
			struc[cpos] = '\0';
		}
	}

	return seq;
}

char *StrandOrdering::getStructure(void) {
	if (struc != NULL)
		return struc;

	int totallength = 0, index = 0, cpos = 0;
	orderingList *traverse = first;
	for (index = 0; index < count; index++, traverse = traverse->next)
		totallength += traverse->size;
	totallength += count - 1;
	//  printf("Total sequence length w/breaks: %d\n",totallength);

	seq = new char[totallength + 1];
	struc = new char[totallength + 1];

	for (index = 0, cpos = 0, traverse = first; index < count; index++, traverse = traverse->next) {
		strncpy(&(seq[cpos]), traverse->thisSeq, traverse->size);
		strncpy(&(struc[cpos]), traverse->thisStruct, traverse->size);

		cpos += traverse->size;

		if (index != count - 1) {
			seq[cpos] = '+';
			struc[cpos] = '+';
			cpos++;
		} else {
			seq[cpos] = '\0';
			struc[cpos] = '\0';
		}
	}

	return struc;
}

char *StrandOrdering::getStrandNames(void) {
	orderingList *traverse = NULL;
	int size = 0;
	if (strandnames != NULL)
		return strandnames;

	for (traverse = first; traverse != NULL; traverse = traverse->next) {
		size += strlen(traverse->thisTag) + 1;
		size += 8;  // temp adjustment for unique ids.
	}

	size += 3; // trailing \0 plus extra space.

	strandnames = new char[size];
	strandnames[0] = '\0';

	char tmpstr[8];
	for (traverse = first; traverse != NULL; traverse = traverse->next) {
		snprintf(tmpstr, 7, "%d:", traverse->uid);
		strcat(strandnames, tmpstr);
		strcat(strandnames, traverse->thisTag);
		if (traverse->next != NULL)
			strcat(strandnames, ",");
	}

	return strandnames;
}

Loop *StrandOrdering::getLoop(void) {
	return first->thisLoop;
}

void StrandOrdering::replaceOpenLoop(Loop *oldLoop, Loop *newLoop) {

	openInfo.upToDate = false;

	orderingList *traverse = NULL;
	for (traverse = first; traverse != NULL; traverse = traverse->next) {
		if (traverse->thisLoop == (OpenLoop *) oldLoop) {
			traverse->thisLoop = (OpenLoop *) newLoop;
			return;
		}
	}

	assert(0); // no loop matched, that's bad.
}

BaseCount& StrandOrdering::getExteriorBases(bool useArr) {
	orderingList *traverse = NULL;

	exteriorBases.clear();

	for (traverse = first; traverse != NULL; traverse = traverse->next) {

		assert(traverse->thisLoop != NULL);

		BaseCount& free_bases = traverse->thisLoop->getFreeBases(useArr);

		exteriorBases.increment(free_bases);

	}

	return exteriorBases;
}

OpenInfo& StrandOrdering::getOpenInfo(void) {

	if (openInfo.upToDate) {

		return openInfo;

	} // else

	openInfo.clear();

	for (orderingList * traverse = first; traverse != NULL; traverse = traverse->next) {

		OpenInfo info = traverse->thisLoop->getOpenInfo();

		openInfo.increment(info);

	}

	openInfo.upToDate = true;

	return openInfo;

}

string StrandOrdering::toString(void) {

	std::stringstream ss;

	ss << "Ordering: ";

	orderingList *traverse = NULL;

	for (traverse = first; traverse != NULL; traverse = traverse->next) {

		assert(traverse->thisLoop != NULL);

		// now also print openInfos
		ss << traverse->thisLoop->typeInternalsToString();

	}

	// now print the combined openInfo
	ss << "Summed OpenInfo: \n";
	ss << getOpenInfo();

	return ss.str();

}

void StrandOrdering::addBasepair(char *first_bp, char *second_bp) {

	char *id[2] = { NULL, NULL };
	char *temp;
	orderingList *traverse = NULL;
	int iflag = 0;

	openInfo.upToDate = false;

	for (traverse = first; traverse != NULL; traverse = traverse->next, iflag = 0) {
		if (((first_bp - traverse->thisCodeSeq) < traverse->size) && ((first_bp - traverse->thisCodeSeq) >= 0)) {
			if (id[0] == NULL)
				id[0] = &traverse->thisStruct[first_bp - traverse->thisCodeSeq];
			else
				id[1] = &traverse->thisStruct[first_bp - traverse->thisCodeSeq];
			iflag = 1;
		}
		if (((second_bp - traverse->thisCodeSeq) < traverse->size) && ((second_bp - traverse->thisCodeSeq) >= 0)) {
			if (id[0] == NULL)
				id[0] = &traverse->thisStruct[second_bp - traverse->thisCodeSeq];
			else {
				temp = &traverse->thisStruct[second_bp - traverse->thisCodeSeq];
				if (iflag == 1 && (temp < id[0])) {
					id[1] = id[0];
					id[0] = temp;
				} else
					id[1] = temp;
			}
		}
	}
	assert(*id[0] == '.' && *id[1] == '.');
	*id[0] = '(';
	*id[1] = ')';
	if (seq != NULL) {
		delete[] seq;
		seq = NULL;
	}
	if (struc != NULL) {
		delete[] struc;
		struc = NULL;
	}

	return;
}

//
void StrandOrdering::breakBasepair(char *first_bp, char *second_bp) {

	char *id[2] = { NULL, NULL };
	char *temp = NULL;
	orderingList *traverse = NULL;
	int iflag = 0;

	openInfo.upToDate = false;

	for (traverse = first; traverse != NULL; traverse = traverse->next, iflag = 0) {

		if (((first_bp - traverse->thisCodeSeq) < traverse->size) && ((first_bp - traverse->thisCodeSeq) >= 0)) {
			if (id[0] == NULL)
				id[0] = &traverse->thisStruct[first_bp - traverse->thisCodeSeq];
			else
				id[1] = &traverse->thisStruct[first_bp - traverse->thisCodeSeq];
			iflag = 1;
		}
		if (((second_bp - traverse->thisCodeSeq) < traverse->size) && ((second_bp - traverse->thisCodeSeq) >= 0)) {
			if (id[0] == NULL)
				id[0] = &traverse->thisStruct[second_bp - traverse->thisCodeSeq];
			else {
				temp = &traverse->thisStruct[second_bp - traverse->thisCodeSeq];
				if (iflag == 1 && (temp < id[0])) {
					id[1] = id[0];
					id[0] = temp;
				} else
					id[1] = temp;
			}
		}
	}

// FD: id points to the characters in thisStruct that will change from ( and ) to . and .
	assert((*id[0] == '(' && *id[1] == ')'));
	*id[0] = '.';
	*id[1] = '.';
	if (seq != NULL) {
		delete[] seq;
		seq = NULL;
	}
	if (struc != NULL) {
		delete[] struc;
		struc = NULL;
	}

	return;
}

int StrandOrdering::getStrandCount(void) {
	return count;
}
