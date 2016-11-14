/*
 Copyright (c) 2007-2008 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 */

#include "scomplexlist.h"
#include <assert.h>
#include <math.h>

#include <vector>
#include <iostream>
#include <simoptions.h>
#include <utility.h>

typedef std::vector<int> intvec;
typedef std::vector<int>::iterator intvec_it;

/*

 SComplexListEntry Constructor/Destructor

 */

SComplexListEntry::SComplexListEntry(StrandComplex *newComplex, int newid) {
	thisComplex = newComplex;
	energy = 0.0;
	rate = 0.0;
	ee_energy.dH = 0;
	ee_energy.nTdS = 0;
	next = NULL;
	id = newid;
}

SComplexListEntry::~SComplexListEntry(void) {
	delete thisComplex;
	if (next != NULL)
		delete next;
}

/*

 SComplexListEntry - InitializeComplex and FillData

 */

void SComplexListEntry::initializeComplex(void) {
	thisComplex->generateLoops();
	thisComplex->generateMoves();
}

void SComplexListEntry::regenerateMoves(void) {
	thisComplex->generateMoves();
}

void SComplexListEntry::fillData(EnergyModel *em) {
	energy = thisComplex->getEnergy() + (em->getVolumeEnergy() + em->getAssocEnergy()) * (thisComplex->getStrandCount() - 1);
	rate = thisComplex->getTotalFlux();
	// thisComplex->fillVisData( &visiblebases )
}

/*
 SComplexListEntry::printComplex
 */

string SComplexListEntry::toString(EnergyModel *em) {

	// print types are depreciated.

	std::stringstream ss;

	ss << "Complex  : " << id << " \n";
	ss << "seq      : " << thisComplex->getSequence() << " \n";
	ss << "struct   : " << thisComplex->getStructure() << " \n";
	ss << "energy   : " << energy - (em->getVolumeEnergy() + em->getAssocEnergy()) * (thisComplex->getStrandCount() - 1) << "\n";
	ss << "totalFlux: " << rate;

	ss << "\n";

	// also print info on the openloop datastructures

	ss << thisComplex->printStrandOrdering();

	return ss.str();

}

void SComplexListEntry::dumpComplexEntryToPython(int *our_id, char **names, char **sequence, char **structure, double *our_energy) {
	*our_id = id;
	*names = thisComplex->getStrandNames();
	*sequence = thisComplex->getSequence();
	*structure = thisComplex->getStructure();

//	cout << "structure = " << string(*structure) << "\n";
	*our_energy = energy;

}

///////////////////////////////////////////////////////////////////////
//                                                                   //
//                       SComplexList                                //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/* 
 SComplexList Constructor and Destructor

 */

SComplexList::SComplexList(EnergyModel *energyModel) {
	numentries = 0;
	first = NULL;
	eModel = energyModel;
	joinRate = 0.0;
	idcounter = 0;
}

SComplexList::~SComplexList(void) {
	SComplexListEntry* traverse;
	traverse = first;
	while (traverse != NULL) {
		traverse->thisComplex->cleanup();
		traverse = traverse->next;
	}
	if (first != NULL)
		delete first;
}

/* 
 SComplexList::addComplex( StrandComplex *newComplex );
 */

SComplexListEntry *SComplexList::addComplex(StrandComplex *newComplex) {
	if (first == NULL)
		first = new SComplexListEntry(newComplex, idcounter);
	else {
		SComplexListEntry *temp = new SComplexListEntry(newComplex, idcounter);
		temp->next = first;
		first = temp;
	}
	numentries++;
	idcounter++;

	return first;
}

/*
 SComplexList::initializeList
 */

void SComplexList::initializeList(void) {
	SComplexListEntry *temp = first;
	while (temp != NULL) {
		temp->initializeComplex();
		temp->fillData(eModel);
		temp = temp->next;
	}
}

void SComplexList::regenerateMoves(void) {
	SComplexListEntry *temp = first;
	while (temp != NULL) {
		temp->regenerateMoves();
		temp->fillData(eModel);
		temp = temp->next;
	}
}

/*
 SComplexList::getTotalFlux
 */

double SComplexList::getTotalFlux(void) {

	double total = 0.0;

	SComplexListEntry *temp = first;
	while (temp != NULL) {

		total += temp->rate;

		temp = temp->next;
	}

	joinRate = getJoinFlux();
	total += joinRate;

	return total;
}

/*
 double SComplexList::getJoinFlux( void )

 Computes the total flux of moves which join pairs of complexes.

 Algorithm:
 1. Sum all exterior bases in all complexes.
 2. For each complex, subtract that complex's exterior bases from the sum.
 2a. Using the new sum, compute sum.A * complex.T + sum.T * complex.A * sum.G * complex.C + sum.C * complex.G
 2b. Add this amount * rate per join move to total.
 */

double SComplexList::getJoinFlux(SimOptions* sOptions) {

	double output = 0.0;
	bool useArr = eModel->useArrhenius();
	BaseCounter total_bases;
	int moveCount = 0;

//	if (sOptions != NULL && sOptions->usingArrhenius()) {
//
//		return getJoinFluxArr();
//
//	}

// just use the regular code for now.

	SComplexListEntry *temp = first;

	if (numentries <= 1)
		return 0.0;

	while (temp != NULL) {

		BaseCounter& ext_bases = temp->thisComplex->getExteriorBases(useArr);
		total_bases.increment(ext_bases);

		temp = temp->next;
	}

	temp = first;
	while (temp != NULL) {

		BaseCounter& ext_bases = temp->thisComplex->getExteriorBases(useArr);
		total_bases.decrement(ext_bases);

		moveCount += total_bases.multiCount(ext_bases);

		temp = temp->next;
	}

	// There are plenty of multi-complex structures with no total moves.
	if (moveCount > 0) {

		output = (double) moveCount * eModel->getJoinRate();
		output = eModel->applyPrefactors(output, loopMove, loopMove);

	}

	return output;

}

// FD: We need to re-write this for the Arrhenius routine

// General strategy: A quick summation of all possible rates.
// On hit, we work back which transition we needed.
// This strategy can be faster, because interior nucleotides of openloops
// always have the same structure
double SComplexList::getJoinFluxArr(void) {

// there are 4x6x6 options, so these rates we have to tally and then sum.

	double rate = 0.0;

	SComplexListEntry* temp = first;

	// The trick is to compute all rates between the first complex,
	// and the remaining complexes. Then, compute the rate between the second complex,
	// and the remaining complexes (set minus first and second), and so on.

	while (temp != NULL) {

		StrandOrdering* order = temp->thisComplex->getOrdering();

		rate += computeArrBiRate(first, order);

		temp = temp->next;

	}

	return rate;

}

double SComplexList::computeArrBiRate(SComplexListEntry* input, StrandOrdering* order) {

	double output = 0.0;

	SComplexListEntry* temp = first;

	while (temp != NULL) {

		if (temp->thisComplex->getOrdering() == order) {
			break;
		}

		temp = temp->next;
	}

	// post: temp is pointing to the entry that matches order

	temp = temp->next;

	// now start computing rates with the remaining entries
	while (temp != NULL) {

		StrandOrdering* otherOrder = temp->thisComplex->getOrdering();
		output += cycleCrossRateArr(order, otherOrder);

		temp = temp->next;
	}

	return output;

}

// Given two strand orderings, compute the bimolecular rate
// according to the arrhenius model.
double SComplexList::cycleCrossRateArr(StrandOrdering* input1, StrandOrdering* input2) {

	double output = 1.0;

	orderinglist* temp1 = input1->first;
	orderinglist* temp2 = input2->first;

	while (temp1 != NULL) {

		OpenLoop* loop1 = temp1->thisLoop;

		while (temp2 != NULL) {

			OpenLoop* loop2 = temp2->thisLoop;

			output += computeCrossRateArr(loop1, loop2);

			temp2 = temp2->next;

		}

		temp1 = temp1->next;
	}

	// post: we've cycled all openloops in input1 over all openloops of input2

	return output;

}

// Given two openloops, compute the bimolecular rate between them\
// there are two phases:
// The internal exposed nucleotides react with a loop x loop local context.
// The external exposed nucleotides react with one of six possibilities
double SComplexList::computeCrossRateArr(OpenLoop* input1, OpenLoop* input2) {

	int adjacent1 = input1->numAdjacent;
	int adjacent2 = input2->numAdjacent;

	int output = 1.0;

	// there are no internal nucleotides
	if (input1->context.numExposedInternal == 0) {

		return output;

	} else {

		// loop, loop local context

	}

	// the only non-(loop, loop) local structures are in the first and last position.

	return output;

}

/*
 SComplexList::getEnergy( int volume_flag )
 */

double *SComplexList::getEnergy(int volume_flag) {
	SComplexListEntry *temp = first;
	double *energies = new double[numentries];
	int index = 0;
	while (temp != NULL) {
		energies[index] = temp->energy;

		if (!(volume_flag & 0x01))
			energies[index] -= (eModel->getVolumeEnergy() * (temp->thisComplex->getStrandCount() - 1));
		if (!(volume_flag & 0x02))
			energies[index] -= (eModel->getAssocEnergy() * (temp->thisComplex->getStrandCount() - 1));

		temp = temp->next;
		index = index + 1;
	}
	return energies;
}

/*
 SComplexList::printComplexList
 */

void SComplexList::printComplexList() {
	SComplexListEntry *temp = first;

	while (temp != NULL) {
		cout << temp->toString(eModel);
		temp = temp->next;
	}

}

SComplexListEntry *SComplexList::getFirst(void) {
	return first;
}

int SComplexList::getCount(void) {
	return numentries;
}
/*
 SComplexListEntry *SComplexList::doBasicChoice( double choice, double newtime )
 */

//SComplexListEntry *SComplexList::doBasicChoice(double choice, double newtime) {
int SComplexList::doBasicChoice(double choice, double newtime) {

	double rchoice = choice, moverate;
	int type;
	SComplexListEntry *temp, *temp2 = first;
	StrandComplex *pickedComplex = NULL, *newComplex = NULL;
	Move *tempmove;
	char *struc;

	if (rchoice < joinRate) {
		doJoinChoice(rchoice);

		RateEnv env = RateEnv(1.0, eModel, loopMove, loopMove);
		return env.arrType;

	} else {
		rchoice -= joinRate;
	}

	temp = first;
	while (temp != NULL) {
		if (rchoice < temp->rate && pickedComplex == NULL) {
			pickedComplex = temp->thisComplex;
			temp2 = temp;
		}
		if (pickedComplex == NULL) {
			//	  assert( rchoice != temp->rate && temp->next == NULL);
			rchoice -= temp->rate;

		}
		temp = temp->next;
	}

	assert(pickedComplex != NULL);

	tempmove = pickedComplex->getChoice(&rchoice);
	moverate = tempmove->getRate();
	type = tempmove->getType();
	newComplex = pickedComplex->doChoice(tempmove);
	if (newComplex != NULL) {
		temp = addComplex(newComplex);
		temp->fillData(eModel);
		temp2->fillData(eModel);

		return NULL;
	}

	temp2->fillData(eModel);
//	return temp2;

	return tempmove->getArrType();

}

/*
 SComplexList::doJoinChoice( double choice )
 */

void SComplexList::doJoinChoice(double choice) {

	SComplexListEntry *temp = first, *temp2 = NULL;

	BaseCounter baseSum;
	bool useArr = eModel->useArrhenius();

	StrandComplex *deleted;
	StrandComplex *picked[2] = { NULL, NULL };
	char types[2] = { 0, 0 };
	int index[2] = { 0, 0 };
	int int_choice;
	int total_move_count = 0;

	int_choice = (int) floor(choice / eModel->applyPrefactors(eModel->getJoinRate(), loopMove, loopMove));

	if (numentries <= 1)
		return;

	while (temp != NULL) {

		BaseCounter& ext_bases = temp->thisComplex->getExteriorBases(useArr);
		baseSum.increment(ext_bases);

		temp = temp->next;
	}

	temp = first;
	while (temp != NULL) {

		BaseCounter& external = temp->thisComplex->getExteriorBases(useArr);
		baseSum.decrement(external);

		if (int_choice < baseSum.A() * external.T()) {

			picked[0] = temp->thisComplex;
			types[0] = 4;
			types[1] = 1;
			temp = temp->next;

			while (temp != NULL) {

				BaseCounter& ext_bases_temp = temp->thisComplex->getExteriorBases(useArr);

				if (int_choice < ext_bases_temp.A() * external.T()) {

					picked[1] = temp->thisComplex;
					index[0] = (int) floor(int_choice / ext_bases_temp.A());
					index[1] = int_choice - index[0] * ext_bases_temp.A();
					temp = NULL;

				} else {

					temp = temp->next;
					int_choice -= ext_bases_temp.A() * external.T();

				}
			}
			continue; // We must have picked something, thus temp must be NULL and we need to exit the loop.
		} else
			int_choice -= baseSum.A() * external.T();

		if (int_choice < baseSum.T() * external.A()) {

			picked[0] = temp->thisComplex;
			types[0] = 1;
			types[1] = 4;
			temp = temp->next;

			while (temp != NULL) {

				BaseCounter& ext_bases_temp = temp->thisComplex->getExteriorBases(useArr);

				if (int_choice < ext_bases_temp.T() * external.A()) {

					picked[1] = temp->thisComplex;
					index[0] = (int) floor(int_choice / ext_bases_temp.T());
					index[1] = int_choice - index[0] * ext_bases_temp.T();
					temp = NULL;

				} else {

					temp = temp->next;
					int_choice -= ext_bases_temp.T() * external.A();

				}
			}
			continue;
		} else

			int_choice -= baseSum.T() * external.A();

		if (int_choice < baseSum.G() * external.C()) {

			picked[0] = temp->thisComplex;
			types[0] = 2;
			types[1] = 3;
			temp = temp->next;

			while (temp != NULL) {

				BaseCounter& ext_bases_temp = temp->thisComplex->getExteriorBases(useArr);

				if (int_choice < ext_bases_temp.G() * external.C()) {

					picked[1] = temp->thisComplex;
					index[0] = (int) floor(int_choice / ext_bases_temp.G());
					index[1] = int_choice - index[0] * ext_bases_temp.G();
					temp = NULL;

				} else {

					temp = temp->next;
					int_choice -= ext_bases_temp.G() * external.C();

				}
			}
			continue;
		} else

			int_choice -= baseSum.G() * external.C();

		if (int_choice < baseSum.C() * external.G()) {
			picked[0] = temp->thisComplex;
			types[0] = 3;
			types[1] = 2;
			temp = temp->next;
			while (temp != NULL) {

				BaseCounter& ext_bases_temp = temp->thisComplex->getExteriorBases(useArr);

				if (int_choice < ext_bases_temp.C() * external.G()) {

					picked[1] = temp->thisComplex;
					index[0] = (int) floor(int_choice / ext_bases_temp.C());
					index[1] = int_choice - index[0] * ext_bases_temp.C();
					temp = NULL;

				} else {

					temp = temp->next;
					int_choice -= ext_bases_temp.C() * external.G();

				}
			}
			continue;
		} else

			int_choice -= baseSum.C() * external.G();

		if (temp != NULL)
			temp = temp->next;
	}

	deleted = StrandComplex::performComplexJoin(picked, types, index, useArr);

	for (temp = first; temp != NULL; temp = temp->next) {

		if (temp->thisComplex == picked[0]) {
			temp->fillData(eModel);
		}

		if (temp->next != NULL) {

			if (temp->next->thisComplex == deleted) {
				temp2 = temp->next;
				temp->next = temp2->next;
				temp2->next = NULL;
				delete temp2;
			}

		}

	}

	if (first->thisComplex == deleted) {
		temp2 = first;
		first = first->next;
		temp2->next = NULL;
		delete temp2;

	}
	numentries--;

	return;
}

/*

 bool SComplexList::checkStopComplexList( class complex_item *stoplist )

 */
bool SComplexList::checkStopComplexList(class complexItem *stoplist) {

	if (stoplist->type == STOPTYPE_BOUND) {

		return checkStopComplexList_Bound(stoplist);

	} else {

		return checkStopComplexList_Structure_Disassoc(stoplist);
	}

}

string SComplexList::toString() {

	return first->toString( NULL);

}

void SComplexList::updateLocalContext(void) {

	SComplexListEntry *temp = first;
	while (temp != NULL) {
		temp->thisComplex->updateLocalContext();
		temp = temp->next;
	}

}

/*

 bool SComplexList::checkStopComplexList_Bound( class complex_item *stoplist )

 */
bool SComplexList::checkStopComplexList_Bound(class complexItem *stoplist) {
	class identList *id_traverse = stoplist->strand_ids;
	class SComplexListEntry *entry_traverse = first;
	int k_flag;
	if (stoplist->next != NULL) {
		fprintf(stderr, "ERROR: (scomplexlist.cc) Attempting to check for multiple complexes being bound, not currently supported.\n");
		return false;  // ERROR: can only check for a single complex/group being bound in current version.
	}

// Check for each listed strand ID, and whether it's bound.
// Note that we iterate through each complex in the list until we
// find one where it's bound, and then move on to the next.
//
	while (id_traverse != NULL) {
		entry_traverse = first;
		k_flag = 0;
		while (entry_traverse != NULL && k_flag == 0) {
			k_flag += entry_traverse->thisComplex->checkIDBound(id_traverse->id);
			entry_traverse = entry_traverse->next;
		}
		if (k_flag == 0)
			return false;
		id_traverse = id_traverse->next;
	}

	return true;
}

bool SComplexList::checkStopComplexList_Structure_Disassoc(class complexItem *stoplist) {
	class SComplexListEntry *entry_traverse = first;
	class complexItem *traverse = stoplist;
	int id_count = 0, max_complexes = 0;
	bool successflag = false;
	class identList *id_traverse = stoplist->strand_ids;

	while (traverse != NULL) {
		max_complexes++;
		traverse = traverse->next;
	}
	if (max_complexes > numentries)
		return 0; // can't match more entries than we have complexes.

	traverse = stoplist;
	while (traverse != NULL) {
// We are checking each entry in the list of stop complexes, verifying that it exists within our list of complexes. So the outer iteration is over the stop complexes, and the inner iteration is over the complexes existant in our system.
// If we reach the end of the iteration successfully, we have every stop complex represented by at least one complex within the system.
// WARNING:: this intentionally allows a single complex within the system to satisfy two (or more) stop complexes. We need to think further about such complicated stop conditions and what they might represent.
// WARNING (cont): One exception is we can't match a list of stop complexes that has more complexes than the system currently does. This is somewhat arbitrary and could be removed.

// count how many strands are in the stop complex. This gets used as a fast check later.
		id_count = 0;
		id_traverse = traverse->strand_ids;

		while (id_traverse != NULL) {
			id_count++;
			id_traverse = id_traverse->next;
		}

		entry_traverse = first;
		successflag = false;
		while (entry_traverse != NULL && successflag == 0) {
			// iterate check for current stop complex (traverse) in our list of system complexes (entry_traverse)
			if (entry_traverse->thisComplex->checkIDList(traverse->strand_ids, id_count) > 0) {
				// if the system complex being checked has the correct circular permutation of strand ids, continue with our checks, otherwise it doesn't match.
				if (traverse->type == STOPTYPE_STRUCTURE) {
					if (strcmp(entry_traverse->thisComplex->getStructure(), traverse->structure) == 0) {
						// if the structures match exactly, we have a successful match.
						successflag = true;
					}
				} else if (traverse->type == STOPTYPE_DISASSOC) {
					// for DISASSOC type checking, we only need the strand id lists to match correctly.
					successflag = true;
				} else if (traverse->type == STOPTYPE_LOOSE_STRUCTURE) {
					successflag = checkLooseStructure(entry_traverse->thisComplex->getStructure(), traverse->structure, traverse->count);
					// the structure matches loosely (see definitions)
				} else if (traverse->type == STOPTYPE_PERCENT_OR_COUNT_STRUCTURE) {
					successflag = checkCountStructure(entry_traverse->thisComplex->getStructure(), traverse->structure, traverse->count);
					// this structure matches to within a % of the correct base pairs, note that %'s are converted to raw base counts by the IO system.
				}
			}
			entry_traverse = entry_traverse->next;
		}
		if (!successflag)
			return false;
// we did not find a successful match for this stop complex in any system complex.

// otherwise, we did, try checking the next stop complex.
		traverse = traverse->next;
	}
	return true;
}

/*
 Methods used for checking loose structure definitions and counting structure defs.

 int SComplexList::checkLooseStructure( char *our_struc, char *stop_struc, int count );
 int SComplexList::checkCountStructure( char *our_struc, char *stop_struc, int count );

 */

bool SComplexList::checkLooseStructure(char *our_struc, char *stop_struc, int count) {
	int loop, len;
	intvec our_pairs, stop_pairs;
	int remaining_distance = count;

	len = strlen(our_struc);
	if (len != strlen(stop_struc))
		return false;  // something weird happened, as it should have the
// same ID list...

	for (loop = 0; loop < len; loop++) {
		if (stop_struc[loop] != '*')
			if (our_struc[loop] != stop_struc[loop]) {
				remaining_distance--;
			}

		if (our_struc[loop] == '(')
			our_pairs.push_back(loop);
		if (stop_struc[loop] == '(')
			stop_pairs.push_back(loop);

		if (our_struc[loop] == ')' && stop_struc[loop] == ')') {
			if (our_pairs.back() != stop_pairs.back()) {
				remaining_distance--; // for position loop, which had
									  // ),) but they were paired wrong
				if (our_struc[stop_pairs.back()] == '(')
					remaining_distance--;								  // for the position we were
																		  // paired with in stop_struc,
																		  // because it was ( in our_struc
																		  // as well, but paired wrong
			}
			our_pairs.pop_back();
			stop_pairs.pop_back();
		} else {
			if (our_struc[loop] == ')')
				our_pairs.pop_back();
			if (stop_struc[loop] == ')') {
				if (our_struc[stop_pairs.back()] == '(')
					remaining_distance--;  // for the position we were
										   // paired with in stop_struc,
										   // because it was ( in our
										   // struc but paired wrong. Note
										   // we have already subtracted
										   // for current position loop,
										   // as our_struc[loop] !=
										   // stop_struc[loop] in this
										   // conditional block.
				stop_pairs.pop_back();
			}
		}

		if (remaining_distance < 0)
			return false;
	}
	return true;
}

bool SComplexList::checkCountStructure(char *our_struc, char *stop_struc, int count) {

	int loop, len;
	intvec our_pairs, stop_pairs;
	int remaining_distance = count;

	len = strlen(our_struc);
	if (len != strlen(stop_struc))
		return false;  // something weird happened, as it should have the
// same ID list...

	for (loop = 0; loop < len; loop++) {
		if (our_struc[loop] != stop_struc[loop]) {
			remaining_distance--;
		}

		if (our_struc[loop] == '(')
			our_pairs.push_back(loop);
		if (stop_struc[loop] == '(')
			stop_pairs.push_back(loop);

		if (our_struc[loop] == ')' && stop_struc[loop] == ')') {
			if (our_pairs.back() != stop_pairs.back()) {
				remaining_distance--; // for position loop, which had
									  // ),) but they were paired wrong
				if (our_struc[stop_pairs.back()] == '(')
					remaining_distance--;								  // for the position we were
																		  // paired with in stop_struc,
																		  // because it was ( in our_struc
																		  // as well, but paired wrong
			}
			our_pairs.pop_back();
			stop_pairs.pop_back();
		} else {
			if (our_struc[loop] == ')')
				our_pairs.pop_back();
			if (stop_struc[loop] == ')') {
				if (our_struc[stop_pairs.back()] == '(')
					remaining_distance--;  // for the position we were
										   // paired with in stop_struc,
										   // because it was ( in our
										   // struc but paired wrong. Note
										   // we have already subtracted
										   // for current position loop,
										   // as our_struc[loop] !=
										   // stop_struc[loop] in this
										   // conditional block.
				stop_pairs.pop_back();
			}
		}

		if (remaining_distance < 0)
			return false;
	}
	return true;
}
