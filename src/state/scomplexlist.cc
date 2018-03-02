/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */

// FD  March 31, 2017. The following include fixes OSX compatibility issues.
#include <Python.h>

#include "scomplexlist.h"
#include <assert.h>
#include <math.h>

#include <vector>
#include <iostream>
#include <simoptions.h>
#include <utility.h>
#include <moveutil.h>
#include <assert.h>

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
	if (utility::debugTraces) {
		cout << "Done generating loops!" << endl;
	}

	thisComplex->generateMoves();

	if (utility::debugTraces) {
		cout << "Done generating moves!" << endl;
	}

}

void SComplexListEntry::regenerateMoves(void) {
	thisComplex->generateMoves();
}

void SComplexListEntry::fillData(EnergyModel *em) {

	energy = thisComplex->getEnergy() + (em->getVolumeEnergy() + em->getAssocEnergy()) * (thisComplex->getStrandCount() - 1);
	rate = thisComplex->getTotalFlux();

}

/*
 SComplexListEntry::printComplex
 */

string SComplexListEntry::toString(EnergyModel *em) {

	std::stringstream ss;

	// more comparable to NUPACK  - - FD replace assoEnergy with 2.44
	double printEnergy = (energy - (em->getVolumeEnergy() + em->getAssocEnergy()) * (thisComplex->getStrandCount() - 1));

	ss << "Complex         : " << id << " \n";
	ss << "seq, struc      : " << thisComplex->getSequence() << " - " << thisComplex->getStructure() << " \n";
	ss << "energy-ms       : " << energy;
	ss << "energy-nu,rate  : " << printEnergy << " - " << rate << "     (T=" << em->simOptions->energyOptions->getTemperature() << ")";
	ss << "\n";

	// also print info on the openloop datastructures
	// FD: removing this from the print because it was used for debugging mainly. (april 2017)
	// ss << thisComplex->ordering->toString();

	// also print moves components.
//	thisComplex->printAllMoves();

	return ss.str();

}

void SComplexListEntry::dumpComplexEntryToPython(ExportData& data) {

	data.id = id;
	data.names = string(thisComplex->getStrandNames());
	data.sequence = thisComplex->getSequence();
	data.structure = thisComplex->getStructure();
	data.energy = energy;
	data.enthalpy = thisComplex->getEnthalpy();

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

	eModel = energyModel;

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
	numOfComplexes++;
	idcounter++;

	return first;
}

/*
 SComplexList::initializeList
 */

void SComplexList::initializeList(void) {

	for (SComplexListEntry* temp = first; temp != NULL; temp = temp->next) {

		temp->initializeComplex();

		if (utility::debugTraces) {
			cout << "Done initializing a complex!" << endl;
		}

		temp->fillData(eModel);

	}

	if (utility::debugTraces) {
		cout << "Done initializing List!" << endl;
	}

}

void SComplexList::regenerateMoves(void) {

	for (SComplexListEntry* temp = first; temp != NULL; temp = temp->next) {

		temp->regenerateMoves();
		temp->fillData(eModel);

	}

}

/*
 SComplexList::getTotalFlux
 */

double SComplexList::getTotalFlux(void) {

	double total = 0.0;

	for (SComplexListEntry *temp = first; temp != NULL; temp = temp->next) {

		total += temp->rate;

	}

	joinRate = getJoinFlux();
	total += joinRate;

	return total;
}

BaseCount SComplexList::getExposedBases() {

	BaseCount output;

	for (SComplexListEntry* it = first; it != NULL; it = it->next) {

		BaseCount& ext_bases = it->thisComplex->getExteriorBases();
		output.increment(ext_bases);

	}

	return output;
}

OpenInfo SComplexList::getOpenInfo() {

	OpenInfo output;

	for (SComplexListEntry* it = first; it != NULL; it = it->next) {

		OpenInfo& ext_bases = it->thisComplex->ordering->getOpenInfo();
		output.increment(ext_bases);

	}

	return output;
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

double SComplexList::getJoinFlux(void) {

// We now compute the exterior nucleotide moves.
	if (numOfComplexes <= 1) {
		return 0.0;
	}

	if (eModel->useArrhenius()) {

		return getJoinFluxArr();

	}

	double output = 0.0;
	int moveCount = 0;

	BaseCount totalBases = getExposedBases();

	for (SComplexListEntry* temp = first; temp != NULL; temp = temp->next) {

		BaseCount& ext_bases = temp->thisComplex->getExteriorBases();
		totalBases.decrement(ext_bases);

		moveCount += totalBases.multiCount(ext_bases);

	}

// There are plenty of multi-complex structures with no moves.
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

// The trick is to compute all rates between the first complex,
// and the remaining complexes. Then, compute the rate between the second complex,
// and the remaining complexes (set minus first and second), and so on.

	double rate = 0.0;

	for (SComplexListEntry* temp = first; temp != NULL; temp = temp->next) {

		rate += computeArrBiRate(temp);

	}

	return rate;

}

double SComplexList::computeArrBiRate(SComplexListEntry* input) {

	double output = 0.0;

	SComplexListEntry* temp = input->next;
	// post: temp is pointing to the next entry

	StrandOrdering* orderIn = input->thisComplex->getOrdering();

	// now start computing rates with the remaining entries
	while (temp != NULL) {

		StrandOrdering* otherOrder = temp->thisComplex->getOrdering();
		output += cycleCrossRateArr(orderIn, otherOrder);

		temp = temp->next;
	}

	return output;

}

// Given two strand orderings, compute the bimolecular rate
// according to the arrhenius model.
double SComplexList::cycleCrossRateArr(StrandOrdering* input1, StrandOrdering* input2) {

	OpenInfo& info1 = input1->getOpenInfo();
	OpenInfo& info2 = input2->getOpenInfo();

	return info1.crossRate(info2, *eModel);

}

/*
 SComplexList::getEnergy( int volume_flag )
 */

double *SComplexList::getEnergy(int volume_flag) {

	SComplexListEntry *temp = first;
	double *energies = new double[numOfComplexes];
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
	return numOfComplexes;
}

double SComplexList::doBasicChoice(SimTimer& myTimer) {

	SComplexListEntry *temp, *temp2 = first;
	StrandComplex* newComplex = NULL;
	Move *tempmove;
	double arrType;

	if (utility::debugTraces) {

		cout << "Doing a basic choice " << endl;
		cout << myTimer;

	}

	if (myTimer.rchoice < joinRate) {

		return doJoinChoice(myTimer.rchoice);

	} else {

		myTimer.rchoice -= joinRate;

	}

	temp = first;
	StrandComplex *pickedComplex = NULL;

	while (temp != NULL) {
		if (myTimer.rchoice < temp->rate && pickedComplex == NULL) {
			pickedComplex = temp->thisComplex;
			temp2 = temp;
		}
		if (pickedComplex == NULL) {

			myTimer.rchoice -= temp->rate;

		}
		temp = temp->next;
	}
// POST: pickedComplex points to the complex that contains the executable move

	assert(pickedComplex != NULL);

	tempmove = pickedComplex->getChoice(&myTimer.rchoice);
	arrType = tempmove->getArrType();

	newComplex = pickedComplex->doChoice(tempmove);

	if (newComplex != NULL) {

		temp = addComplex(newComplex);
		temp->fillData(eModel);

	}

	temp2->fillData(eModel);

	// FD Oct 20, 2017.
	// If co-transcriptional mode is activated, and the time indicates a new nucleotide has been added,
	// then ask all  loops to regenerated their transitions -- taking the new time into account.
	if (myTimer.checkForNewNucleotide()) {

		eModel->numActiveNT = eModel->simOptions->initialActiveNT + myTimer.nuclAdded;

		if (eModel->numActiveNT % 25 == 0) {
			cout << "Nucleotide count = " << eModel->numActiveNT << "   time = " << myTimer.stime << " sec" << endl;
		}
		regenerateMoves();

	}

	return arrType;

}

/*
 SComplexList::doJoinChoice( double choice )
 */

double SComplexList::doJoinChoice(double choice) {

	// this function will return the arrType move;

	assert(numOfComplexes > 1);

	JoinCriteria crit;

	// before we do anything, print crit (this is for debugging!)
	if (utility::debugTraces) {
		cout << "For the current state: \n";
		cout << toString();
	}

	if (!eModel->useArrhenius()) {

		crit = cycleForJoinChoice(choice);

	} else {

		crit = cycleForJoinChoiceArr(choice);

	}

	// before we do anything, print crit (this is for debugging!)
	if (utility::debugTraces) {
//	if (true) {
		cout << "Found a criteria to join: \n";
		cout << crit.arrType;
	}

	assert(crit.complexes[0]!=NULL);
	assert(crit.complexes[1]!=NULL);

// here we actually perform the complex join, using criteria as input.

	SComplexListEntry *temp2 = NULL;
	StrandComplex *deleted;

	deleted = StrandComplex::performComplexJoin(crit, eModel->useArrhenius());
	for (SComplexListEntry* temp = first; temp != NULL; temp = temp->next) {

		if (temp->thisComplex == crit.complexes[0]) {
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
	numOfComplexes--;

	return crit.arrType;

}

JoinCriteria SComplexList::cycleForJoinChoice(double choice) {

	int int_choice = (int) floor(choice / eModel->applyPrefactors(eModel->getJoinRate(), loopMove, loopMove));
	BaseCount baseSum = getExposedBases();

	for (SComplexListEntry* temp = first; temp != NULL; temp = temp->next) {

		BaseCount& external = temp->thisComplex->getExteriorBases();
		baseSum.decrement(external);

		for (BaseType base : { baseA, baseT, baseG, baseC }) {

			int combinations = baseSum.count[base] * external.count[5 - base];

			if (int_choice < combinations) {

				// break both loops, because the right bases are identified.
				return findJoinNucleotides(base, int_choice, external, temp);

			} else {
				int_choice -= combinations;
			}

		}

	}

	assert(0);
	return JoinCriteria();

}

// FD: crit is an export variable, but the bool return signifies if a pair has been selected or not.
JoinCriteria SComplexList::findJoinNucleotides(BaseType base, int choice, BaseCount& external, SComplexListEntry* temp, HalfContext* lowerHalf) {

	JoinCriteria crit;

	int otherBase = 5 - (int) base;

	crit.complexes[0] = temp->thisComplex;
	crit.types[0] = otherBase;
	crit.types[1] = base;

	temp = temp->next;

	while (temp != NULL) {

		BaseCount externOther = temp->thisComplex->getExteriorBases(lowerHalf);

		if (choice < externOther.count[base] * external.count[otherBase]) {

			crit.complexes[1] = temp->thisComplex;
			crit.index[0] = (int) floor(choice / externOther.count[base]);
			crit.index[1] = choice - crit.index[0] * externOther.count[base];
			temp = NULL; // this exits the loop.

		} else {

			temp = temp->next;
			choice -= externOther.count[base] * external.count[otherBase];

		}
	}

	return crit;
}

JoinCriteria SComplexList::cycleForJoinChoiceArr(double choice) {

	// Like the non-arrhenius version, but this time, we have to cycle over all the
	// possible local structures.

	OpenInfo baseSum = getOpenInfo();

	for (SComplexListEntry* temp = first; temp != NULL; temp = temp->next) {

		OpenInfo& external = temp->thisComplex->ordering->getOpenInfo();

		baseSum.decrement(external);

		if (baseSum.numExposed > 0) {

			for (std::pair<HalfContext, BaseCount> con : baseSum.tally) {

				for (std::pair<HalfContext, BaseCount> ton : external.tally) {

					int combinations = con.second.multiCount(ton.second);

					if (combinations > 0) {

						MoveType left = moveutil::combineBi(con.first.left, ton.first.right);
						MoveType right = moveutil::combineBi(con.first.right, ton.first.left);

						double joinRate = eModel->applyPrefactors(eModel->getJoinRate(), left, right);

						double rate = joinRate * combinations;

						if (choice < rate) {

							// we have determined the HalfContexts for the upper and lower strand.
							int choice_int = floor(choice / joinRate);

							for (BaseType base : { baseA, baseT, baseG, baseC }) {

								int combinations = con.second.count[base] * ton.second.count[5 - base];

								if (choice_int < combinations) {

									// return the joining criteria;

									JoinCriteria crit = findJoinNucleotides(base, choice_int, ton.second, temp, &con.first);

									crit.half[0] = ton.first;
									crit.half[1] = con.first;

									crit.arrType = (double) moveutil::getPrimeCode(left, right);

									cout << "Joining two complexes, found the arrtype to be " << crit.arrType << "\n";
									cout << "Found the contexts to be " << moveutil::MoveToString[left] << "  " << moveutil::MoveToString[right] << "\n";
									cout.flush();

									return crit;

								} else {
									choice_int -= combinations;
								}

							}

						} else {

							choice = choice - rate;

						}

					}

				}

			}

		}
	}

	cout << "Failing. Remaining rate= " << choice << "\n";

	assert(0);
	return JoinCriteria();

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

	string output = "";

	for (SComplexListEntry *temp = first; temp != NULL; temp = temp->next) {

		output += temp->toString(eModel) + "\n";

	}

	return output;

}

void SComplexList::updateOpenInfo(void) {

	SComplexListEntry *temp = first;

	while (temp != NULL) {

		temp->thisComplex->getOpenInfo();
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
	if (max_complexes > numOfComplexes)
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
					if (strcmp(entry_traverse->thisComplex->getStructure().c_str(), traverse->structure) == 0) {
						// if the structures match exactly, we have a successful match.
						successflag = true;
					}
				} else if (traverse->type == STOPTYPE_DISASSOC) {
					// for DISASSOC type checking, we only need the strand id lists to match correctly.
					successflag = true;
				} else if (traverse->type == STOPTYPE_LOOSE_STRUCTURE) {
					successflag = checkLooseStructure(entry_traverse->thisComplex->getStructure().c_str(), traverse->structure, traverse->count);
					// the structure matches loosely (see definitions)
				} else if (traverse->type == STOPTYPE_PERCENT_OR_COUNT_STRUCTURE) {
					successflag = checkCountStructure(entry_traverse->thisComplex->getStructure().c_str(), traverse->structure, traverse->count);
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

bool SComplexList::checkLooseStructure(const char *our_struc, const char *stop_struc, int count) {
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
					remaining_distance--; // for the position we were
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

bool SComplexList::checkCountStructure(const char *our_struc, const char *stop_struc, int count) {

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
					remaining_distance--; // for the position we were
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
