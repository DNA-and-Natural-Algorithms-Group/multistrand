/*
 Copyright (c) 2007-2010 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 */

#include "energymodel.h"
#include "simoptions.h"
#include "loop.h"
#include "moveutil.h"
#include "sequtil.h"

#include <iostream>
#include <fstream>

bool printedRates = false; // to print the constants file once

const double INIT_PENALTY = 2.0; //kcal / mol

EnergyModel::EnergyModel(PyObject *options) {
	// nothing yet

}

EnergyModel::EnergyModel(void) {
	// nothing yet either
}

EnergyModel::~EnergyModel(void) {
	// nothing
}

// Return state of the energyOption toggle
bool EnergyModel::useArrhenius(void) {

	return simOptions->energyOptions->usingArrhenius();

}

double expRate(double A, double E, double temperature) {

	return exp(A - E / (gasConstant * temperature));

}

void EnergyModel::setArrheniusRate(double ratesArray[], EnergyOptions* options, double temperature, int left, int right) {

	double ALeft = options->AValues[left];
	double ELeft = options->EValues[left];

	double ARight = options->AValues[right];
	double ERight = options->EValues[right];

	double kLeft = expRate(ALeft, ELeft, temperature);
	double kRight = expRate(ARight, ERight, temperature);

	ratesArray[left * MOVETYPE_SIZE + right] = kLeft * kRight;

}

void EnergyModel::computeArrheniusRates(double temperature) {

	for (int i = 0; i < MOVETYPE_SIZE; i++) {

		for (int j = 0; j < MOVETYPE_SIZE; j++) {

			arrheniusRates[i * MOVETYPE_SIZE + j] = NULL;

			setArrheniusRate(arrheniusRates, simOptions->energyOptions, temperature, i, j);

		}

	}

}

void EnergyModel::writeConstantsToFile() {

	// Print constants to file.

	std::stringstream ss;

	ss << "Sodium      :  " << simOptions->energyOptions->sodium << " M \n";
	ss << "Magnesium   :  " << simOptions->energyOptions->magnesium << " M \n";
	ss << "Temperature :  " << simOptions->energyOptions->getTemperature() << endl;
	ss << "useArr      :  " << simOptions->energyOptions->usingArrhenius() << endl;

	if (!simOptions->energyOptions->usingArrhenius()) {

		ss << "    biScale     kUni    \n";
		ss << "     " << simOptions->energyOptions->getBiScale();
		ss << "     " << simOptions->energyOptions->getUniScale();

		ss << "\n";

	} else {

		ss << "type      ";

		for (int i = 0; i < MOVETYPE_SIZE; i++) {

			ss << moveutil::MoveToString[i];
			ss << moveutil::MoveToString2[i] << " ";

		}

		ss << setprecision(8);

		ss << "\nA         ";

		for (int i = 0; i < MOVETYPE_SIZE; i++) {

			ss << simOptions->energyOptions->AValues[i] << "      ";

		}

		ss << "\nE         ";

		for (int i = 0; i < MOVETYPE_SIZE; i++) {

			ss << simOptions->energyOptions->EValues[i] << "      ";

		}

		ss << "\nR         ";

		for (int i = 0; i < MOVETYPE_SIZE; i++) {

			ss << arrheniusRates[MOVETYPE_SIZE * i + i] << "  ";

		}

		ss << " \n \n";
		ss << "    dS_A     dH_A        biScale     kUni    \n";
		ss << "    " << simOptions->energyOptions->dSA;
		ss << "     " << simOptions->energyOptions->dHA;
		ss << "     " << simOptions->energyOptions->getBiScale();
		ss << "     " << simOptions->energyOptions->getUniScale();

		ss << "\n";

		// now dumping the rate matrix too,

		for (int i = 0; i < MOVETYPE_SIZE; i++) {

			for (int j = 0; j < MOVETYPE_SIZE; j++) {

				ss << getJoinRate() * arrheniusRates[MOVETYPE_SIZE * i + j] << "  ";

			}

			ss << " \n";

		}

	}


	if (!printedRates) {

		ofstream myfile;
		myfile.open("multistrandRun.log");

		myfile << ss.str();

		myfile.close();

		cout << "Wrote constants to multistrandRun.log" << endl;

		printedRates = true;

	}

}

double EnergyModel::applyPrefactors(double tempRate, MoveType left, MoveType right) {

	if (useArrhenius()) {

		return tempRate * arrheniusRates[left * MOVETYPE_SIZE + right];

	} else {

		return tempRate;
	}

}

// FD: A base pair is present between a stacking loop and a multi loop.
// FD: We query the local context of the middle pair;
// FD: this can be either a loop, stack+loop, or stack+stack situation.
MoveType EnergyModel::getPrefactorsMulti(int index, int numAdjacent, int sideLengths[]) {

// FD: if adjacent are not neighbored in this order, then we need to replace this code with the correct mapping.
	int rightStrand = (index + 1) % numAdjacent;

	return this->prefactorInternal(sideLengths[index], sideLengths[rightStrand]);

}

MoveType EnergyModel::prefactorInternal(int sideLength1, int sideLength2) {
// FD: A base pair is present between a stacking loop and a multi loop.
// FD: We query the local context of the middle pair;
// FD: this can be either a loop, stack+loop, or stack+stack situation.

	if (sideLength1 > 0 && sideLength2 > 0) {	 // at least one unpaired base on each side.

		return loopMove;

	} else if (sideLength1 > 0 || sideLength2 > 0) {  // at least one unpaired base on at least one side

		return stackLoopMove;

	} else { // two nucleotides on each side;

		return stackStackMove;

	}

}

MoveType prefactorEndBothOpen(int left, int right) {

	if (left > 0 && right > 0) {

		return loopMove;

	} else if (left > 0 || right > 0) {

		return loopEndMove;

	} else {

		return endMove;

	}

}

MoveType prefactorEndSingleOpen(int openSide, int closedSide) {

	bool dangleOpen = (openSide > 0);
	bool dangleClosed = (closedSide > 0);

	if (dangleOpen && dangleClosed) {

		return loopMove;

	} else if (!dangleOpen && !dangleClosed) {

		return stackEndMove;

	} else if (dangleOpen) {

		return stackLoopMove;

	} else {

		return loopEndMove;

	}

}

MoveType EnergyModel::prefactorOpen(int index, int numOfSides, int sideLengths[]) {

	assert(index < (numOfSides + 1));

	bool openOnLeft = (index == 0);
	bool openOnRight = ((index + 1) == numOfSides - 1);

// each strand can either be non-existing, or single stranded, or double stranded.

	if (!openOnLeft && !openOnRight) { // this is the multi-loop case

		return prefactorInternal(sideLengths[index], sideLengths[index + 1]);

	} else if (openOnLeft && openOnRight) { // this is the "end" case

		return prefactorEndBothOpen(sideLengths[index], sideLengths[index + 1]);

	} else {
		// there is exactly one side exposed to the open
		// this side is always single stranded, or non-existing.

		if (openOnLeft) {

			return prefactorEndSingleOpen(sideLengths[index], sideLengths[index + 1]);

		} else {

			return prefactorEndSingleOpen(sideLengths[index + 1], sideLengths[index]);

		}

	}

}

double EnergyModel::singleStrandedStacking(char* sequence, int length) {

	if (simOptions->energyOptions->usingArrhenius() && length > 4) {

		return arrheniusLoopEnergy(sequence, length);

	} else {

		return 0.0;

	}

}

// FD: April 28 2017
// FD: Adding initialization penalty when side length is zero, and
// Fd: only when there is an extension (single stranded or stack) on either side.
double EnergyModel::initializationPenalty(int length, int loop, int size) {

	double output = 0.0;

	// FD -- note that an openloop containing   -:CC:T    A:CC:A      T:CC   counts as size TWO not THREE (three sides, two joins)
	if ((loop > 0) && (loop < size)) {

		// not adjusting for temperature, hardcoded for now, etc.

		if (length == 0) {

			if (debugTraces) {
				cout << "Adding initalization penalty " << INIT_PENALTY << " -- length = " << length << " --loop = " << loop << " --   size " << size << endl;

			}
			return INIT_PENALTY;
		}

		if (length == 1) {

			if (debugTraces) {
				cout << "Adding initalization penalty/2- length = " << length << " --loop = " << loop << " --   size" << size << endl;
			}
			return (INIT_PENALTY / 2.0);
		}

	}

	return output;

}

double EnergyModel::arrheniusLoopEnergy(char* seq, int length) {

	double output = 0.0;

	for (int i = 0; i < (length - 1); i++) {

		int myMult = seq[i] * seq[i + 1];

		switch (myMult) {

		case baseA * baseA:
			output += (simOptions->energyOptions->dHA - simOptions->energyOptions->getTemperature() * (simOptions->energyOptions->dSA / 1000.0));
			break;
		}

	}

	return output;

}

double EnergyModel::saltCorrection(int size) {

// FD: Nupack makes a distinction between long (>20nt) and short domains.
// FD: For short domains, magnesium correction is not used. See computeSaltCorrection in utils/init.c for NUPACK 3.0.4.
// FD: In multistrand we don't set this distinction.

	return 0.368 * (size - 1) * log(simOptions->energyOptions->sodium + 3.3 * sqrt(simOptions->energyOptions->magnesium));

}

//double EnergyModel::ArrheniusLoopEnergy(char* seq, int size) {
//
//	double output = 0.0;
//
//	for (int i = 0; i < size; i++) {
//
//		switch (seq[i]) {
//
//		case BASE_A:
//			output += (simOptions->energyOptions->dSA);
//			break;
//		case BASE_C:
//			output += (simOptions->energyOptions->dSC);
//			break;
//		case BASE_G:
//			output += (simOptions->energyOptions->dSG);
//			break;
//		case BASE_T:
//			output += (simOptions->energyOptions->dST);
//			break;
//		}
//
//	}
//
//	return -output * simOptions->energyOptions->getTemperature() / 1000.0;
//
//}

int pairs[5] = { 0, 0, 0, 0, 0 };
int pairtypes[5][5] = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 } };
int basepair_sw[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

int lookuphelper[26] = { 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0 };		// A C G T    1 2 3 4
//                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z

// // helper function to convert to numerical base format.
int baseLookup(char base) {
	char temp = toupper(base);
	if (temp < 'A' || temp > 'Z')
		return base;
	return lookuphelper[temp - 'A'];
}

