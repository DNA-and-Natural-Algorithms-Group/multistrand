/*
 Copyright (c) 2007-2010 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 */

#include "energymodel.h"
#include "simoptions.h"
#include "loop.h"

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

	return A * exp(E / (gasConstant * temperature));

}

void setArrheniusRate(double ratesArray[], EnergyOptions* options, double temperature, int left, int right) {

	double ALeft = options->AValues[left];
	double ELeft = options->EValues[left];

	double ARight = options->AValues[right];
	double ERight = options->EValues[right];

	double kLeft = expRate(ALeft, ELeft, temperature);
	double kRight = expRate(ARight, ERight, temperature);

	// DEBUG! Set all rates to be 1.0 for debugging.

//	ratesArray[left * MOVETYPE_SIZE + right] = kLeft * kRight;

	ratesArray[left * MOVETYPE_SIZE + right] = 1.0;

}

void EnergyModel::computeArrheniusRates(double temperature) {

	for (int i = 0; i < MOVETYPE_SIZE; i++) {

		for (int j = 0; j < MOVETYPE_SIZE; j++) {

			arrheniusRates[i * MOVETYPE_SIZE + j] = NULL;

			setArrheniusRate(arrheniusRates, simOptions->energyOptions, temperature, i, j);

		}

	}

}

double EnergyModel::applyPrefactors(MoveType left, MoveType right) {

	return arrheniusRates[left * MOVETYPE_SIZE + right];

}

//MoveType EnergyModel::getPrefactorsMulti(int index, Loop* loop, int sideLengths[]) {
//
//	return this->getPrefactorsMulti(index, loop->numAdjacent, sideLengths);
//
//}

// FD: A base pair is present between a stacking loop and a multi loop.
// FD: We query the local context of the middle pair;
// FD: this can be either a loop, stack+loop, or stack+stack situation.
MoveType EnergyModel::getPrefactorsMulti(int index, int numAdjacent, int sideLengths[]) {

// FD: if adjacent are not neighbored in this order, then we need to replace this code with the correct mapping.
	int rightStrand = (index + 1) % numAdjacent;

	return this->prefactorMulti(sideLengths[index], sideLengths[rightStrand]);

}

MoveType EnergyModel::prefactorMulti(int sideLength1, int sideLength2) {
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

		return stackLoopMove;

	} else if (!dangleOpen && !dangleClosed) {

		return endMove;

	} else if (dangleOpen) {

		return loopMove;

	} else {

		return stackEndMove;

	}

}

MoveType EnergyModel::prefactorOpen(int index, int numAdjacent, int sideLengths[]) {

	assert(index < (numAdjacent + 1));

	bool openOnLeft = (index == 0);
	bool openOnRight = (index + 1 == numAdjacent);

// each strand can either be non-existing, or single stranded, or double stranded.

	if (!openOnLeft && !openOnRight) { // this is the multi-loop case

		return prefactorMulti(sideLengths[index], sideLengths[index + 1]);

	} else if (openOnLeft && openOnLeft) { // this is the "end" case

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

double EnergyModel::ArrheniusLoopEnergy(char* seq, int size) {

	double output = 0.0;

	for (int i = 0; i < size; i++) {

		switch (seq[i]) {

		case BASE_A:
			output += (simOptions->energyOptions->dS_A);
			break;
		case BASE_C:
			output += (simOptions->energyOptions->dS_C);
			break;
		case BASE_G:
			output += (simOptions->energyOptions->dS_G);
			break;
		case BASE_T:
			output += (simOptions->energyOptions->dS_T);
			break;
		}

	}

	return -output * simOptions->energyOptions->getTemperature() / 1000.0;

}

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

