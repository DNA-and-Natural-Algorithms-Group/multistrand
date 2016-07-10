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

double EnergyModel::applyPrefactors(MoveType left, MoveType right) {

	return arrheniusRates[left * MOVETYPE_SIZE + right];

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

