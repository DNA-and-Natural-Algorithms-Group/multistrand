/*
 * energyOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

#include "options.h"	 	// python options helper
#include "energyoptions.h"
#include "moveutil.h"

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>



using std::vector;
using std::string;



//int EnergyOptions::ARRTYPEF(MoveType left, MoveType right) {
//
//	return (int) (EnergyOptions::valuesPrime[left] * EnergyOptions::valuesPrime[right]);
//
//}


void EnergyOptions::initializeArrheniusConstants(void) {

	// unused, for now.
	// this is used for import/export of the constants.

}

double EnergyOptions::getTemperature(void) {

	return temperature;

}

long EnergyOptions::getDangles(void) {

	return dangles;

}

long EnergyOptions::getLogml(void) {

	return logml;

}
bool EnergyOptions::getGtenable(void) {

	return gtenable;

}

long EnergyOptions::getKineticRateMethod(void) {

	return kinetic_rate_method;

}

double EnergyOptions::getJoinConcentration(void) {

	return joinConcentration;

}

bool EnergyOptions::usingArrhenius(void) {

	return useArrRates;

}

double EnergyOptions::getBiScale(void) {

	return biScale;

}

double EnergyOptions::getUniScale(void) {

	if (useArrRates) {
		return 1.0;
	} else {
		return uniScale;
	}

}

string EnergyOptions::toString(void) {

// not sure if these are long
	long substrate_type = NULL;

	std::stringstream ss;

	ss << "temperature = " << temperature << " \n";
	ss << "dangles = " << dangles << " \n";
	ss << "logml = " << logml << " \n";
	ss << "gtenable = " << gtenable << " \n";
	ss << "kinetic_rate_method = " << kinetic_rate_method << " \n";
	ss << "joinConcentration = " << joinConcentration << " \n";
	ss << "biScale = " << biScale << " \n";
	ss << "uniScale = " << uniScale << " \n";
	ss << " substrate_type = " << substrate_type << " \n";

	string output = ss.str();

	return output;

}

string EnergyOptions::primeRateToString(double rate) {

	std::stringstream ss;

	if (fmod(rate, 1.0) > 0.49) {		// dealing with rate / 2 case

		rate = 2 * rate;

	}

	int myRate = round(rate);

	for (int i = 0; i < MOVETYPE_SIZE; i++) {

		int myPrime = moveutil::valuesPrime[i];

		if ((myRate % myPrime) == 0) {

			ss << moveutil::MoveToString[i] << ", ";

			if ((myRate % (myPrime * myPrime) == 0)) {

				ss << moveutil::MoveToString[i] << ", ";

			}

		}

	}

	return ss.str();

}

PEnergyOptions::PEnergyOptions(PyObject* input) :
		EnergyOptions() {

// extended constructor, inherits from regular energyOptions

	python_settings = input;

	getDoubleAttr(python_settings, temperature, &temperature);
	getLongAttr(python_settings, dangles, &dangles);
	getLongAttr(python_settings, log_ml, &logml);
	getBoolAttr(python_settings, gt_enable, &gtenable);
	getLongAttr(python_settings, rate_method, &kinetic_rate_method);

	getDoubleAttr(python_settings, bimolecular_scaling, &biScale);
	getDoubleAttr(python_settings, unimolecular_scaling, &uniScale);

// not sure if these are for both models or not;
	getDoubleAttr(python_settings, join_concentration, &joinConcentration);

	getBoolAttr(python_settings, useArrRates, &useArrRates);

	if (useArrRates) {

		uniScale = 1.0;

//		cout << "Loading Arrhenius constants -- do not forget to set AEnd, Eend, and so on. \n";

		getDoubleAttr(python_settings, AEnd, &AEnd);
		AValues[endMove] = AEnd;

		getDoubleAttr(python_settings, ALoop, &ALoop);
		AValues[loopMove] = ALoop;

		getDoubleAttr(python_settings, AStack, &AStack);
		AValues[stackMove] = AStack;

		getDoubleAttr(python_settings, AStackStack, &AStackStack);
		AValues[stackStackMove] = AStackStack;

		getDoubleAttr(python_settings, ALoopEnd, &ALoopEnd);
		AValues[loopEndMove] = ALoopEnd;

		getDoubleAttr(python_settings, AStackEnd, &AStackEnd);
		AValues[stackEndMove] = AStackEnd;

		getDoubleAttr(python_settings, AStackLoop, &AStackLoop);
		AValues[stackLoopMove] = AStackLoop;

		getDoubleAttr(python_settings, EEnd, &EEnd);
		EValues[endMove] = EEnd;

		getDoubleAttr(python_settings, ELoop, &ELoop);
		EValues[loopMove] = ELoop;

		getDoubleAttr(python_settings, EStack, &EStack);
		EValues[stackMove] = EStack;

		getDoubleAttr(python_settings, EStackStack, &EStackStack);
		EValues[stackStackMove] = EStackStack;

		getDoubleAttr(python_settings, ELoopEnd, &ELoopEnd);
		EValues[loopEndMove] = ELoopEnd;

		getDoubleAttr(python_settings, EStackEnd, &EStackEnd);
		EValues[stackEndMove] = EStackEnd;

		getDoubleAttr(python_settings, EStackLoop, &EStackLoop);
		EValues[stackLoopMove] = EStackLoop;

		getDoubleAttr(python_settings, dSA, &dSA);
		getDoubleAttr(python_settings, dST, &dST);
		getDoubleAttr(python_settings, dSC, &dSC);
		getDoubleAttr(python_settings, dSG, &dSG);

		getDoubleAttr(python_settings, alpha, &alpha);

		// also loading four constants for entropy of nucleotide chain

	}

}

bool PEnergyOptions::compareSubstrateType(long type) {

	return testLongAttr(python_settings, substrate_type, =, type);

}

void PEnergyOptions::getParameterFile(char* input, PyObject* tempString) {

	input = (char *) getStringAttr(python_settings, parameter_file, tempString);

}

// CENERGYOPTIONS

CEnergyOptions::CEnergyOptions() :
		EnergyOptions() {

// extended constructor, inherits from regular energyOptions

// some basic default values;
	temperature = 310.15;
	dangles = 1;
	logml = 0;
	gtenable = 0;
	kinetic_rate_method = 2;

	biScale = 1.38e+06;
	uniScale = 1.5e+08;

	joinConcentration = 1e-06;

}

bool CEnergyOptions::compareSubstrateType(long type) {

	return (SUBSTRATE_DNA == type); // only test TRUE for when substrate is DNA.

}

// functionality to specify an energy file in a custom place for non-DNA/RNA substrates. unused atm.
void CEnergyOptions::getParameterFile(char* input, PyObject* tempString) {

	input = "";
}

