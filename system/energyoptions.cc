/*
 * energyOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

#include "options.h"	 	// python options helper
#include "energyoptions.h"

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

using std::vector;
using std::string;

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

	return useArrhenius;

}

double EnergyOptions::getBiScale(void) {

	return biScale;

}

double EnergyOptions::getUniScale(void) {

	return uniScale;

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

	ss << rate << " - or: ";

	if (fmod(rate, 1.0) > 0.49) {

		rate = 2 * rate;  // dealing with rate / 2 case

	}

	int myRate = round(rate);

	for (int i = 0; i < MOVETYPE_SIZE; i++) {

		if (rate == round(AValues[i])) {

			myRate = myRate * 2; // the combination contains prime factor 2; fix this.

		}

	}

	for (int i = 0; i < MOVETYPE_SIZE; i++) {

		int myPrime = round(AValues[i]);

		if ((myRate % myPrime) == 0) {

			ss << MoveToString[i] << ", ";

			if ((myRate % (myPrime * myPrime) == 0)) {

				ss << MoveToString[i] << ", ";

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

