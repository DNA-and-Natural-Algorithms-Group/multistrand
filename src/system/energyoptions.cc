/*
Copyright (c) 2017 California Institute of Technology. All rights reserved.
Multistrand nucleic acid kinetic simulator
help@multistrand.org
*/

/*
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

#include "options.h"	 	// python options helper
#include "energyoptions.h"
#include "moveutil.h"
#include "sequtil.h"

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

	return (3 == kinetic_rate_method);

}

double EnergyOptions::getBiScale(void) {

	return biScale;

}

double EnergyOptions::getUniScale(void) {

	if (usingArrhenius()) {
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

	getDoubleAttr(python_settings, join_concentration, &joinConcentration);


	if(!usingArrhenius() && biScale < 0.0){
		cout << "Warning! bimolecular_scaling is unset or negative!" << endl;
		cout <<	"Please set them using multistrand::utils::XPMetropolis37() or similar,"<< endl;
		cout <<	"or set Options.bimolecular_scaling directly."<< endl;
		cout << "Reverting to bimolecular_scaling = 5.0e6"<< endl;

		biScale = 1400000.0;

	}

	if(!usingArrhenius() && uniScale < 0.0){
		cout << "Warning! unimolecular_scaling is unset or negative!"<< endl;
		cout <<	"Please set them using multistrand::utils::XPMetropolis37() or similar,"<< endl;
		cout <<	"or set Options.bimolecular_scaling directly."<< endl;
		cout << "Reverting to unimolecular_scaling = 5.0e6"<< endl;

		uniScale = 5000000.0;

	}

	if (usingArrhenius()){

		uniScale = 1.0;

		getDoubleAttr(python_settings, lnAEnd, &AEnd);
		AValues[endMove] = AEnd;

		getDoubleAttr(python_settings, lnALoop, &ALoop);
		AValues[loopMove] = ALoop;

		getDoubleAttr(python_settings, lnAStack, &AStack);
		AValues[stackMove] = AStack;

		getDoubleAttr(python_settings, lnAStackStack, &AStackStack);
		AValues[stackStackMove] = AStackStack;

		getDoubleAttr(python_settings, lnALoopEnd, &ALoopEnd);
		AValues[loopEndMove] = ALoopEnd;

		getDoubleAttr(python_settings, lnAStackEnd, &AStackEnd);
		AValues[stackEndMove] = AStackEnd;

		getDoubleAttr(python_settings, lnAStackLoop, &AStackLoop);
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
		getDoubleAttr(python_settings, dHA, &dHA);

		// also loading four constants for entropy of nucleotide chain

	}

	// ionic conditions
	getDoubleAttr(python_settings, sodium, &sodium);
	getDoubleAttr(python_settings, magnesium, &magnesium);

	if (magnesium < 0.00 || magnesium > 0.2) {

		cout << "Magnesium concentration (" << magnesium << " M) is out of bounds (0.0 M - 0.2 M). Setting Na+/Mg2+ to 1.0 M / 0.0 M" << endl;
		sodium = 1.0;
		magnesium = 0.0;

	}

	if (sodium < 0.01 || sodium > 1.3) {

		cout << "Sodium concentration (" << sodium << " M) is out of bounds (0.01 M - 1.3 M). Setting Na+/Mg2+ to 1.0 M / 0.0 M" << endl;
		sodium = 1.0;
		magnesium = 0.0;

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
	gtenable = 1;
	kinetic_rate_method = 1;

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

