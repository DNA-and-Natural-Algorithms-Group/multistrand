/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include "options.h"	 	// python options helper
#include "energyoptions.h"
#include "moveutil.h"
#include "sequtil.h"

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

	return (kinetic_rate_method == RATE_METHOD_ARRHENIUS);

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

	std::stringstream ss;

	ss << "temperature = " << temperature << " \n";
	ss << "dangles = " << dangles << " \n";
	ss << "logml = " << logml << " \n";
	ss << "gtenable = " << gtenable << " \n";
	ss << "kinetic_rate_method = " << kinetic_rate_method << " \n";
	ss << "joinConcentration = " << joinConcentration << " \n";
	ss << "biScale = " << biScale << " \n";
	ss << "uniScale = " << uniScale << " \n";

	string output = ss.str();

	return output;

}

PEnergyOptions::PEnergyOptions(PyObject* input) :
		EnergyOptions() {

	// extended constructor, inherits from regular energyOptions

	python_settings = input;

	getDoubleAttrReify(python_settings, temperature, &temperature);
	getLongAttrReify(python_settings, dangles, &dangles);
	getLongAttr(python_settings, log_ml, &logml);
	getBoolAttr(python_settings, gt_enable, &gtenable);
	getLongAttr(python_settings, rate_method, &kinetic_rate_method);

	getDoubleAttrReify(python_settings, bimolecular_scaling, &biScale);
	getDoubleAttrReify(python_settings, unimolecular_scaling, &uniScale);
	getDoubleAttr(python_settings, join_concentration, &joinConcentration);


	if(!usingArrhenius() && biScale < 0.0){
		cout << "Warning! `bimolecular_scaling` is unset or negative!" << endl;
		cout << "  Please set all parameters using `Options.DNA23Metropolis()` or similar,"<< endl;
		cout << "  or set `Options.bimolecular_scaling` directly."<< endl;
		cout << "  Reverting to `bimolecular_scaling = 1.38e+6` from `Options.JSDefault()`."<< endl;

		kinetic_rate_method = RATE_METHOD_KAWASAKI;
		biScale = 1.38e+6;
	}
	if(!usingArrhenius() && uniScale < 0.0){
		cout << "Warning! `unimolecular_scaling` is unset or negative!"<< endl;
		cout << "  Please set all parameters using `Options.DNA23Metropolis()` or similar,"<< endl;
		cout << "  or set `Options.unimolecular_scaling` directly."<< endl;
		cout << "  Reverting to `unimolecular_scaling = 1.50e+08` from `Options.JSDefault()`."<< endl;

		kinetic_rate_method = RATE_METHOD_KAWASAKI;
		uniScale = 1.50e+08;
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
	getDoubleAttrReify(python_settings, sodium, &sodium);
	getDoubleAttrReify(python_settings, magnesium, &magnesium);

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

bool PEnergyOptions::compareSubstrateType(int type) {

	return testLongAttr(python_settings, substrate_type, =, type);

}

void PEnergyOptions::getParameterFile(char* input, PyObject* tempString) {

	input = (char *) getStringAttrReify(python_settings, parameter_file, tempString);

}
