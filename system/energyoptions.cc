/*
 * energyOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

#include "options.h"	 // python options helper
#include "energyoptions.h"

#include <vector>
#include <iostream>
#include <string>

using std::vector;
using std::string;

EnergyOptions::EnergyOptions(void) {

	temperature = NULL;
	dangles = NULL;
	logml = NULL;
	gtenable = NULL;
	kinetic_rate_method = NULL;
	substrate_type = NULL;

}

EnergyOptions::~EnergyOptions(void) {

// empty deconstructor

}
// easy getters

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

double EnergyOptions::getBiScale(void) {

	return biScale;

}

double EnergyOptions::getUniScale(void) {

	return uniScale;

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

void PEnergyOptions::getParameterFile(char* input) {

	input = (char *) getStringAttr(python_settings, parameter_file, NULL);

}

