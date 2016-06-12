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
}

EnergyOptions::~EnergyOptions(void) {

// empty deconstructor

}

PEnergyOptions::PEnergyOptions(void) :
		PEnergyOptions(NULL) {

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

}

