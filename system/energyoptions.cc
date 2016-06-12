/*
 * energyOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

#include "options.h"	 // python options helper
#include "energyoptions.h"
#include "energyoptions.h"

#include <vector>
#include <iostream>
#include <string>

using std::vector;
using std::string;

EnergyOptions::EnergyOptions(void) {

// empty constructor

}

EnergyOptions::~EnergyOptions(void) {

// empty deconstructor

}

PEnergyOptions::PEnergyOptions(void) :
		PEnergyOptions(NULL) {

// passed-on constructor

}

PEnergyOptions::PEnergyOptions(PyObject* input) {

// empty constructor

	python_settings = input;

	temperature = NULL;
	dangles  = NULL;
	logml = NULL;
	gtenable = NULL;
	kinetic_rate_method = NULL;

	getDoubleAttr(python_settings, temperature, &temperature);
	getLongAttr(python_settings, dangles, &dangles);
	getLongAttr(python_settings, log_ml, &logml);
	getBoolAttr(python_settings, gt_enable, &gtenable);
	getLongAttr(python_settings, rate_method, &kinetic_rate_method);

}

