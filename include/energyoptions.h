/*
 * m_options.h
 *
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

/* EnergyOption class header. This is currently a wrapper for the options set in python. */

#ifndef __ENERGYOPTIONS_H_
#define __ENERGYOPTIONS_H_

#include <python2.7/Python.h>
#include "ssystem.h"
#include <vector>
#include <string>
#include <iostream>

#include "utility.h"

using std::vector;
using std::string;
using namespace utility;

class EnergyOptions {
public:

	// Constructors
	EnergyOptions(void);
	virtual ~EnergyOptions(void);

	// non-virtual getters
//	long getSimulationMode(void);
	// IO Methods
//	string toString(void);

	// actual option values

protected:

	double temperature = NULL;
	long dangles  = NULL;
	long logml = NULL;
	bool gtenable = NULL;
	long kinetic_rate_method = NULL;

};

class PEnergyOptions: public EnergyOptions {
public:
	//constructors
	PEnergyOptions(void);
	PEnergyOptions(PyObject*);

protected:
	PyObject* python_settings;


};

#endif
