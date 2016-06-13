/*
 * m_options.h
 *
 *  Created on: Jun 12, 2016
 *      Author: Frits Dannenberg
 */

/* EnergyOption class header. This is currently a wrapper for the options set in python. */

#ifndef __ENERGYOPTIONS_H_
#define __ENERGYOPTIONS_H_

#include <python2.7/Python.h>
#include "ssystem.h"
#include "simoptions.h"
#include <vector>
#include <string>
#include <iostream>

#include "utility.h"
//
//using std::vector;
//using std::string;
//using namespace utility;

class EnergyOptions {
public:
	// Constructors
	EnergyOptions(void);
	~EnergyOptions(void);

	// non-virtual getters
	double getTemperature(void);
	long getDangles(void);
	long getLogml(void);
	bool getGtenable(void);
	long getKineticRateMethod(void);

	// virtual
	virtual bool compareSubstrateType(long) =0;

protected:

	double temperature = NULL;
	long dangles = NULL;
	long logml = NULL;
	bool gtenable = NULL;
	long kinetic_rate_method = NULL;

	// not sure if these are long
	long substrate_type = NULL;

};

class PEnergyOptions: public EnergyOptions {
public:
	//constructors
	//PEnergyOptions(void);
	PEnergyOptions(PyObject*);

	// implemented virtual
	bool compareSubstrateType(long);


protected:
	PyObject* python_settings;

};

#endif
