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

class EnergyOptions {
public:
	// Constructors
	EnergyOptions(void);
	virtual ~EnergyOptions(void) = 0;

	// non-virtual getters
	double getTemperature(void);
	long getDangles(void);
	long getLogml(void);
	bool getGtenable(void);
	long getKineticRateMethod(void);
	double getJoinConcentration(void);

	double getBiScale(void);
	double getUniScale(void);

	string toString(void);

	// virtual
	virtual bool compareSubstrateType(long) =0;
	virtual void getParameterFile(char*, PyObject*) = 0;

protected:

	double temperature = NULL;
	long dangles = NULL;
	long logml = NULL;
	bool gtenable = NULL;
	long kinetic_rate_method = NULL;
	double joinConcentration = NULL;

	double biScale = NULL;
	double uniScale = NULL;

	// not sure if these are long
	long substrate_type = NULL;

};

class PEnergyOptions: public EnergyOptions {
public:
	//constructors
	PEnergyOptions(PyObject*);
	~PEnergyOptions();

	// implemented virtual
	bool compareSubstrateType(long);
	void getParameterFile(char*, PyObject*);

protected:
	PyObject* python_settings;
};

class CEnergyOptions: public EnergyOptions {
public:
	// constructors
	CEnergyOptions();
	~CEnergyOptions();

	// implemented virtual
	bool compareSubstrateType(long);
	void getParameterFile(char*, PyObject*);

protected:
	// empty

};

#endif
