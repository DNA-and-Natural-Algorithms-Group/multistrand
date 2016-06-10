/*
 * m_options.h
 *
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

/* SimOption class header. This is currently a wrapper for the options set in python. */

#ifndef __SIMOPTIONS_H_
#define __SIMOPTIONS_H_

#include <python2.7/Python.h>
#include "ssystem.h"
#include <vector>
#include <string>
#include <iostream>

#include "utility.h"

using namespace std;
using namespace utility;




class SimOptions {
public:

	// Constructors
	SimOptions(void);
	virtual ~SimOptions(void);

	// Virtual methods
	virtual long getSimulationMode(void) = 0;
	virtual long getSimulationCount(void) = 0;
	virtual long getOInterval(void) = 0;
	virtual double getOTime(void) = 0;
	virtual void incrementTrajectoryCount(void) = 0;	// PyObject compliance
	virtual long getStopOptions(void) = 0;
	virtual long getStopCount(void) = 0;
	virtual double getMaxSimTime(void) = 0;
	virtual void sendTransitionInfo(PyObject*) = 0; // PyObject compliance
	//virtual void storeStrandComplex(char*, char*, identlist*) = 0; // backwards way of refactoring
	virtual PyObject* getPythonSettings(void) = 0;
	virtual void getComplexes(
			vector<complex_input>&, PyObject*, long) = 0;

	// actual option values
protected:
	long simulation_mode;
	long simulation_count;
	long o_interval;
	double o_time;
	long stop_options;
	long stop_count;
	double max_sim_time;
	vector<complex_input>* myComplexes;
	long seed;
};

class PSimOptions: public SimOptions {
public:
	//constructors
	PSimOptions(void);
	PSimOptions(PyObject *system_options);

	// Implemented virtual methods
	long getSimulationMode(void);
	long getSimulationCount(void);
	long getOInterval(void);
	double getOTime(void);
	void incrementTrajectoryCount(void);		// PyObject compliance
	long getStopOptions(void);
	long getStopCount(void);
	double getMaxSimTime(void);
	void sendTransitionInfo(PyObject *transitions);
	//void storeStrandComplex(char* input1, char* input2, identlist* input3);
	PyObject* getPythonSettings(void);
	void getComplexes(
			vector<complex_input>& myComplexes, PyObject *alternate_start,
			long current_seed);

protected:
	bool debug;
	PyObject *python_settings;

};

#endif

