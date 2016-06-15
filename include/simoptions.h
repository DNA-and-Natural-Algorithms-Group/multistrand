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

#include "energyoptions.h"
#include "utility.h"

using std::vector;
using std::string;
using namespace utility;

class EnergyOptions;
// forward declare

class SimOptions {
public:

	// Constructors
	SimOptions(void);
	virtual ~SimOptions(void) =0;

	// Non-virtual
	bool useFixedRandomSeed();
	long getInitialSeed();
	EnergyOptions* getEnergyOptions();
	long getSimulationMode();
	long getSimulationCount();
	long getOInterval(void);
	double getOTime(void);
	long getStopOptions(void);
	long getStopCount(void);
	double getMaxSimTime(void);

	// Virtual methods

	virtual PyObject* getPythonSettings(void) = 0;
	virtual void generateComplexes(PyObject*, long) = 0;
	virtual stopComplexes* getStopComplexes(int) = 0;

	// Exit signalling
	virtual void stopResultError(long) = 0;
	virtual void stopResultNan(long) = 0;
	virtual void stopResultNormal(long, double, char*) = 0;
	virtual void stopResultTime(long, double) = 0;
	virtual void stopResultBimolecular(string, long, double, double, char*) = 0;

	// Pushing info back to Python
//	virtual void incrementTrajectoryCount(void) = 0;
//	virtual void sendTransitionInfo(PyObject*) = 0;
//	virtual void pushTrajectory(long, int, char*, char*, char*, double) = 0;
//	virtual void pushTrajectoryInf(double)=0;

// IO Methods
	string toString(void);

	// actual option values
	vector<complex_input>* myComplexes;

protected:
	long simulation_mode;
	long simulation_count;
	long o_interval;
	double o_time;
	long stop_options;
	long stop_count;
	double max_sim_time;
	long seed;
	stopComplexes* myStopComplexes;
	bool fixedRandomSeed;

	EnergyOptions* myEnergyOptions;
};

class PSimOptions: public SimOptions {
public:
	//constructors
	PSimOptions(void);
	PSimOptions(PyObject *system_options);

	// Implemented virtual methods
	//long getSimulationMode(void);
	//long getSimulationCount(void);
	//long getOInterval(void);
	//double getOTime(void);
//	long getStopOptions(void);
//	long getStopCount(void);
//	double getMaxSimTime(void);
	PyObject* getPythonSettings(void);
	void generateComplexes(PyObject *alternate_start, long current_seed);
	stopComplexes* getStopComplexes(int);

	// Error signaling
	void stopResultError(long);
	void stopResultNan(long);
	void stopResultNormal(long, double, char*);
	void stopResultTime(long, double);
	void stopResultBimolecular(string, long, double, double, char*);

	// Push back to Python
//	void incrementTrajectoryCount(void);		// PyObject compliance
//	void sendTransitionInfo(PyObject *transitions);
//	void pushTrajectory(long, int, char*, char*, char*, double);
//	void pushTrajectoryInf(double);

protected:
	bool debug;
	PyObject *python_settings;

};

class CSimOptions: public SimOptions {
public:
	//constructors
	CSimOptions();

	PyObject* getPythonSettings(void);
	void generateComplexes(PyObject *alternate_start, long current_seed);
	stopComplexes* getStopComplexes(int);

	// Error signaling
	void stopResultError(long);
	void stopResultNan(long);
	void stopResultNormal(long, double, char*);
	void stopResultTime(long, double);
	void stopResultBimolecular(string, long, double, double, char*);

protected:
	bool debug;
	PyObject *python_settings = NULL;

};

#endif

