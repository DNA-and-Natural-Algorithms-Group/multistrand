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

class SimOptions {
public:

	// Constructors
	SimOptions(void);

	// Structs
	struct complex_input {

		std::string sequence;
		std::string structure;
		std::string list;

		complex_input(char* string1, char* string2, char* list1) {

			std::string other1(string1);
			std::string other2(string2);
			std::string other3(list1);

			sequence = other1;
			structure = other2;
			list = other3;

		}
	};

	// Virtual methods
	virtual ~SimOptions(void);
	virtual long getSimulationMode(void) = 0;
	virtual long getSimulationCount(void) = 0;
	virtual long getOInterval(void) = 0;
	virtual double getOTime(void) = 0;
	virtual void incrementTrajectoryCount(void) = 0;	// PyObject compliance
	virtual long getStopOptions(void) = 0;
	virtual long getStopCount(void) = 0;
	virtual double getMaxSimTime(void) = 0;


	virtual void sendTransitionInfo(PyObject*) = 0; // PyObject compliance
	virtual void storeStrandComplex(char*, char*, char*) = 0; 	// backwards way of refactoring
	virtual PyObject* getPythonSettings(void) = 0;

	// actual option values
protected:
	long simulation_mode;
	long simulation_count;
	long o_interval;
	double o_time;
	long stop_options;
	long stop_count;
	double max_sim_time;
	std::vector<complex_input>* myComplexes;
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
	void storeStrandComplex(char* input1, char* input2, char* input3 );
	PyObject* getPythonSettings(void);

protected:
	bool debug;
	PyObject *python_settings;

};

#endif

