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

class SimOptions {
public:

	// Constructors
	SimOptions(void);

	// Structs
	struct complex_input {
		char *sequence;
		char *structure;
		identlist *list;
		complex_input(char *string1, char *string2, identlist *list1) {

			printf("Making new struct \n");
//			printf(string1);
//			printf(string2);

			sequence = string1;
			structure = string2;
			list = list1;

			printf("Done making new struct \n");

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
	virtual std::vector<complex_input> getComplexes(PyObject *, long) = 0;



	// actual option values
protected:
	long simulation_mode;
	long simulation_count;
	long o_interval;
	double o_time;
	long stop_options;
	long stop_count;
	double max_sim_time;
	std::vector<complex_input> *myComplexes;
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
	std::vector<complex_input> getComplexes(PyObject *alternate_start,
			long current_seed);

protected:
	bool debug;
	PyObject *python_settings;

};

#endif

