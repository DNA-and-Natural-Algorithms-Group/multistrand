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

class SimOptions {
public:

	// Constructors
	SimOptions(void);

	// Virtual methods
	virtual ~SimOptions(void);
	virtual long getSimulationMode(void) = 0;
	virtual long getSimulationCount(void) = 0;

};



class PSimOptions: public SimOptions {
public:
	//constructors
	PSimOptions(void);
	PSimOptions(PyObject *system_options);

	long getSimulationMode(void);
	long getSimulationCount(void);

protected:
	PyObject *python_settings;
	long simulation_mode;
	long simulation_count;
};

#endif
