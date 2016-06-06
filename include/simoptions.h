/*
 * m_options.h
 *
 *  Created on: Jun 5, 2016
 *      Author: hazel
 */

#ifndef SYSTEM_SIMOPTIONS_H_
#define SYSTEM_SIMOPTIONS_H_

#include <python2.7/Python.h>

class SimOptions {
public:

	// Constructors
	SimOptions(void);

	// Virtual methods
	virtual ~SimOptions(void);
	virtual long getSimulationMode(void) = 0;
	virtual void functionTwo(void) = 0;

};



class PSimOptions: public SimOptions {
public:
	//constructors
	PSimOptions(void);
	PSimOptions(PyObject);

	long getSimulationMode(void);
	void functionTwo(void);

protected:
	PyObject *python_settings;
	long simulation_mode;
};

#endif /* SYSTEM_SIMOPTIONS_H_ */
