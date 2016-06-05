/*
 * m_options.h
 *
 *  Created on: Jun 5, 2016
 *      Author: hazel
 */

#ifndef SYSTEM_SIMOPTIONS_H_
#define SYSTEM_SIMOPTIONS_H_

class SimOptions {
public:

	// Constructors
	SimOptions(void);

	// Virtual methods
	virtual ~SimOptions(void);
	virtual void getSimulationMode(void) = 0;
	virtual void functionTwo(void) = 0;

};



class PySimOptions: public SimOptions {
public:
	void PySimOptions(void);

	//void ~PySimOptions(void);
	void getSimulationMode(void);
	void functionTwo(void);

protected:
	PyObject python_settings;
};

#endif /* SYSTEM_SIMOPTIONS_H_ */
