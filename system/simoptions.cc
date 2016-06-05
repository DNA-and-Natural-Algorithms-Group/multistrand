/*
 * simOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: hazel
 */

#include "simoptions.h"
#include "options.h"

SimOptions::SimOptions(void){

	// empty constructor

}


SimOptions::~SimOptions(void){

	// empty deconstructor

}



PySimOptions::PySimOptions(PyObject input){

	python_settings = input;

}



long PySimOptions::getSimulationMode(void) {

		long output; // = getLongAttr(system_options, simulation_mode, &simulation_mode );

		output = 0.0;

		return output;

}



//#endif /* SYSTEM_SIMOPTIONS_H_ */
