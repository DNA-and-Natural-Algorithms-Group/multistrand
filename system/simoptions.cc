/*
 * simOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: hazel
 */

#include "simOptions.h"
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

		return getLongAttr(system_options, simulation_mode, &simulation_mode );

};



#endif /* SYSTEM_SIMOPTIONS_H_ */
