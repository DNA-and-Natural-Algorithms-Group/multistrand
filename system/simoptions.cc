/*
 * simOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: hazel
 */

#include "options.h"
#include "simoptions.h"


SimOptions::SimOptions(void){

	// empty constructor

}


SimOptions::~SimOptions(void){

	// empty deconstructor

}




PSimOptions::PSimOptions(void) : PSimOptions(NULL){
	// Delegated constructor

}


PSimOptions::PSimOptions(PyObject *input){

	python_settings = input;
	simulation_mode = NULL;
	simulation_count = NULL;

}



long PSimOptions::getSimulationMode(void) {

	if(simulation_mode==NULL){

	  getLongAttr(python_settings, simulation_mode, &simulation_mode);


	}

	return simulation_mode;

}


long PSimOptions::getSimulationCount(void) {


	if(simulation_count==NULL){

		getLongAttr(python_settings, num_simulations, &simulation_count);


	}

	return simulation_count;

}




//#endif /* SYSTEM_SIMOPTIONS_H_ */
