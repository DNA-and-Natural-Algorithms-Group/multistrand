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




PSimOptions::PSimOptions(void){

	python_settings = NULL;
	simulation_mode = NULL;
	// empty constructor

}


PSimOptions::PSimOptions(PyObject input){

	python_settings = *input;
	simulation_mode = NULL;

}



long PSimOptions::getSimulationMode(void) {

	if(simulation_mode==NULL){

	  //getLongAttr(python_settings, simulation_mode, &simulation_mode);
	  //getLongAttr(python_settings, simulation_mode, &simulation_mode );


	}

	return simulation_mode;

}


void PSimOptions::functionTwo(void) {



}




//#endif /* SYSTEM_SIMOPTIONS_H_ */
