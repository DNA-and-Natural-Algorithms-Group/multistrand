/*
 * simOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: hazel
 */

#include "options.h"
#include "ssystem.h"
#include "simoptions.h"

SimOptions::SimOptions(void) {

	// empty constructor

}

SimOptions::~SimOptions(void) {

	// empty deconstructor

}

PSimOptions::PSimOptions(void) :
		PSimOptions(NULL) {
	// Delegated constructor

}

PSimOptions::PSimOptions(PyObject *input) {

	python_settings = input;
	simulation_mode = NULL;
	simulation_count = NULL;
	debug = true;	// this is the main switch for simOptions debug, for now.

}

long PSimOptions::getSimulationMode(void) {

	if (simulation_mode == NULL) {

		getLongAttr(python_settings, simulation_mode, &simulation_mode);

		if (debug) {

			printf("The simulation mode is %li \n", simulation_mode);

		}

	}

	return simulation_mode;

}

long PSimOptions::getSimulationCount(void) {

	if (simulation_count == NULL) {

		getLongAttr(python_settings, num_simulations, &simulation_count);

		if (debug) {

			printf("The simulation count is %li \n", simulation_count);

		}

	}

	return simulation_count;

}

long PSimOptions::getOInterval(void) {

	if (o_interval == NULL) {

		getLongAttr(python_settings, output_interval,&o_interval);

		if (debug) {

			printf("The o interval is %li \n", o_interval);

		}

	}

	return simulation_count;

}


double PSimOptions::getOTime(void) {

	if (o_time == NULL) {

		getDoubleAttr(python_settings, output_time,&o_time);

		if (debug) {

			printf("The output time is %li \n", o_time);

		}

	}

	return simulation_count;

}


void PSimOptions::incrementTrajectoryCount(void) {

	if (python_settings != NULL) {
		//pingAttr(python_settings, increment_trajectory_count);
	}

}

//#endif /* SYSTEM_SIMOPTIONS_H_ */
