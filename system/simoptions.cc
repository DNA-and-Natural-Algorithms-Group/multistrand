/*
 * simOptions.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

#include "options.h"
#include "ssystem.h"
#include "simoptions.h"
#include "scomplex.h"

#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

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
	o_interval = NULL;
	o_time = NULL;
	stop_options = NULL;
	stop_count = NULL;
	max_sim_time = NULL;
	myComplexes = NULL;

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

		getLongAttr(python_settings, output_interval, &o_interval);

		if (debug) {

			printf("The o interval is %li \n", o_interval);

		}

	}

	return o_interval;

}

double PSimOptions::getOTime(void) {

	if (o_time == NULL) {

		getDoubleAttr(python_settings, output_time, &o_time);

		if (debug) {

			printf("The output time is %li \n", o_time);

		}

	}

	return o_time;

}

void PSimOptions::incrementTrajectoryCount(void) {

	if (python_settings != NULL) {
		pingAttr(python_settings, increment_trajectory_count);
	}

}

long PSimOptions::getStopOptions(void) {

	if (stop_options == NULL) {

		getLongAttr(python_settings, use_stop_conditions, &stop_options);

		if (debug) {

			printf("The stop option is %li \n", stop_options);

		}

	}

	return stop_options;

}

long PSimOptions::getStopCount(void) {

	if (stop_count == NULL) {

		getLongAttr(python_settings, stop_count, &stop_count);

		if (debug) {

			printf("The stop count is %li \n", stop_count);

		}

	}

	return stop_count;

}

double PSimOptions::getMaxSimTime(void) {

	if (max_sim_time == NULL) {

		getDoubleAttr(python_settings, simulation_time, &max_sim_time);

		if (debug) {

			printf("The max sim time is %d \n", max_sim_time);

		}

	}

	return max_sim_time;

}

void PSimOptions::sendTransitionInfo(PyObject *transition_tuple) {

	if (python_settings != NULL) {

		pushTransitionInfo(python_settings, transition_tuple);

	}

}

void PSimOptions::storeStrandComplex(char* input1, char* input2, char* input3) {

	if (myComplexes == NULL) {

		myComplexes = new std::vector<complex_input>();

	}

	complex_input myInfo = complex_input(input1, input2, input3);

	myComplexes->push_back(myInfo);

}

PyObject* PSimOptions::getPythonSettings() {

	return python_settings;

}

