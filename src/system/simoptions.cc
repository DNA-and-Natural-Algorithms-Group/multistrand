/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */

/*
 *
 *  Created on: Jun 5, 2016
 *      Author: Frits Dannenberg
 */

#include "options.h"	 // python options helper
#include "ssystem.h"
#include "simoptions.h"
#include "energyoptions.h"
#include "scomplex.h"

#include <time.h>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>

using std::vector;
using std::string;

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

PSimOptions::PSimOptions(PyObject* input) :
		SimOptions() {
	// Inherited constructor.

	python_settings = input;

	// initializers calling python object -- these can use a super object getter.
	// Not clear at the moment if calling all settings is possible without crashing.
	getBoolAttr(python_settings, initial_seed_flag, &fixedRandomSeed);

	if (fixedRandomSeed) {

		getLongAttr(python_settings, initial_seed, &seed);

	}

	energyOptions = new PEnergyOptions(python_settings);

	getLongAttr(python_settings, simulation_mode, &simulation_mode);
	getLongAttr(python_settings, num_simulations, &simulation_count);
	getLongAttr(python_settings, output_interval, &o_interval);
	getDoubleAttr(python_settings, output_time, &o_time);
	getLongAttr(python_settings, stop_count, &stop_count);
	getLongAttr(python_settings, use_stop_conditions, &stop_options);
	getDoubleAttr(python_settings, simulation_time, &max_sim_time);

	// to print the initial state, used for statespace building
	getBoolAttr(python_settings, print_initial_first_step, &printInitialFirstStep);
	// to enable cotranscriptional folding, assumes a single strand
	getBoolAttr(python_settings, cotranscriptional, &cotranscriptional);
	getDoubleAttr(python_settings, cotranscriptional_rate, &cotranscriptional_rate);

	getLongAttr(python_settings, verbosity, &verbosity);
	getBoolAttr(python_settings, activestatespace, &statespaceActive);
	getBoolAttr(python_settings, reuse_energymodel, &reuseEnergyModel);
	getDoubleAttr(python_settings, ms_version, &ms_version);

	debug = false;	// this is the main switch for simOptions debug, for now.

}

string SimOptions::toString() {

	std::stringstream ss;

	ss << "simulation_mode = " << simulation_mode << " \n";
	ss << "simulation_count = " << simulation_count << " \n";
	ss << "o_interval = " << o_interval << " \n";
	ss << "o_time = " << o_time << " \n";
	ss << "stop_options = " << stop_options << " \n";
	ss << "stop_count = " << stop_count << " \n";
	ss << "max_sim_time = " << max_sim_time << " \n";
	ss << "seed = " << seed << " \n";

//	ss << "myComplexes = { ";
//
//	for (int i = 0; i < myComplexes->size(); i++) {
//
//		ss << "{ " << myComplexes->at(i).sequence << ", " << myComplexes->at(i).structure << " }";
//
//	}
//
//	ss << "} \n";
//
//	ss << " myStopComplexes = { ";
//
//	// linked list iterator
//	stopComplexes* myStopComplex = myStopComplexes; // copying pointer so we can iterate.
//
//	ss << "} \n";
//
	string output = ss.str();

	output += energyOptions->toString();

	return output;

}

bool SimOptions::useFixedRandomSeed() {

	return fixedRandomSeed;

}

long SimOptions::getSeed() {

	return seed;

}

EnergyOptions* SimOptions::getEnergyOptions() {

	return energyOptions;

}

long SimOptions::getSimulationMode(void) {

	return simulation_mode;

}

long SimOptions::getSimulationCount(void) {

	return simulation_count;

}

long SimOptions::getOInterval(void) {

	return o_interval;

}

double SimOptions::getOTime(void) {

	return o_time;

}

long SimOptions::getStopOptions(void) {

	return stop_options;

}

long SimOptions::getStopCount(void) {

	return stop_count;

}

double SimOptions::getMaxSimTime(void) {

	return max_sim_time;

}

bool SimOptions::getPrintIntialFirstStep() {

	return printInitialFirstStep;
}

bool SimOptions::usingArrhenius(void) {

	return energyOptions->usingArrhenius();

}

PyObject* PSimOptions::getPythonSettings() {

	return python_settings;

}

void PSimOptions::generateComplexes(PyObject *alternate_start, long current_seed) {

	myComplexes = new vector<complex_input>(0); // wipe the pointer to the previous object;

	PyObject *py_start_state = NULL, *py_complex = NULL;
	PyObject *py_seq = NULL, *py_struc = NULL;
	PyObject *py_err = NULL;

	if (myComplexes->size() == 0) {

		complex_input *tempcomplex = NULL;
		const char *sequence, *structure;
		class identList *id;
		int start_count;
		PyObject *py_start_state = NULL, *py_complex = NULL;
		PyObject *py_seq = NULL, *py_struc = NULL;
		PyObject *py_err = NULL;

		if (alternate_start != NULL)
			py_start_state = alternate_start;
		else
			py_start_state = getListAttrReify(python_settings, start_state);

		if (py_start_state != Py_None)
			// doesn't need reference counting for this size call.
			// the getlistattr call we decref later.
			start_count = PyList_GET_SIZE(py_start_state);
		else
			start_count = 0;

		if (start_count == 0) {	// FD Jun 2018: adding throw if no initial state is set.
			throw std::invalid_argument("Initial state was not set.");
		}

		for (int index = 0; index < start_count; index++) {

			// #ifndef DEBUG_MACROS
			py_complex = PyList_GET_ITEM(py_start_state, index);
			// Borrowed reference, we do NOT decref it at end of loop.

#ifdef DEBUG_MACROS
			printPyError_withLineNumber();
#endif

			sequence = getStringAttrReify(py_complex, sequence, py_seq);
			// new reference

			structure = getStringAttrReify(py_complex, structure, py_struc);
			// new reference
			// Need to check if an error occurred, specifically, it could be an OSError due to sample failing. If so, we need to get the heck out of dodge right now.
			py_err = PyErr_Occurred();
			// py_err is a borrowed reference

			if (py_err != NULL) { // then an error occurred while getting the structure. Test for OSError (sample failure):
				if (PyErr_ExceptionMatches(PyExc_OSError)) {
					fprintf(stderr,
							"MULTISTRAND: Starting Structure could not be retrieved for index %d in your options object's start_state. This is likely due to Boltzmann sampling failing: please check that the program 'sample' exists and points correctly to the NUPACK sample binary. Or try 'print o.start_state[%d].structure' where 'o' is your options object and refer to that error message (if any).\n",
							index, index);
				} else {
					fprintf(stderr, "MULTISTRAND: An unidentified exception occurred while trying to initialize the system.\n");
				}
				return;
			}

			id = getID_list(python_settings, index, alternate_start);

			complex_input myTempComplex = complex_input(sequence, structure, id);

			// StrandComplex does make its own copy of the seq/structure, so we can now decref.
			myComplexes->push_back(myTempComplex);

			Py_DECREF(py_seq);
			Py_DECREF(py_struc);

		}
		Py_DECREF(py_start_state);

		// Update the current seed and store the starting structures
		//   note: only if we actually have a system_options, e.g. no alternate start
		if (alternate_start == NULL && python_settings != NULL) {
			setLongAttr(python_settings, interface_current_seed, current_seed);
		}
		seed = current_seed;

	}

	return;
}

stopComplexes* PSimOptions::getStopComplexes(int) {

	myStopComplexes = getStopComplexList(python_settings, 0);

	return myStopComplexes;

}

void PSimOptions::stopResultError(long seed) {

	if (!statespaceActive) {
		printStatusLine(python_settings, seed, STOPRESULT_ERROR, 0.0, result_type::STR_ERROR.c_str());
	}

}

void PSimOptions::stopResultNan(long seed) {

	if (!statespaceActive) {
		printStatusLine(python_settings, seed, STOPRESULT_NAN, 0.0, result_type::STR_NAN.c_str());
	}

}

void PSimOptions::stopResultNormal(long seed, double time, char* message) {

	if (!statespaceActive) {
		printStatusLine(python_settings, seed, STOPRESULT_NORMAL, time, message);
	}

}

void PSimOptions::stopResultTime(long seed, double time) {

	if (!statespaceActive) {
		printStatusLine(python_settings, seed, STOPRESULT_TIME, time, result_type::STR_TIMEOUT.c_str());
	}

}

void PSimOptions::stopResultFirstStep(long seed, double stopTime, double rate, const char* message) {

	if (!statespaceActive) {
		printStatusLine_First_Bimolecular(python_settings, seed, STOPRESULT_NORMAL, stopTime, rate, message);
	}
}

///// CSIMOPTIONS
CSimOptions::CSimOptions(void) {

	// initializers calling python object -- these can use a super object getter.
	// Not clear at the moment if calling all settings is possible without crashing.
	fixedRandomSeed = true;

	seed = 7777;

	energyOptions = new CEnergyOptions();

	simulation_mode = 16;
	simulation_count = 1000;
	o_time = 0;
	o_interval = 0;
	stop_count = 1;
	stop_options = 1;
	max_sim_time = 0.1;

	debug = false;	// this is the main switch for simOptions debug, for now.

}

PyObject* CSimOptions::getPythonSettings() {

	cout << "getPythonSettings, cannot proceed \n";
	abort();
	return NULL;

}

void CSimOptions::generateComplexes(PyObject *alternate_start, long current_seed) {

	myComplexes = new vector<complex_input>(0); // wipe the pointer to the previous object;

//	// setting default value
//	char* mySeq = "GTTAGACTCGGAGGTGGTAGCAATGGATCAG+CTGATCCATTGCTACCACCTCCGAGTCTAACCATATC+GATATGGTTAGACTCGGAGGTGGTAGCAATG";
//	char* myStructure = ".........................((((((+))))))(((((((((((((((((((((((((((((((+)))))))))))))))))))))))))))))))";
//	identList* myIdentity1 = new identList(1337, "myID-1", NULL);
//	identList* myIdentity2 = new identList(1338, "myID-2", myIdentity1);
//	identList* myIdentity3 = new identList(1339, "myID-3", myIdentity2);
//
//	complex_input myInput = complex_input(mySeq, myStructure, NULL);
//
//	myComplexes->push_back(myInput);

}

stopComplexes* CSimOptions::getStopComplexes(int) {

	cout << "getStopComplexes, cannot proceed \n";
	abort();
	return NULL;
}

void CSimOptions::stopResultError(long seed) {

	cout << "stopResultError, cannot send to python \n";

}

void CSimOptions::stopResultNan(long seed) {

	cout << "stopResultNan, cannot send to python \n";

}

void CSimOptions::stopResultNormal(long seed, double time, char* message) {

	cout << "stopResultNormal, cannot send to python \n";

}

void CSimOptions::stopResultTime(long seed, double time) {

	cout << "stopResultTime, cannot send to python \n";

}

void CSimOptions::stopResultFirstStep(long seed, double stopTime, double rate, const char* message) {

	cout << "stopResultBimolecular, cannot send to python \n";

}

