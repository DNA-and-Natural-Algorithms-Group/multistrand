/*
 * simOptions.cc
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

	ss << "myComplexes = { ";

	for (int i = 0; i < myComplexes->size(); i++) {

		ss << "{ " << myComplexes->at(i).sequence << ", " << myComplexes->at(i).structure << " }";

	}

	ss << "} \n";

	ss << " myStopComplexes = { ";

	// linked list iterator
	stopComplexes* myStopComplex = myStopComplexes; // copying pointer so we can iterate.

	ss << "} \n";

	string output = ss.str();

	output += energyOptions->toString();

	return output;

}

bool SimOptions::useFixedRandomSeed() {

	return fixedRandomSeed;

}

long SimOptions::getInitialSeed() {

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
		char *sequence, *structure;
		class identList *id;
		int start_count;
		PyObject *py_start_state = NULL, *py_complex = NULL;
		PyObject *py_seq = NULL, *py_struc = NULL;
		PyObject *py_err = NULL;

		if (alternate_start != NULL)
			py_start_state = alternate_start;
		else
			py_start_state = getListAttr(python_settings, start_state);

		start_count = PyList_GET_SIZE(py_start_state);
		// doesn't need reference counting for this size call.
		// the getlistattr call we decref later.

		for (int index = 0; index < start_count; index++) {

			// #ifndef DEBUG_MACROS
			py_complex = PyList_GET_ITEM(py_start_state, index);
			// Borrowed reference, we do NOT decref it at end of loop.

#ifdef DEBUG_MACROS
			printPyError_withLineNumber();
#endif

			sequence = getStringAttr(py_complex, sequence, py_seq);
			// new reference

			structure = getStringAttr(py_complex, structure, py_struc);
			// new reference
			// Need to check if an error occurred, specifically, it could be an IOError due to sample failing. If so, we need to get the heck out of dodge right now.
			py_err = PyErr_Occurred();
			// py_err is a borrowed reference

			if (py_err != NULL) { // then an error occurred while getting the structure. Test for IOError (sample failure):
				if (PyErr_ExceptionMatches(PyExc_IOError)) {
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

	printStatusLine(python_settings, seed, STOPRESULT_ERROR, 0.0, NULL);
	return;

}

void PSimOptions::stopResultNan(long seed) {

	printStatusLine(python_settings, seed, STOPRESULT_NAN, 0.0, NULL);
	return;

}

void PSimOptions::stopResultNormal(long seed, double time, char* message) {

	printStatusLine(python_settings, seed, STOPRESULT_NORMAL, time, message);
	return;

}

void PSimOptions::stopResultTime(long seed, double time) {

	printStatusLine(python_settings, seed, STOPRESULT_TIME, time, NULL);
	return;

}

void PSimOptions::stopResultBimolecular(string type, long seed, double stopTime, double rate, char* message) {

	if (type.compare("Reverse")) {

		printStatusLine_First_Bimolecular(python_settings, seed, STOPRESULT_REVERSE, stopTime, rate, message);

	} else if (type.compare("Forward")) {

		printStatusLine_First_Bimolecular(python_settings, seed, STOPRESULT_FORWARD, stopTime, rate, message);

	} else if (type.compare("FTime")) {

		printStatusLine_First_Bimolecular(python_settings, seed, STOPRESULT_FTIME, stopTime, rate, NULL);

	} else if (type.compare("NoMoves")) {
		printStatusLine_First_Bimolecular(python_settings, seed, STOPRESULT_NOMOVES, stopTime, rate, NULL);

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

	// setting default value
	char* mySeq = "GTTAGACTCGGAGGTGGTAGCAATGGATCAG+CTGATCCATTGCTACCACCTCCGAGTCTAACCATATC+GATATGGTTAGACTCGGAGGTGGTAGCAATG";
	char* myStructure = ".........................((((((+))))))(((((((((((((((((((((((((((((((+)))))))))))))))))))))))))))))))";
	identList* myIdentity1 = new identList(1337, "myID-1", NULL);
	identList* myIdentity2 = new identList(1338, "myID-2", myIdentity1);
	identList* myIdentity3 = new identList(1339, "myID-3", myIdentity2);

	complex_input myInput = complex_input(mySeq, myStructure, NULL);

	myComplexes->push_back(myInput);

	return;
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

void CSimOptions::stopResultBimolecular(string type, long seed, double stopTime, double rate, char* message) {

	cout << "stopResultBimolecular, cannot send to python \n";

}

