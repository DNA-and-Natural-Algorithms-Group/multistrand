/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */

/*
 *  Created on: Feb 23, 2018
 *      Author: Frits Dannenberg
 *
 *      This class collects every visited state, and transitions between visited states.
 *      After the simulation is done, the set is exported to a text file.
 *
 */

#include <statespace.h>
#include <simoptions.h>

#include <iostream>
#include <fstream>

Builder::Builder(void) {

}

Builder::Builder(SimOptions* options) {

	simOptions = options;

}

// Put the statespace in memory
// Note: this copies the contents of the input.
// Could use shared pointer to prevent this,
// but statespace building will only be used during inference anyhow.
// also note the copy only occurs if the state / transition is not found already
// this is a relatively rare event for a typicall simulation.

void Builder::addState(ExportData& data, const double arrType) {

	protoSpace.insert(data);

	// also record the transition itself.
	if (lastState.complex_count > 0) {

		ExportTransition trans;
		trans.state1 = lastState;
		trans.state2 = data;
		trans.type = arrType;

		protoTransitions.insert(std::move(trans));

	} else { // set the initial state

		auto element = protoInitialStates.find(data);

		if (element == protoInitialStates.end()) {

			ExportInitial newEntry = ExportInitial();
			newEntry.join_rate = arrType; // overloading arrType to be join rate
			newEntry.observation_count++;

			protoInitialStates[data] = std::move(newEntry);

		} else {

			(*element).second.observation_count++;

		}

	}

	lastState = std::move(data);

}

// export the final state to the appropriate map.
void Builder::stopResultNormal(double endtime, string tag) {

	if (lastState.complex_count > 0) {

		auto element = protoFinalStates.find(lastState);

		if (element == protoFinalStates.end()) {

			ExportFinal newEntry = ExportFinal();
			newEntry.tag = tag;
			newEntry.observation_count++;

			protoFinalStates[std::move(lastState)] = std::move(newEntry);

		} else {

			(*element).second.observation_count++;

		}

	}

	// now wipe the last state because we have moved the contents
	lastState = ExportData();

}

string Builder::filename(string input) {

	return string(to_string(simOptions->getSeed())) + "/" + input + ".txt";

}

// Write the statespaces to file and reset the hashmaps.
//
void Builder::writeToFile(void) {

	// create dir
	// system call to create a directory  TODO fix this to use experimental/filesystem
	system((string("mkdir -p ") + to_string(simOptions->getSeed())).c_str());

	// states
	unordered_set<ExportData>::iterator itr;

	std::ofstream myfile;
	myfile.open(filename("protospace"));

	for (auto element : protoSpace) {

		myfile << element;
	}

	myfile.close();

	// transitions
	myfile.open(filename("prototransitions"));

	for (auto element : protoTransitions) {

		myfile << element;

	}

	myfile.close();

	// init
	myfile.open(filename("protoinitialstates"));

	for (auto element : protoInitialStates) {

		myfile << element.first;
		myfile << element.second;

	}

	myfile.close();

	// final states
	myfile.open(filename("protofinalstates"));

	for (auto element : protoFinalStates) {

		myfile << element.first;
		myfile << element.second;

	}

	myfile.close();

	// now clear the maps
	protoSpace.clear();
	protoTransitions.clear();
	protoFinalStates.clear();
	protoInitialStates.clear();

}

