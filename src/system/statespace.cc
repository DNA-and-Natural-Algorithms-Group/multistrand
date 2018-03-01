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

void Builder::addState(ExportData& data, int arrType) {

	protoSpace.insert(data);

	// also record the transition itself.
	if (lastState.complex_count > 0) {

		ExportTransition trans;
		trans.state1 = lastState;
		trans.state2 = data;
		trans.type = arrType;

		protoTransitions.insert(std::move(trans));

	}

	lastState = data;

}

// Write the statespaces to file and reset the hashmaps.
//
void Builder::writeToFile(void) {

	// states

	unordered_set<ExportData>::iterator itr;
	string filename = "protospace-" + string(to_string(simOptions->getSeed())) + ".txt";

	std::ofstream myfile;
	myfile.open(filename);

	for (auto element : protoSpace) {

		myfile << element;
	}

	myfile.close();

	// now clear the map
	protoSpace.clear();

	// transitions

	filename = "prototransitions-" + string(to_string(simOptions->getSeed())) + ".txt";
	myfile.open(filename);

	for (auto element : protoTransitions) {

		myfile << element;

	}

	myfile.close();

}

