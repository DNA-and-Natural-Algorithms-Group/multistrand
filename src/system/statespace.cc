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

Builder::Builder(void) {

}

Builder::Builder(SimOptions* options) {

	simOptions = options;

}

// Put the statespace in memory
//

void Builder::addState(ExportData data) {

	protoSpace.insert(data);

}

// Write the statespaces to file and reset the hashmaps.
//
void Builder::writeToFile(void) {

	string filename = string(to_string(simOptions->getSeed()));

}

