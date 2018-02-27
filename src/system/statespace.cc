/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */

/*
 *  Created on: Feb 23, 2018
 *      Author: Frits Dannenberg
 *
 *      This class collects every visisted state, and transitions between visited states.
 *      After the simulation is done, the set is exported to a text file.
 *
 */

#include <statespace.h>

Builder::Builder(SimOptions& options) {

	simOptions = &options;

}

// Write the statespaces to file and reset the hashmaps.
//
void addState(SComplexList* state, int arrType) {

}

// Write the statespaces to file and reset the hashmaps.
//
void Builder::writeToFile(void) {

}

