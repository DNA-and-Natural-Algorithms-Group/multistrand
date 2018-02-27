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

#include <simoptions.h>
#include <scomplexlist.h>
//#include <>

class Builder {
public:

	Builder(SimOptions& options);

	void addState(SComplexList* state, int arrType);
	void writeToFile(void);

private:

	SimOptions* simOptions;



};
