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
#ifndef __STATESPACE_H__
#define __STATESPACE_H__

//#include <simoptions.h>
#include <scomplexlist.h>

#include <unordered_map>
#include <unordered_set>

using std::unordered_map;
using std::unordered_set;


class SimOptions;


struct ExportTransition {

	ExportData state1, state2;

};



/*
 *		key: states. Value: energy
 * 		self.protoSpace = dict()
 *
 *		key: transitions. Value: ArrheniusType (negative if it is a bimolecular transition)
 *		self.protoTransitions = dict()
 *
 *		key: states. Value: a initCountFlux object that tells how many times the state has been the initial state and the join flux (rate)
 *		self.protoInitialStates = dict()
 *
 *		key: states: Value: the result of this final state can be SUCCES or FAILURE
 *		self.protoFinalStates = dict()
 */
class Builder {
public:

	Builder();
	Builder(SimOptions* options);

	void addState(ExportData);
	void writeToFile(void);

private:

	ExportData lastState;
	SimOptions* simOptions = NULL;

	unordered_set<ExportData,  ExportDataHasher> protoSpace;
//	unordered_map<transition, double> protoTransitions;


};


#endif
