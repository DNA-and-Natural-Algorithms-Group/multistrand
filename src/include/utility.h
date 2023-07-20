/*
Copyright (c) 2017 California Institute of Technology. All rights reserved.
Multistrand nucleic acid kinetic simulator
help@multistrand.org
*/

/*
 *  Created on: Jun 9, 2016
 *      Author: Frits Dannenberg
 */

#ifndef INCLUDE_UTILITY_CC_
#define INCLUDE_UTILITY_CC_

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

#include "optionlists.h"

using namespace std;


namespace utility {

//helper functions

string sequenceToString(char*, int);

char* copyToCharArray(string&);

string copyToString(char*);

string moveType(int);

}

#endif /* INCLUDE_UTILITY_CC_ */
