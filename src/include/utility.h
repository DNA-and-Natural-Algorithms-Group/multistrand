/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef INCLUDE_UTILITY_CC_
#define INCLUDE_UTILITY_CC_

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

#include "optionlists.h"
#include "basetype.h"

using namespace std;


namespace utility {

//helper functions

string sequenceToString(BaseType*, int);

char* copyToCharArray(string&);

string copyToString(char*);

string moveType(int);

}

#endif /* INCLUDE_UTILITY_CC_ */
