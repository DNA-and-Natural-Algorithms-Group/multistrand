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

const bool debugTraces = true;


// Structs
struct complex_input {

	std::string sequence;
	std::string structure;
	identList* list;

	complex_input() { // empty constuctor
		sequence = "default";
		structure = "default";
		list = NULL;
	}

	complex_input(const char* string1, const char* string2, identList* list1) {

		char *tempseq = (char *) new char[strlen(string1) + 1];
		char *tempstruct = (char *) new char[strlen(string2) + 1];
		strcpy(tempseq, string1);
		strcpy(tempstruct, string2);

		sequence = string(tempseq);
		structure = string(tempstruct);
		list = list1;

	}
};


//helper functions

string sequenceToString(char*, int);

char* copyToCharArray(string&);

string copyToString(char*);

string moveType(int);



void printIntegers(int[], int);
void printDouble(double);
void printDoubleArray(double[], int);
void printDoubleMatrix(double[], int, int, int);


}





#endif /* INCLUDE_UTILITY_CC_ */
