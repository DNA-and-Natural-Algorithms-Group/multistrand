/*
 * utility.h
 *
 *  Created on: Jun 9, 2016
 *      Author: Frits Dannenberg
 */

#ifndef INCLUDE_UTILITY_CC_
#define INCLUDE_UTILITY_CC_

#include <vector>
#include <string>
#include <iostream>
#include "optionlists.h"

using namespace std;

// no implementation file for utilities.
namespace utility {

// Structs
struct complex_input {

	std::string sequence;
	std::string structure;
	identlist* list;

	complex_input() { // empty constuctor
		sequence = "default";
		structure = "default";
		list = NULL;
	}

	complex_input(char* string1, char* string2, identlist* list1) {

		char *tempseq = (char *) new char[strlen(string1) + 1];
		char *tempstruct = (char *) new char[strlen(string2) + 1];
		strcpy(tempseq, string1);
		strcpy(tempstruct, string2);

		sequence = string(tempseq);
		structure = string(tempstruct);
		list = list1;

	}
};

char* copyToCharArray(string&);

string copyToString(char*);

}

#endif /* INCLUDE_UTILITY_CC_ */
