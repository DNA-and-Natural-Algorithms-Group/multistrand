/*
 * utility.h
 *
 *  Created on: Jun 9, 2016
 *      Author: Frits Dannenberg
 */

#ifndef INCLUDE_UTILITY_H_
#define INCLUDE_UTILITY_H_

#include <vector>
#include <string>
#include <iostream>

using namespace std;


// no implementation file for utilities.
namespace utility{


// Structs
struct complex_input {

	std::string sequence;
	std::string structure;
	identlist* list;

	complex_input() {
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


//void copyToCharArray222(char* myArray, string& myString){
//
////	char* newArray = (char *) new char[myString.length() + 1];
////	strcpy(newArray, myString.c_str());
////
////	myArray = newArray;
//
//};

}


#endif /* INCLUDE_UTILITY_H_ */
