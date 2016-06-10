/*
 * utility.cc

 *
 *  Created on: Jun 9, 2016
 *      Author: Frits Dannenberg
 */

#include "utility.h"
#include <string>
#include <cstring>

char* utility::copyToCharArray(string& myString) {

	char* newArray = (char *) new char[myString.length() + 1];
	strcpy(newArray, myString.c_str());

	return newArray;

}


string utility::copyToString(char* myCharArray){

//	char* newArray = new char[myCharArray-> + 1];
////	strcpy(newArray, myString.c_str());

	return string(myCharArray);

}
