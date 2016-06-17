/*
 * utility.cc

 *
 *  Created on: Jun 9, 2016
 *      Author: Frits Dannenberg
 */

#include "utility.h"
#include <string>
#include <sstream>
#include "move.h"

char* utility::copyToCharArray(string& myString) {

	char* newArray = (char *) new char[myString.length() + 1];
	strcpy(newArray, myString.c_str());

	return newArray;

}

string utility::copyToString(char* inputCharArray) {

	char* newArray = (char *) new char[strlen(inputCharArray) + 1];
	strcpy(newArray, inputCharArray);

	return string(newArray);

}

string utility::moveType(int type) {

	std::stringstream ss;

	if (type & MOVE_INVALID) {

		ss << "invalid, ";

	}

	if (type & MOVE_CREATE) {

		ss << "create, ";

	}

	if (type & MOVE_DELETE) {

		ss << "delete, ";

	}

	if (type & MOVE_SHIFT) {

		ss << "shift, ";

	}

	if (type & MOVE_1) {

		ss << "Move-1, ";

	}

	if (type & MOVE_2) {

		ss << "Move-2, ";

	}

	if (type & MOVE_3) {

		ss << "Move-3, ";

	}

	string output = ss.str();

	return output;

}

void utility::printDouble(double input) {

	cout << "Printing DOUBLE:" << input << " \n";

}

