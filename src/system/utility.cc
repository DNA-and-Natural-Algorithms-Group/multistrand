/*
Copyright (c) 2017 California Institute of Technology. All rights reserved.
Multistrand nucleic acid kinetic simulator
help@multistrand.org
*/

/*
 *
 *  Created on: Jun 9, 2016
 *      Author: Frits Dannenberg
 */

#include "utility.h"
#include <string>
#include <sstream>
#include "move.h"
#include <energymodel.h>



char* utility::copyToCharArray(string& myString) {

	char* newArray = (char *) new char[myString.length() + 1];
	strcpy(newArray, myString.c_str());

	return newArray;

}

string utility::sequenceToString(char* sequence, int size) {

	// the first and final character is the paired base -- adjust the print for this.

	std::stringstream ss;

	int preBase = (int) sequence[0];
	int postBase = (int) sequence[size+1];

	if( preBase < 0 || preBase > 5){
		cout << "Warning! prebase is outside of range" << endl;
	}

	if( postBase < 0 || postBase > 5){
		cout << "Warning! postbase is outside of range" << endl;
	}



	ss << "";
	ss << baseTypeString[(int) sequence[0]];
	ss << ":";

	for (int i = 1; i < size + 1; i++) {

		ss << baseTypeString[(int) sequence[i]];

	}

	ss << ":";
	ss << baseTypeString[(int) sequence[size + 1]];
	ss << "";

	return ss.str();

}

string utility::copyToString(char* inputCharArray) {

	char* newArray = (char *) new char[strlen(inputCharArray) + 1];
	strcpy(newArray, inputCharArray);

	return string(newArray);

}

string utility::moveType(int type) {

	std::stringstream ss;

	if (type & MOVE_INVALID) {

		ss << "invalid";

	}

	if (type & MOVE_CREATE) {

		ss << "create";

	}

	if (type & MOVE_DELETE) {

		ss << "delete";

	}

	if (type & MOVE_SHIFT) {

		ss << "shift";

	}

	if (type & MOVE_1) {

		ss << "_1, ";

	}

	if (type & MOVE_2) {

		ss << "_2, ";

	}

	if (type & MOVE_3) {

		ss << "_3, ";

	}

	string output = ss.str();

	return output;

}

void utility::printIntegers(int array[], int size) {

	cout << "Integers:";

	for (int i = 0; i < size; i++) {

		cout << array[i];

	}

	cout << "\n";

}

void utility::printDouble(double input) {

	cout << "Printing DOUBLE:" << input << " \n";

}

void utility::printDoubleArray(double input[], int length) {

	cout << "Double array: ";

	for (int i = 0; i < length; i++) {

		cout << input[i] << ",";

	}

}

void utility::printDoubleMatrix(double input[], int length, int depth, int precision) {

	for (int j = 0; j < depth; j++) {

		for (int i = 0; i < length; i++) {

			cout << std::setprecision(precision) << input[j * depth + i] << ",";

		}

		cout << " \n";
		cout.flush();

	}

}

