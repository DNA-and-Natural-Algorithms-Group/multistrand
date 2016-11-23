// nothing here FD nov 8 2016

#include <sequtil.h>
#include <stdio.h>
#include <iostream>

// base constructor
BaseCount::BaseCount() {

	// empty

}

// compatibility constructor
BaseCount::BaseCount(int* input) {

	count[1] = input[1];
	count[2] = input[2];
	count[3] = input[3];
	count[4] = input[4];

}

std::ostream& operator<<(std::ostream& ss, BaseCount& input) {

//	ss << "A/C/G/T  ";

	for (int i : { baseA, baseC, baseG, baseT }) {

		ss << input.count[i];

		if (i != baseT) {

			ss << "/";

		}

	}

//	ss << "    \n";

	return ss;

}

void BaseCount::clear(void) {

	count = {0,0,0,0,0};

}

void BaseCount::increment(BaseCount& other) {

	for (int i : { baseA, baseC, baseG, baseT }) {
		count[i] += other.count[i];
	}
}

void BaseCount::decrement(BaseCount& other) {

	for (int i : { baseA, baseC, baseG, baseT }) {
		count[i] -= other.count[i];
	}
}

int BaseCount::multiCount(BaseCount& other) {

	int output = count[baseA] * other.count[baseT];
	output += count[baseC] * other.count[baseG];
	output += count[baseG] * other.count[baseC];
	output += count[baseT] * other.count[baseA];

	return output;

}

int BaseCount::countFromChar(char c) {

	BaseType type = BaseType(c);

	return count[type];

}

int BaseCount::A(void) {

	return count[baseA];

}

int BaseCount::T(void) {

	return count[baseT];

}

int BaseCount::G(void) {

	return count[baseG];

}

int BaseCount::C(void) {

	return count[baseC];

}
