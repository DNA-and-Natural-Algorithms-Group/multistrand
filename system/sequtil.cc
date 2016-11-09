// nothing here FD nov 8 2016

#include <sequtil.h>
#include <stdio.h>
#include <iostream>

//using std::cout;

// base constructor
BaseCounter::BaseCounter() {

	// empty

}

// compatibility constructor
BaseCounter::BaseCounter(exterior_bases* input) {

	count[1] = input->A;
	count[2] = input->C;
	count[3] = input->G;
	count[4] = input->T;

}

void BaseCounter::increment(BaseCounter* other) {

	std::cout << "Attempting to increment";
	std::cout.flush();

	for (int i : { baseA, baseC, baseG, baseT }) {
		count[i] += other->count[i];
	}
}

void BaseCounter::decrement(BaseCounter* other) {

	std::cout << "Attempting to increment";
	std::cout.flush();

	for (int i : { baseA, baseC, baseG, baseT }) {
		count[i] -= other->count[i];
	}
}

int BaseCounter::multiCount(BaseCounter* other) {

	int output = 0;

	output += count[baseA] * other->count[baseT];
	output += count[baseC] * other->count[baseG];
	output += count[baseG] * other->count[baseC];
	output += count[baseT] * other->count[baseA];

	return output;

}

int BaseCounter::countFromChar(char c) {

	BaseType type = BaseType(c);

	return count[type];

}

int BaseCounter::A(void) {

	return count[baseA];

}

int BaseCounter::T(void) {

	return count[baseT];

}

int BaseCounter::G(void) {

	return count[baseG];

}

int BaseCounter::C(void) {

	return count[baseC];

}
