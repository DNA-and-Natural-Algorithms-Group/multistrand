#ifndef __SEQUTIL_H__
#define __SEQUTIL_H__

#include <stdio.h>
#include <string>
#include <vector>

using std::string; using std::vector;


// Vienna: 0 is invalid, then CG, GC, GU, UG, AU, UA, and Special are 1-7
// MFold/Nupack:  0 is AT, then CG, GC, TA, GT, TG
// 0 is invalid, then A, C, G, U

const int NUM_BASEPAIRS_VIENNA = 8;
const int NUM_BASEPAIRS_NUPACK = 6;
const int NUM_BASES = 5;

const static string basepairString[NUM_BASEPAIRS_NUPACK] = { "A/T", "C/G", "G/C", "T/A", "G/T", "T/G" };

// structure containing information about bases exterior to the complex, IE bases that could pair with other complexes. First incarnation of such.a structure, prolly will change as I work out multiple-complex issues.
// FD: refactoring to stop using exterior_bases, in favor of BaseCounter -- nov 9 2016
struct exterior_bases {
	int A, T, C, G;
};

enum BaseType {

	baseNone, baseA, baseC, baseG, baseT, BASETYPE_SIZE

};

const string baseTypeString[BASETYPE_SIZE] = { "·", "A", "C", "G", "T" };


// structure containing information about bases exterior to the complex, IE bases that could pair with other complexes. First incarnation of such.a structure, prolly will change as I work out multiple-complex issues.
struct BaseCounter {

	vector<int> count = {0, 0, 0, 0, 0}; // use baseType as access

	// Constructor to help interface with existing code
	// This also helps in refactoring
	BaseCounter(exterior_bases* );

	void increment(BaseCounter* other);
	void decrement(BaseCounter* other);
	int multiCount(BaseCounter* other);

	int countFromChar(char c);

//public:
	// Convenience gets
	int A(void);
	int T(void);
	int G(void);
	int C(void);

};

#endif
