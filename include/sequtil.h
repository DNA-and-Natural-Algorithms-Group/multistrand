#ifndef __SEQUTIL_H__
#define __SEQUTIL_H__

#include <stdio.h>
#include <string>
#include <vector>

using std::string;
using std::vector;

const bool OLD_LOOP_ADJUSTMENT_NASIM = false;


// Vienna: 0 is invalid, then CG, GC, GU, UG, AU, UA, and Special are 1-7
// MFold/Nupack:  0 is AT, then CG, GC, TA, GT, TG
// 0 is invalid, then A, C, G, U

const int NUM_BASEPAIRS_VIENNA = 8;
const int NUM_BASEPAIRS_NUPACK = 6;
const int NUM_BASES = 5;

// FD: I keep having to remove 1 from the pairtype internal representation.
// FD: In this regard, multistrand follows Vienna notation.
// FD: So I added an offset to the tostring array.
const static string basepairString[NUM_BASEPAIRS_NUPACK + 1] = { "VOID", "A/T", "C/G", "G/C", "T/A", "G/T", "T/G" };

enum BaseType {

	baseNone, baseA, baseC, baseG, baseT, BASETYPE_SIZE

};

const string baseTypeString[BASETYPE_SIZE] = { "·", "A", "C", "G", "T" };

// structure containing information about bases exterior to the complex, IE bases that could pair with other complexes. First incarnation of such.a structure, prolly will change as I work out multiple-complex issues.
struct BaseCount {

	// Constructor to help interface with existing code
	// This also helps in refactoring
	// also declare base constructor because auto constructor has been voided
	BaseCount();
	BaseCount(int*);	 // second friendly constructor for refactoring

	void clear(void);

	friend std::ostream& operator<<(std::ostream&, BaseCount&);

	void increment(BaseCount& other);
	void decrement(BaseCount& other);
	int multiCount(BaseCount& other);

	int countFromChar(char c);

//public:
	// Convenience gets
	int A(void);
	int T(void);
	int G(void);
	int C(void);

	// the actual data structure; rest is convienience
	vector<int> count = { 0, 0, 0, 0, 0 }; // use baseType as access
};

#endif
