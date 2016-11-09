#ifndef __SEQUTIL_H__
#define __SEQUTIL_H__

#include <stdio.h>
#include <string>

using std::string;

// Vienna: 0 is invalid, then CG, GC, GU, UG, AU, UA, and Special are 1-7
// MFold/Nupack:  0 is AT, then CG, GC, TA, GT, TG
// 0 is invalid, then A, C, G, U

const int NUM_BASEPAIRS_VIENNA = 8;
const int NUM_BASEPAIRS_NUPACK = 6;
const int NUM_BASES = 5;

enum BaseType {

	baseNone, baseA, baseC, baseG, baseT, BASETYPE_SIZE

};

const static string basepairString[NUM_BASEPAIRS_NUPACK] = { "A/T", "C/G", "G/C", "T/A", "G/T", "T/G" };

const string baseTypeString[BASETYPE_SIZE] = { "·", "A", "C", "G", "T" };

#endif
