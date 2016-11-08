#ifndef __SEQUTIL_H__
#define __SEQEUTIL_H__

#include <stdio.h>
#include <string>

const int NUM_BASEPAIRS_VIENNA = 8;
// Vienna: 0 is invalid, then CG, GC, GU, UG, AU, UA, and Special are 1-7
const int NUM_BASEPAIRS_NUPACK = 6;
// MFold/Nupack:  0 is AT, then CG, GC, TA, GT, TG
const int NUM_BASES = 5;
// 0 is invalid, then A, C, G, U

enum BaseType {

	baseA, baseC, baseG, baseT, BASETYPE_SIZE

};

const static string basepairString[NUM_BASEPAIRS_NUPACK] = { "A/T", "C/G", "G/C", "T/A", "G/T", "T/G" };

const string baseTypeString[BASETYPE_SIZE + 1] = { "·", "A", "C", "G", "T" };

#endif
