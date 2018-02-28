/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */

// Contains all move-related utility as a base-class.
#ifndef __MOVEUTIL_H__
#define __MOVEUTIL_H__

#include <string>
#include <vector>
#include <map>
#include <sequtil.h>

using std::vector;
using std::string;
using std::map;

class StrandComplex;
class EnergyModel;

enum MoveType {
	endMove, loopMove, stackMove, stackStackMove, loopEndMove, stackEndMove, stackLoopMove, MOVETYPE_SIZE
};

enum QuartContext {
	endC, strandC, stackC, HALFCONTEXT_SIZE
};

namespace moveutil {

MoveType combineBi(QuartContext&, QuartContext&);

const static double valuesPrime[MOVETYPE_SIZE + 1] = { 3, 5, 7, 11, 13, 17, 19, 999 };
const static string MoveToString[MOVETYPE_SIZE] = { "End", "Loop", "Stack", "StackStack", "LoopEnd", "StackEnd", "StackLoop" };
const static string MoveToString2[MOVETYPE_SIZE] = { "       ", "      ", "     ", "", "   ", "  ", " " };

int getPrimeCode(MoveType, MoveType);
string primeToDesc(int);

QuartContext getContext(char input);

int typeMult(MoveType left, MoveType right);

}

// UTILITY STRUCTS

// structs
struct ExportData {

	int id = 0;
//	char* names = NULL;
	string names;
	string sequence;
	string structure;
	double energy = 0.0;
	double enthalpy = 0.0;

	bool operator==(const ExportData &other) const {
		return (id == other.id && sequence == other.sequence && structure == other.structure);
	}

	void merge(ExportData& other);

};

inline std::ostream & operator<<(std::ostream& str, const ExportData& k) {

	str << k.names + " ";
	str << k.sequence + " ";
	str << k.structure + " ";
	str << std::to_string(k.energy) + " ";
	str << std::to_string(k.enthalpy) + "\n";

	return str;
}

struct StateData {

	vector<ExportData> state;

	bool operator==(const StateData &other) const {

		if (state.size() != other.state.size()) {
			return false;
		}

		auto itA = state.begin();
		auto itB = other.state.begin();

		while (itA != state.end() || itB != other.state.end()) {

			if (!((*itA) == (*itB))) {

				return false;

			}

			itA++;
			itB++;

		}

		return true;

	}

};

struct ExportDataHasher {

	std::size_t operator()(const ExportData k) const {
		using std::size_t;
		using std::hash;

		// Compute individual hash values for first, second and third
		// https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
		size_t res = 17;
		res = res * 31 + hash<string>()(k.sequence);
		res = res * 31 + hash<string>()(k.structure);
		res = res * 31 + hash<int>()(k.id);
		return res;

	}

};

struct HalfContext {

	HalfContext();
	HalfContext(QuartContext, QuartContext);
	friend std::ostream& operator<<(std::ostream&, HalfContext&);
	bool operator==(const HalfContext& other) const;
	bool operator<(const HalfContext&) const;

	QuartContext left = endC;
	QuartContext right = endC;

};

struct JoinCriteria {

	JoinCriteria();
	friend std::ostream& operator<<(std::ostream&, JoinCriteria&);

	StrandComplex* complexes[2] = { NULL, NULL };
	char types[2] = { 0, 0 };
	int index[2] = { 0, 0 };

	// arrhenius rates only
	HalfContext half[2] = { HalfContext(), HalfContext() };
	int arrType = 0; // used for returning the chosen movetype.

};

// This struct contains info computed
// at-time-of-creation for the OpenLoop object.

struct OpenInfo {

public:
	friend std::ostream& operator<<(std::ostream&, OpenInfo&);
	void clear(void);
	void increment(QuartContext, char, QuartContext);
	void increment(HalfContext, BaseCount&);
	void increment(OpenInfo&);

	void decrement(HalfContext, BaseCount&);
	void decrement(OpenInfo&);

	double crossRate(OpenInfo&, EnergyModel&);

	map<HalfContext, BaseCount> tally;

	int numExposedInternal = 0;
	int numExposed = 0;

	bool upToDate = false;

};

#endif
