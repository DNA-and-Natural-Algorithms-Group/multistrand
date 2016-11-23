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
bool isPair(BaseType, BaseType);

const static double valuesPrime[MOVETYPE_SIZE] = { 3, 5, 7, 11, 13, 17, 19 };
const static string MoveToString[MOVETYPE_SIZE] = { "End", "Loop", "Stack", "StackStack", "LoopEnd", "StackEnd", "StackLoop" };
const static string MoveToString2[MOVETYPE_SIZE] = { "       ", "      ", "     ", "", "   ", "  ", " " };

QuartContext getContext(char input);

int typeMult(MoveType left, MoveType right);

}

// UTILITY STRUCTS

struct JoinCriteria {

	JoinCriteria();
	friend std::ostream& operator<<(std::ostream&, JoinCriteria&);

	StrandComplex* picked[2] = { NULL, NULL };
	char types[2] = { 0, 0 };
	int index[2] = { 0, 0 };

};

struct HalfContext {

	HalfContext();
	HalfContext(QuartContext, QuartContext);
	friend std::ostream& operator<<(std::ostream&, HalfContext&);
	bool operator<(const HalfContext&) const;

	QuartContext left = endC;
	QuartContext right = endC;

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

//class ContextList {
//
//	ContextList();
//	void increment(OpenInfo&);
//	void clear(void);
//
//	// contains a BaseCount for each possible half-context.
//	map<HalfContext, BaseCount> tally;
//
//};

// UTILITY CLASSES

// general-purpose class that describes transitions and a
// general-purpose transition container class TransitionList
// FD: these are very similar to the MoveContainer, MoveList, and Move classes,
// FD: except that for this class, we do not include any program logic.
// FD: They only store rates and pointers to actual data.
// FD: Nov 18 2016: these will be depreciated soon.
//class Transition {
//
//public:
//
//	Transition(double, char);
//	friend std::ostream& operator<<(std::ostream&, Transition&);
//
//private:
//	double rate = 0.0;
//	char pairType;	 // a 2-set of which nucleotides are binding in this move
//
//};
//
//class TransitionList {
//
//public:
//	void push(double, char);
//	void clear(void);
//
//	friend std::ostream& operator<<(std::ostream&, TransitionList&);
//
//	double rateSum = 0.0;
//
//private:
//	vector<Transition> list;
//
//};

#endif
