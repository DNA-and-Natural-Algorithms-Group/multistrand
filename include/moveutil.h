// Contains all move-related utility as a base-class.
#ifndef __MOVEUTIL_H__
#define __MOVEUTIL_H__

#include <string>
#include <vector>
#include <sequtil.h>

using std::vector;
using std::string;

enum MoveType {
	endMove, loopMove, stackMove, stackStackMove, loopEndMove, stackEndMove, stackLoopMove, MOVETYPE_SIZE
};

enum QuartContext {
	endC, strandC, stackC, HALFCONTEXT_SIZE
};

namespace moveutil {

MoveType combine(QuartContext&, QuartContext&);
bool isPair(BaseType, BaseType);

const static double valuesPrime[MOVETYPE_SIZE] = { 3, 5, 7, 11, 13, 17, 19 };
const static string MoveToString[MOVETYPE_SIZE] = { "End", "Loop", "Stack", "StackStack", "LoopEnd", "StackEnd", "StackLoop" };
const static string MoveToString2[MOVETYPE_SIZE] = { "      ", "     ", "   ", "", "  ", " ", " " };

QuartContext getContext(char input);

int typeMult(MoveType left, MoveType right);

}

// general-purpose class that describes transitions and a
// general-purpose transition container class TransitionList
// FD: these are very similar to the MoveContainer, MoveList, and Move classes,
// FD: except that for this class, we do not include any program logic.
// FD: They only store rates and pointers to actual data.
class Transition {

public:

	Transition(double, char);
	friend std::ostream& operator<<(std::ostream&, Transition&);

private:
	double rate = 0.0;
	char pairType;	 // a 2-set of which nucleotides are binding in this move

};

class TransitionList {

public:
//	void push(double, BaseType, BaseType);
	void push(double, char);
	void clear(void);

	friend std::ostream& operator<<(std::ostream&, TransitionList&);

	double rateSum = 0.0;

private:
	vector<Transition> list;

};

struct HalfContext {

	HalfContext(char base);

//public:

	friend std::ostream& operator<<(std::ostream&, HalfContext&);

	BaseType base;
	QuartContext left = endC;
	QuartContext right = endC;

};

// This struct contains info computed
// at-time-of-creation for the OpenLoop object.

struct OpenInfo {

public:
	friend std::ostream& operator<<(std::ostream&, OpenInfo&);
	void clear(void);
	void push(vector<HalfContext>&);

	vector<vector<HalfContext>> context;
	vector<int> exposedInternalNucl = { 0, 0, 0, 0, 0 };
	int numExposedInternal;
	int numExposed;

};

#endif
