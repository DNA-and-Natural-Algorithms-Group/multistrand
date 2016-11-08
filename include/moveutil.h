// Contains all move-related utility as a base-class.
#ifndef __MOVEUTIL_H__
#define __MOVEUTIL_H__

#include <string>
#include <vector>

using std::vector;
using std::string;


enum MoveType {
	endMove, loopMove, stackMove, stackStackMove, loopEndMove, stackEndMove, stackLoopMove, MOVETYPE_SIZE
};

enum quartContext {
	endC, strandC, stackC, HALFCONTEXT_SIZE
};

struct halfContext {

public:
	friend std::ostream& operator<<(std::ostream&, halfContext&);

	quartContext left = endC;
	quartContext right = endC;

};

// This struct contains info computed
// at-time-of-creation for the OpenLoop object.

struct openInfo {

public:
	friend std::ostream& operator<<(std::ostream&, openInfo&);

	vector<vector<halfContext>> context;
	int exposedInternalNucl[5];
	int numExposedInternal;

};

namespace moveutil {

const static double valuesPrime[MOVETYPE_SIZE] = { 3, 5, 7, 11, 13, 17, 19 };
const static string MoveToString[MOVETYPE_SIZE] = { "End", "Loop", "Stack", "StackStack", "LoopEnd", "StackEnd", "StackLoop" };
const static string MoveToString2[MOVETYPE_SIZE] = { "      ", "     ", "   ", "", "  ", " ", " " };

int typeMult(MoveType left, MoveType right);

}

#endif
