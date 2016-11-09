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

struct HalfContext {

	HalfContext(char base);

public:
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
	BaseCounter exposedInternalNucl = BaseCounter();
	int numExposedInternal;
	int numExposed;

};

namespace moveutil {

const static double valuesPrime[MOVETYPE_SIZE] = { 3, 5, 7, 11, 13, 17, 19 };
const static string MoveToString[MOVETYPE_SIZE] = { "End", "Loop", "Stack", "StackStack", "LoopEnd", "StackEnd", "StackLoop" };
const static string MoveToString2[MOVETYPE_SIZE] = { "      ", "     ", "   ", "", "  ", " ", " " };


QuartContext getContext(char input);

int typeMult(MoveType left, MoveType right);

}

#endif
