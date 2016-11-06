
// Contains all move-related utility as a base-class.
#ifndef __MOVEUTIL_H__
#define __MOVEUTIL_H__

#include <string>

using namespace std;

enum MoveType {
	endMove, loopMove, stackMove, stackStackMove, loopEndMove, stackEndMove, stackLoopMove, MOVETYPE_SIZE
};


enum quartContext {
	endC, strandC, stackC, QUARTCONTEXT_SIZE
};

struct halfContext {

	quartContext left;
	quartContext right;

};

namespace moveutil {

//double valuesPrime[MOVETYPE_SIZE];
//string MoveToString[MOVETYPE_SIZE];
//string MoveToString2[MOVETYPE_SIZE];

const static double valuesPrime[MOVETYPE_SIZE] = { 3, 5, 7, 11, 13, 17, 19 };
const static string MoveToString[MOVETYPE_SIZE] = { "End", "Loop", "Stack", "StackStack", "LoopEnd", "StackEnd", "StackLoop" };
const static string MoveToString2[MOVETYPE_SIZE] = { "      ", "     ", "   ", "", "  ", " ", " " };

int typeMult(MoveType left, MoveType right);

}

#endif
