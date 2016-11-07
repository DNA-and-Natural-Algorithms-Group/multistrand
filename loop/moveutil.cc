

#include <moveutil.h>
#include <string>

using namespace std;

//double moveutil::valuesPrime[MOVETYPE_SIZE] = { 3, 5, 7, 11, 13, 17, 19 };
//string moveutil::MoveToString[MOVETYPE_SIZE] = { "End", "Loop", "Stack", "StackStack", "LoopEnd", "StackEnd", "StackLoop" };
//string moveutil::MoveToString2[MOVETYPE_SIZE] = { "      ", "     ", "   ", "", "  ", " ", " " };
//

int moveutil::typeMult(MoveType left, MoveType right) {

	return (int) (valuesPrime[left] * valuesPrime[right]);

}
