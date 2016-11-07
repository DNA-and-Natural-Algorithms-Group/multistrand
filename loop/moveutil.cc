

#include <moveutil.h>
#include <string>
#include <iostream>

//using std::string;



std::string quartContextString[HALFCONTEXT_SIZE] = { "endC  ", "strandC", "stackC" };

int moveutil::typeMult(MoveType left, MoveType right) {

	return (int) (valuesPrime[left] * valuesPrime[right]);

}


std::ostream& operator<< (std::ostream &os, halfContext& m) {

	os << "left quart: " << quartContextString[m.left] << ", left right:  " << quartContextString[m.right] << " \n";
	return os;

}
