#include <moveutil.h>
#include <sequtil.h>
#include <string>
#include <vector>
#include <iostream>

using std::vector;

std::string quartContextString[HALFCONTEXT_SIZE] = { "end", "loop", "stack" };

int moveutil::typeMult(MoveType left, MoveType right) {

	return (int) (valuesPrime[left] * valuesPrime[right]);

}

std::ostream& operator<<(std::ostream &ss, openInfo& m) {

	// prints vector<halfContext>
	for (vector<halfContext> & value : m.context) {
		for (halfContext& context : value) {
			ss << context << ", ";
		}
		ss << " \n";
	}

	ss << "#ExposedInternalNucl= " << m.numExposedInternal;

	ss << "\n";

	for (int i : { 0, 1, 2, 3, 4, 5 }) {

		ss << baseTypeString[ i ] << ": " << m.exposedInternalNucl[i] << " ";

	}


	return ss;

}

std::ostream& operator<<(std::ostream &os, halfContext& m) {

	os << "(" << quartContextString[m.left] << ", " << quartContextString[m.right] << ") ";

	return os;

}

