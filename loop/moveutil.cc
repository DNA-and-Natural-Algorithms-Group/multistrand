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

QuartContext moveutil::getContext(char input) {

	if (input > 0) { // there is a stack on the exteriour

		return stackC;

	} else {

		return endC;

	}

}

std::ostream& operator<<(std::ostream &ss, OpenInfo& m) {

	// prints vector<halfContext>
	for (vector<HalfContext> & value : m.context) {
		for (HalfContext& context : value) {
			ss << context << ", ";
		}
		ss << " \n";	// prints empty line for empty interior sequences
	}

	ss << "#ExposedInternalNucl= " << m.numExposedInternal << "\n";
	ss << "#ExposedNucl=         " << m.numExposed << "\n";

	ss << "\n";

	for (int i : { 0, 1, 2, 3, 4 }) {

		ss << baseTypeString[i] << ": " << m.exposedInternalNucl[i] << " ";

	}

	return ss;

}

void OpenInfo::clear(void) {

	context.clear();
	exposedInternalNucl = {0,0,0,0,0};
	numExposedInternal = 0;
	numExposed = 0;

}

// simply store the vector of halfContext onto the list we already have
void OpenInfo::push(vector<HalfContext>& input) {

	context.push_back(input);

}

// constructor assigns the base
HalfContext::HalfContext(char input) {

	base = BaseType(input);

}

std::ostream& operator<<(std::ostream &os, HalfContext& m) {

	os << "(" << quartContextString[m.left] << ", " << quartContextString[m.right] << ") ";

	return os;

}

