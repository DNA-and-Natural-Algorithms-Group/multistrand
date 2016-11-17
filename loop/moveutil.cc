#include <moveutil.h>
#include <sequtil.h>
#include <assert.h>
#include <string>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;

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

MoveType moveutil::combine(QuartContext& one, QuartContext& two) {

	// c++ doesn't do double variable switch

	if (one == endC) {

		if (two == endC) {

			return endMove;

		}

		if (strandC) {

			return loopEndMove;

		}

		if (two == stackC) {

			return stackEndMove;

		}

	}

	if (one == strandC) {

		if (two == endC) {

			return loopEndMove;

		}

		if (strandC) {

			return loopMove;

		}

		if (two == stackC) {

			return stackLoopMove;

		}

	}

	if (one == stackC) {

		if (two == endC) {

			return stackEndMove;

		}

		if (strandC) {

			return stackLoopMove;

		}

		if (two == stackC) {

			return stackStackMove;

		}

	}

	cout << "Failure to recongize quarter context in external nucleotides";
	assert(0);

	return MOVETYPE_SIZE;

}

bool moveutil::isPair(BaseType one, BaseType two) {

	return (one + two == 5);

}

std::ostream& operator<<(std::ostream &os, HalfContext& m) {

	os << "(" << quartContextString[m.left] << ", " << quartContextString[m.right] << ") ";

	return os;

}

Transition::Transition(double rateIn, char pairTypeIn) {

	rate = rateIn;
	pairType = pairTypeIn;

}

std::ostream& operator<<(std::ostream &ss, Transition& m) {

	ss << "rate: " << m.rate << "   ";
	ss << "Bond: " << basepairString[m.pairType];

	ss << "\n";

	return ss;

}

std::ostream& operator<<(std::ostream &ss, TransitionList& m) {

	ss << "Printing transitionList: \n";

	for (Transition trans : m.list) {

		ss << trans;

	}

	return ss;

}

void TransitionList::push(double rate, char nucleotides) {

	Transition trans = Transition(rate, nucleotides);
	list.push_back(trans);

}

//void TransitionList::push(double rate, BaseType left, BaseType right) {
//
//
//	Transition trans = Transition(rate, (char) left);
//	list.push_back(trans);
//
//}

void TransitionList::clear() {

	rateSum = 0.0;
	list.clear();

}
