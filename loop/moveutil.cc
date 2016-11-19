#include <moveutil.h>
#include <sequtil.h>
#include <assert.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <scomplex.h>

using std::vector;
using std::map;
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

	for (std::pair<HalfContext, BaseCount> value : m.tally) {

		ss << value.first << "  ";
		ss << value.second;
		ss << "\n";

	}

	ss << "Exposed, Intern/Total = " << m.numExposedInternal << " / ";
	ss << m.numExposed << "	\n";

	ss << "\n";

	return ss;

}

void OpenInfo::clear(void) {

	tally.clear();
	numExposedInternal = 0;
	numExposed = 0;

}

// simply store the vector of halfContext onto the list we already have
void OpenInfo::increment(QuartContext left, char base, QuartContext right) {

	HalfContext con = HalfContext(left, right);

	if (tally.count(con) > 0) { // if count > 1 then proceed

		tally.find(con)->second.count[base]++;

	} else {

		BaseCount count = BaseCount();
		count.count[base]++;
		tally.insert(std::pair<HalfContext, BaseCount>(con, count));

	}

}

// constructor assigns the base
LocalContext::LocalContext(char input) {

	base = BaseType(input);

}

JoinCriterea::JoinCriterea() {

	// empty constructor

}

std::ostream& operator<<(std::ostream &ss, JoinCriterea& m) {

	ss << "Types = " << m.types[0] << " " << m.types[1];

	return ss;

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

// constructor assigns the base
HalfContext::HalfContext() {

}

HalfContext::HalfContext(QuartContext in1, QuartContext in2) {

	left = in1;
	right = in2;

}

std::ostream& operator<<(std::ostream &os, HalfContext& m) {

	os << "(" << quartContextString[m.left] << ", ";
	os << quartContextString[m.right] << ") ";

	return os;

}


// c++ expects partial ordering instead of equality for mapping
bool HalfContext::operator<(const HalfContext& other) const {

	return (left < other.left) || (left == other.left && right < other.right);

}

std::ostream& operator<<(std::ostream &os, LocalContext& m) {

	os << "(" << quartContextString[m.half.left] << ", ";
	os << m.base << ", ";
	os << quartContextString[m.half.right] << ") ";

	return os;

}

ContextList::ContextList() {

}

void ContextList::increment(void) {

	//TODO

}

void ContextList::clear() {

	tally.clear();

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

void TransitionList::clear() {

	rateSum = 0.0;
	list.clear();

}
