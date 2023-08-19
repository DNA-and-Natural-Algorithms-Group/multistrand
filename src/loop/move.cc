/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

// TEST
#include <stdio.h>
#include <assert.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "move.h"
#include "loop.h"
#include "utility.h"
#include "energyoptions.h"
#include "energymodel.h"

using std::string;

RateEnv::RateEnv(void) {

	rate = 0;
	arrType = -99;

}

RateEnv::RateEnv(double mRate, EnergyModel* eModel, MoveType left, MoveType right) {

	rate = eModel->applyPrefactors(mRate, left, right);
	arrType = moveutil::typeMult(left, right);

}

std::ostream& operator<<(std::ostream& ss, RateEnv& env) {

	ss << env.arrType << " ";

	if (env.arrType < 100) {

		ss << " ";

	}

	ss << moveutil::primeToDesc(int(env.arrType));
	ss << "  " << env.rate << "   ";

	return ss;

}

string RateEnv::toString(bool useArr) {

	stringstream ss;

	if (useArr) {

		ss << *this;

	} else {

		ss << "  " << rate << "   ";

	}
	return ss.str();

}

Move::Move(void) {

	type = 0;
	rate = RateEnv();
	affected[0] = affected[1] = NULL;

}

Move::Move(int mtype, RateEnv mrate, Loop *affected_1, int index1, int index2) {

	type = mtype;
	rate = mrate;
	affected[0] = affected_1;
	affected[1] = NULL;
	index[0] = index1;
	index[1] = index2;
	index[2] = -1;
	index[3] = -1;

}

Move::Move(int mtype, RateEnv mrate, Loop *affected_1, int index1, int index2, int index3) {
	type = mtype;
	rate = mrate;
	affected[0] = affected_1;
	affected[1] = NULL;
	index[0] = index1;
	index[1] = index2;
	index[2] = index3;
	index[3] = -1;

}

Move::Move(int mtype, RateEnv mrate, Loop *affected_1, Loop *affected_2, int index1, int index2) {
	type = mtype;
	rate = mrate;
	affected[0] = affected_1;
	affected[1] = affected_2;
	index[0] = index1;
	index[1] = index2;
	index[2] = -1;
	index[3] = -1;

}

Move::Move(int mtype, RateEnv mrate, Loop *affected_1, Loop *affected_2, int index1) {
	type = mtype;
	rate = mrate;
	affected[0] = affected_1;
	affected[1] = affected_2;
	index[0] = index1;
	index[1] = -1;
	index[2] = -1;
	index[3] = -1;

}

Move::Move(int mtype, RateEnv mrate, Loop *affected_1, int index1, int index2, int index3, int index4) {
	type = mtype;
	rate = mrate;
	affected[0] = affected_1;
	affected[1] = NULL;
	index[0] = index1;
	index[1] = index2;
	index[2] = index3;
	index[3] = index4;

}

Move::Move(int mtype, RateEnv mrate, Loop *affected_1, int *indexarray) {
	type = mtype;
	rate = mrate;
	affected[0] = affected_1;
	affected[1] = NULL;
	for (int loop = 0; loop < 4; loop++)
		index[loop] = indexarray[loop];

}

Move::~Move(void) {
	/* destruction of a move does not imply destruction of the associated loops */
	affected[0] = affected[1] = NULL;
}

double Move::getRate(void) {
	return rate.rate;
}

int Move::getType(void) {
	return type;
}

double Move::getArrType(void) {
	return rate.arrType;
}

Loop *Move::getAffected(int index) {
	assert(index >= 0 && index <= 1);
	return affected[index];
}

Loop *Move::doChoice(void) {
	if (!(type & MOVE_DELETE)) {
		Loop *newLoop = NULL;
		affected[0]->doChoice(this, &newLoop);
		affected[0]->cleanupAdjacent();
		delete affected[0];
		return newLoop;
	}
	if (type & MOVE_DELETE) {
		return Loop::performDeleteMove(this);
	} else
		return NULL;
}

string Move::toString(bool useArr) {

	std::stringstream ss;

	// FD: print even if the rate equals 0

	ss << "type=" << type << " desc:  ";
	ss << utility::moveType(type) << " ";

	if (affected[0] != NULL) {
		ss << affected[0]->toStringShort() << " ";
	}

	if (affected[1] != NULL) {
		ss << affected[1]->toStringShort();
	}
	ss << ", ";
	ss << "(" << index[0] << ", " << index[1] << ", " << index[2] << ", " << index[3] << "),  ";

	ss << rate.toString(useArr) << "\n";

	string output = ss.str();

	return output;

}

///* MoveTree info */
//MoveTree::~MoveTree(void) {
//	/* destruction of a move tree node does imply destruction of all subnodes */
//	if (left != NULL)
//		delete left;
//	if (right != NULL)
//		delete right;
//}
//
///* MoveList */

MoveList::MoveList(int initial_size) {
	totalrate = 0.0;
	moves_size = initial_size;
	if (moves_size >= 1) {
		moves = new Move *[moves_size];

		for (moves_index = 0; moves_index < moves_size; moves_index++)
			moves[moves_index] = NULL;

		moves_index = 0;
	} else {
		moves = NULL;
		moves_size = 0;
		moves_index = 0;
	}

	del_moves = NULL;
	del_moves_size = 0;
	del_moves_index = 0;
	int_index = 0;
}

MoveList::~MoveList(void) {
	int iter = 0;
	// must remove all moves in the movelist.
	while (iter < moves_index) {
		if (moves[iter] != NULL) {
			delete moves[iter];
			moves[iter] = NULL;
		}
		iter++;
	}
	if (moves != NULL)
		delete[] moves;

	iter = 0;
	while (iter < del_moves_index) {
		if (del_moves[iter] != NULL) {
			delete del_moves[iter];
			del_moves[iter] = NULL;
		}
		iter++;
	}
	if (del_moves != NULL)
		delete[] del_moves;
}

void MoveList::resetDeleteMoves(void) {
	int iter = 0;
	while (iter < del_moves_index) {
		if (del_moves[iter] != NULL) {
			delete del_moves[iter];
			del_moves[iter] = NULL;
		}
		iter++;
	}
	if (del_moves != NULL)
		delete[] del_moves;

	if (del_moves == NULL) {
		del_moves = new Move *[2];
		del_moves_size = 2;
	}
	del_moves_index = 0;
}

void MoveList::printAllMoves(bool useArr) {

	for (int i = 0; i < moves_index; i++) {

		cout << "Move" << i << " ";
		cout << moves[i]->toString(useArr);

	}

	for (int i = 0; i < del_moves_index; i++) {

		cout << "Move" << i + moves_index << " ";
		cout << del_moves[i]->toString(useArr);

	}

}

void MoveList::addMove(Move *newmove) {
	int type = newmove->getType();

	if (!(type & MOVE_DELETE)) {

		moves[moves_index] = newmove;
		totalrate += newmove->getRate();
		moves_index++;
		// next add would be an overflow.
		if (moves_index == moves_size) {
			Move **temp = moves;
			moves = new Move *[moves_size * 2];
			for (int loop = 0; loop < moves_size * 2; loop++) {
				if (loop < moves_size) {
					moves[loop] = temp[loop];
					temp[loop] = NULL;
				} else
					moves[loop] = NULL;
			}
			moves_size = moves_size * 2;
			delete[] temp;
		}
	} else // if type is delete move
	{
		totalrate += newmove->getRate();
		if (del_moves == NULL) {
			del_moves = new Move *[2];
			del_moves_size = 2;
			del_moves_index = 0;
		}
		// we do the bounds checking in the opposite order as before, as we want the size for deletion moves to default to 2 and not expand except when we are assured we need more, as opposed to creation moves where we pretty much always will need more. Though perhaps I should change it there also so as to remove an extra operation in the end position cases.
		if (del_moves_index == del_moves_size) {
			Move **temp = del_moves;
			del_moves = new Move *[del_moves_size * 2];
			for (int loop = 0; loop < del_moves_size * 2; loop++) {
				if (loop < del_moves_size) {
					del_moves[loop] = temp[loop];
					temp[loop] = NULL;
				} else
					del_moves[loop] = NULL;
			}
			del_moves_size = del_moves_size * 2;
			delete[] temp;
		}

		del_moves[del_moves_index] = newmove;
		del_moves_index++;
	}
}

Move *MoveList::getMove(Move *iterator) {
	if (iterator == NULL)
		int_index = 0;

	if (int_index == moves_size)
		return NULL;

	return moves[int_index++];
}

Move *MoveList::getChoice(SimTimer& timer) {

	double tmp;

	for (int index = 0; index < moves_index + del_moves_index; index++) {

		if (index < moves_index) {

			tmp = moves[index]->getRate();

		} else {
			tmp = del_moves[index - moves_index]->getRate();

		}
		if (timer.wouldBeHit(tmp) && index < moves_index) {

			return moves[index];

		} else if (timer.wouldBeHit(tmp)) {

			return del_moves[index - moves_index];

		} else {
			timer.checkHit(tmp);
		}
	}
//	assert(*rnd > 0);
	assert(0); // should never call for a move from a container unless it will get one.
	return NULL;
}
/*

 MoveContainer

 */

MoveContainer::MoveContainer(void) {
	totalrate = 0.0;
}

MoveContainer::~MoveContainer(void) {
	totalrate = 0.0;
}

double MoveContainer::getRate(void) {
	return totalrate;
}

int MoveList::getCount(void) {

	return moves_index + del_moves_index;

}

