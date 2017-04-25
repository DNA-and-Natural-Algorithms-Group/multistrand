/*
 Copyright (c) 2007-20016 Caltech. All rights reserved.
 Programming: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 Frits Danennberg (fdann@dna.caltech.edu)
 */

#include <stdio.h>
#include <assert.h>
#include "loop.h"
#include <typeinfo>

#include "utility.h"
#include "moveutil.h"
#include <simoptions.h>
#include <energyoptions.h>

using std::string;

EnergyModel* Loop::energyModel = NULL;

struct RateArr;

inline double Loop::getEnergy(void) {

	if (energyComputed)
		return energy;
	else {
		calculateEnergy();
		energyComputed = true;
		return energy;
	}
}

char Loop::getType(void) {
	return identity;
}

Loop* Loop::getAdjacent(int index) {
	if (adjacentLoops != NULL)
		return adjacentLoops[index];
	return NULL;
}

int Loop::getCurAdjacent(void) {
	return curAdjacent;
}

int Loop::getNumAdjacent(void) {

	return numAdjacent;

}

void Loop::cleanupAdjacent(void) {
	for (int loop = 0; loop < curAdjacent; loop++) {
		adjacentLoops[loop] = NULL;
	}
	curAdjacent = 0;
}

double Loop::returnEnergies(Loop *comefrom) {

	double total = getEnergy();

	for (int loop = 0; loop < curAdjacent; loop++) {

		if (adjacentLoops[loop] != comefrom) {
			total = total + adjacentLoops[loop]->returnEnergies(this);
		}

		assert(adjacentLoops[loop] != NULL);

	}
	return total;
}

double Loop::returnFlux(Loop *comefrom) {

	double total = 0.0;

	if (moves != NULL) {
		total = moves->getRate();
	}

	for (int loop = 0; loop < curAdjacent; loop++) {

		if (adjacentLoops[loop] != comefrom) {
			total = total + adjacentLoops[loop]->returnFlux(this);
		}

		assert(adjacentLoops[loop] != NULL);
	}
	return total;
}

void Loop::firstGen(Loop *comefrom) {

	generateMoves();

	for (int loop = 0; loop < curAdjacent; loop++) {
		if (adjacentLoops[loop] != comefrom && adjacentLoops[loop] != NULL)
			// shouldn't happen, being careful.
			adjacentLoops[loop]->firstGen(this);
		assert(adjacentLoops[loop] != NULL);
	}

}

inline double Loop::getTotalRate(void) {
	return totalRate;
}

Loop::Loop(void) {
	numAdjacent = 0;
	curAdjacent = 0;
	energy = 0.0;
	totalRate = 0.0;
	add_index = 0;

	adjacentLoops = NULL;
	moves = NULL;
	identity = '\0';

}

Loop::~Loop(void) {
	int counter;
	if (adjacentLoops != NULL) {
		for (counter = 0; counter < curAdjacent; counter++) {
			if (adjacentLoops[counter] != NULL) {
				adjacentLoops[counter]->replaceAdjacent(this, NULL);
				adjacentLoops[counter] = NULL; // Hah, take that!
			}
		}
		delete[] adjacentLoops;
		adjacentLoops = NULL;
	}
	if (moves != NULL) {
		delete moves;
		moves = NULL;
	}
}

void Loop::initAdjacency(int index) {
	add_index = index;
}

void Loop::addAdjacent(Loop *loopToAdd) {
	if (curAdjacent >= numAdjacent) // we have already filled up the loops.
									// add an error condition.
			{
		fprintf(stderr, "Adjacent loop attempted to be added when Loop was already full.\n");
		return;
	}
	if (adjacentLoops != NULL) // space has been allocated for the pointers
	{
		adjacentLoops[(add_index + curAdjacent) % numAdjacent] = loopToAdd;
		curAdjacent++;
	}
}

int Loop::replaceAdjacent(Loop *loopToReplace, Loop *loopToReplaceWith) {
	for (int index = 0; index < curAdjacent; index++) {
		if (adjacentLoops[index] == loopToReplace) {
			adjacentLoops[index] = loopToReplaceWith;
			return 1;
		}
	}
	return 0;
}

void Loop::SetEnergyModel(EnergyModel *newEnergyModel) {
	energyModel = newEnergyModel;
}

EnergyModel *Loop::GetEnergyModel(void) {
	return energyModel;
}

void Loop::performComplexSplit(Move *move, Loop **firstOpen, Loop **secondOpen) {
	//  return;

	OpenLoop *tempLoop[2], *newLoop;
	int index[2], sizes[2], loop, flipflop = 0, temp;

	int *pairtypes = NULL;
	int *sidelen = NULL;
	char **seqs = NULL;

	index[1] = -1;

	tempLoop[0] = (OpenLoop*) move->affected[0];
	tempLoop[1] = (OpenLoop*) move->affected[1];

	index[0] = move->index[0];
	for (loop = 0; loop < tempLoop[1]->numAdjacent; loop++)
		if (tempLoop[1]->adjacentLoops[loop] == tempLoop[0])
			index[1] = loop;

	assert(index[1] != -1);

	sizes[0] = index[0] + tempLoop[1]->numAdjacent - index[1] - 1;
	sizes[1] = index[1] + tempLoop[0]->numAdjacent - index[0] - 1;

	for (flipflop = 0; flipflop < 2; flipflop++) {

		pairtypes = new int[sizes[flipflop]];
		sidelen = new int[sizes[flipflop] + 1];
		seqs = new char *[sizes[flipflop] + 1];

		for (loop = 0; loop < sizes[flipflop] + 1; loop++) {
			if (loop < index[flipflop]) {
				if (loop < sizes[flipflop]) // should never be false?
					pairtypes[loop] = tempLoop[flipflop]->pairtype[loop];
				sidelen[loop] = tempLoop[flipflop]->sidelen[loop];
				seqs[loop] = tempLoop[flipflop]->seqs[loop];

			}
			if (loop == index[flipflop]) {
				if (loop < sizes[flipflop])
					pairtypes[loop] = tempLoop[1 - flipflop]->pairtype[loop - index[flipflop] + index[1 - flipflop] + 1];
				sidelen[loop] = tempLoop[flipflop]->sidelen[loop] + tempLoop[1 - flipflop]->sidelen[index[1 - flipflop] + 1] + 1;
				seqs[loop] = tempLoop[flipflop]->seqs[loop];
			}
			if (loop > index[flipflop]) {
				if (loop < sizes[flipflop])
					pairtypes[loop] = tempLoop[1 - flipflop]->pairtype[loop - index[flipflop] + index[1 - flipflop] + 1];
				sidelen[loop] = tempLoop[1 - flipflop]->sidelen[loop - index[flipflop] + index[1 - flipflop] + 1];
				seqs[loop] = tempLoop[1 - flipflop]->seqs[loop - index[flipflop] + index[1 - flipflop] + 1];
			}
		}
		//initialize the new openloops, and connect them correctly, then initialize their moves, etc.
		newLoop = new OpenLoop(sizes[flipflop], pairtypes, sidelen, seqs);

		for (loop = 0; loop < sizes[flipflop]; loop++) {
			if (loop < index[flipflop]) {
				newLoop->addAdjacent(tempLoop[flipflop]->adjacentLoops[loop]);
				temp = tempLoop[flipflop]->adjacentLoops[loop]->replaceAdjacent(tempLoop[flipflop], newLoop);
				assert(temp > 0);
			} else {
				newLoop->addAdjacent(tempLoop[1 - flipflop]->adjacentLoops[loop - index[flipflop] + index[1 - flipflop] + 1]);

				temp = tempLoop[1 - flipflop]->adjacentLoops[loop - index[flipflop] + index[1 - flipflop] + 1]->replaceAdjacent(tempLoop[1 - flipflop],
						newLoop);
				assert(temp > 0);

			}
		}
		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < sizes[flipflop]; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		if (flipflop == 0)
			*firstOpen = newLoop;
		else
			*secondOpen = newLoop;
	}
	//printf("O/O performed, new energies: %lf  %lf\n", (*firstOpen)->getEnergy(), (*secondOpen)->getEnergy() );
	return;
}

string identityToString(char loop) {

	switch (loop) {

	case 'O':
		return "openloop";
	case 'S':
		return "stackloop";
	case 'M':
		return "multiloop";
	case 'B':
		return "bulgeloop";
	case 'I':
		return "interiorloop";
	case 'H':
		return "hairpinloop";

	}

	return "Could not identify loop";
}

//void Loop::setPrimeRates(bool input) {
//
//	energyModel->simOptions->setPrimeRates(input);
//
//}

string Loop::toString() {

	std::stringstream ss;

	ss << "\n** " << identityToString(identity);

	ss << " adjacent ";

	for (int i = 0; i < numAdjacent; i++) {

		ss << identityToString(adjacentLoops[i]->identity) << ",";

	}

	ss << " dG=" << std::setprecision(3) << energy << "\n";

	ss << "** ";

	ss << this->typeInternalsToString();

	return ss.str();

}

string Loop::toStringShort(void) {

	return identityToString(identity);

}

void Loop::printAllMoves(Loop* from) {

	// Doing a short version of the print here
	std::cout << toStringShort();
	std::cout << "\n";

//	moves->printAllMoves(energyModel->simOptions->usePrimeRates);
	moves->printAllMoves(energyModel->useArrhenius());

	for (int i = 0; i < numAdjacent; i++) {

		if (adjacentLoops[i] != from) {

			adjacentLoops[i]->printAllMoves(this);

		}

	}

}

void Loop::generateAndSaveDeleteMove(Loop* input, int position) {

	RateArr tempRate = Loop::generateDeleteMoveRate(this, input);

	double rate = energyModel->applyPrefactors(tempRate.rate, tempRate.left, tempRate.right);

	// if rate is less than zero, nuke the program.
	assert(rate >= 0.0);

	if (rate > 0.0) {

		RateEnv rateEnv = RateEnv(tempRate.rate, energyModel, tempRate.left, tempRate.right);

		moves->addMove(new Move(MOVE_DELETE | MOVE_1, rateEnv, this, input, position));

	}

}

RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end) {

	double tempRate;
	double new_energy, old_energy;
	MoveType left = stackMove;
	MoveType right = stackMove;

	if (start->identity == 'S' && end->identity == 'S') {

		StackLoop *start_ = (StackLoop *) start;
		StackLoop *end_ = (StackLoop *) end;
		int s_index = 0;
		int e_index = 0;

		for (int loop = 0; loop < 2; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				e_index = loop;
			}
		}
		// resulting will be an interior loop with sidelengths 1 and 1.
		// we must get the mismatches correct for this to come out right
		// as they are special cases in the energy model, for 1,1.

		new_energy = energyModel->InteriorEnergy(start_->seqs[s_index], end_->seqs[e_index], 1, 1);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: Two stack loops are created by three basepairs.
		// FD: We query the local context of the middle pair;
		// FD: This is simply two stack environments.

		right = stackMove;
		left = stackMove;

		return RateArr(tempRate / 2.0, left, right);
	}

	if ((start->identity == 'S' && end->identity == 'I') || (start->identity == 'I' && end->identity == 'S')) {

		StackLoop *start_;
		InteriorLoop *end_;

		int s_index = 0, e_index = 0;

		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (InteriorLoop *) end;
		} else {
			start_ = (StackLoop *) end;
			end_ = (InteriorLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				e_index = loop;
			}
		}
		// resulting will be an interior loop side lengths equal to the length
		// of the 'input' bulge loop, plus one on each side (ie, one side will be B+1, the other 1.
		//      if(s_index == e_index )
		//	fprintf(stderr, "ERROR: Misaligned loops in loop.cc generateDeleteMoveRate S/B, values: %d %d\n",s_index,e_index);

		new_energy = energyModel->InteriorEnergy(start_->seqs[s_index], end_->int_seq[e_index], end_->sizes[1 - e_index] + 1, end_->sizes[e_index] + 1);
		old_energy = start->getEnergy() + end->getEnergy();

		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: A base-pair is present between a stacking loop and an interior loop.
		// FD: We query the local context of the middle pair;
		// FD: This is a stack and loop environment.

		right = stackMove;
		left = loopMove;

		return RateArr(tempRate / 2.0, left, right);
	}

	if ((start->identity == 'S' && end->identity == 'B') || (start->identity == 'B' && end->identity == 'S')) {
		StackLoop *start_;
		BulgeLoop *end_;
		int s_index = 0, e_index = 0;

		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (BulgeLoop *) end;
		} else {
			start_ = (StackLoop *) end;
			end_ = (BulgeLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {

				s_index = loop;

			}
			if (end_->adjacentLoops[loop] != start_) {

				e_index = loop;

			}
		}
		// resulting will be an interior loop side lengths equal to the length
		// of the 'input' bulge loop, plus one on each side (ie, one side will be B+1, the other 1.
		//      if(s_index == e_index )
		//	fprintf(stderr, "ERROR: Misaligned loops in loop.cc generateDeleteMoveRate S/B, values: %d %d\n",s_index,e_index);

		new_energy = energyModel->InteriorEnergy(start_->seqs[s_index], end_->bulge_seq[e_index], end_->bulgesize[1 - e_index] + 1,
				end_->bulgesize[e_index] + 1);
		old_energy = start->getEnergy() + end->getEnergy();

		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: A base pair is present between a stacking loop and a bulge loop.
		// FD: We query the local context of the middle pair;

		right = stackMove;
		left = stackLoopMove;

		return RateArr(tempRate / 2.0, left, right);
	}

	if ((start->identity == 'S' && end->identity == 'H') || (start->identity == 'H' && end->identity == 'S')) {
		StackLoop *start_;
		HairpinLoop *end_;
		int s_index = 0;

		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (HairpinLoop *) end;
		} else {
			start_ = (StackLoop *) end;
			end_ = (HairpinLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
		}
		// end is the hairpin, which has no extra adjacencies.

		// resulting will be a hairpin loop equal to the previous plus an extra base on each side.

		new_energy = energyModel->HairpinEnergy(start_->seqs[s_index], end_->hairpinsize + 2);
		old_energy = start->getEnergy() + end->getEnergy();

		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: A base pair is present between a stacking loop and a hairpin loop.
		// FD: We query the local context of the middle pair;

		right = stackMove;
		left = loopMove;

		return RateArr(tempRate / 2.0, left, right);
	}

	if ((start->identity == 'S' && end->identity == 'M') || (start->identity == 'M' && end->identity == 'S')) {
		StackLoop *start_;
		MultiLoop *end_;
		int s_index = 0, e_index = 0;

		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (StackLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}

		// note e_index has different meaning now for multiloops.
		// FD: e_index is the index of the attached loop for multiloop end_
		// FD: s_index is the index of the attached loop for stackloop start_

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent];
		char **seqs = new char *[end_->numAdjacent];

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				pairtypes[loop] = end_->pairtype[loop];
				if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0)) {
					sidelens[loop] = end_->sidelen[loop];
				} else {
					sidelens[loop] = end_->sidelen[loop] + 1;
				}
				seqs[loop] = end_->seqs[loop];
			} else {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + 1;
				seqs[loop] = start_->seqs[s_index];
			}
		}

		new_energy = energyModel->MultiloopEnergy(end_->numAdjacent, sidelens, seqs);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		left = stackMove;
		right = energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);

		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;

		return RateArr(tempRate / 2.0, left, right);
	}

	if ((start->identity == 'S' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'S')) {

		StackLoop *start_ = (StackLoop *) end;
		OpenLoop *end_ = (OpenLoop *) start;

		int s_index = 0, e_index = 0;

		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (OpenLoop *) end;
		}

		for (int loop = 0; (loop < end_->numAdjacent) || (loop < 2); loop++) {

			if (loop <= 1) {
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			}
			if (loop < end_->numAdjacent && (end_->adjacentLoops != NULL) && end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for openloops.

//		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent + 1];
		char **seqs = new char *[end_->numAdjacent + 1];

		for (int loop = 0; loop < end_->numAdjacent + 1; loop++) {
			if (loop == e_index) {
//				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + 1;
				seqs[loop] = end_->seqs[loop];
			} else if (loop == e_index + 1) {
//				if (loop < end_->numAdjacent)
//					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop] + 1;
				seqs[loop] = start_->seqs[s_index];
			} else {
//				if (loop < end_->numAdjacent)
//					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop];
				seqs[loop] = end_->seqs[loop];
			}
		}

		new_energy = energyModel->OpenloopEnergy(end_->numAdjacent, sidelens, seqs);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: A stack and an open loop. e_index is the location of the stack.

		left = energyModel->prefactorOpen(e_index, (end_->numAdjacent + 1), end_->sidelen);
		right = stackMove;

//		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;

		return RateArr(tempRate / 2.0, left, right);

	}

	if (start->identity == 'I' && end->identity == 'I') {
		InteriorLoop *end_ = (InteriorLoop *) end, *start_ = (InteriorLoop *) start;
		Loop *start_extra, *end_extra;
		int s_index = 0, e_index = 0;

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				start_extra = start_->adjacentLoops[loop];
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				end_extra = end_->adjacentLoops[loop];
				e_index = loop;
			}
		}

		new_energy = energyModel->InteriorEnergy(start_->int_seq[s_index], end_->int_seq[e_index], end_->sizes[1 - e_index] + start_->sizes[s_index] + 1,
				end_->sizes[e_index] + start_->sizes[1 - s_index] + 1);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: interior loop and interior loop, this has to be loopMove and loopMove;

		left = loopMove;
		right = loopMove;

		return RateArr(tempRate / 2.0, left, right);

	}

	if ((start->identity == 'I' && end->identity == 'B') || (start->identity == 'B' && end->identity == 'I')) {
		InteriorLoop *start_;
		BulgeLoop *end_;
		Loop *start_extra, *end_extra;
		int s_index = 0, e_index = 0;

		if (start->identity == 'S') {
			start_ = (InteriorLoop *) start;
			end_ = (BulgeLoop *) end;
		} else {
			start_ = (InteriorLoop *) end;
			end_ = (BulgeLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				start_extra = start_->adjacentLoops[loop];
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				end_extra = end_->adjacentLoops[loop];
				e_index = loop;
			}
		}
		// resulting will be an interior loop side lengths equal to the length
		// of the 'input' bulge loop, plus one on each side (ie, one side will be B+1, the other 1.
		//      if(s_index == e_index )
		//	fprintf(stderr, "ERROR: Misaligned loops in loop.cc generateDeleteMoveRate I/B, values: %d %d\n",s_index,e_index);

		new_energy = energyModel->InteriorEnergy(start_->int_seq[s_index], end_->bulge_seq[e_index], end_->bulgesize[1 - e_index] + 1 + start_->sizes[s_index],
				end_->bulgesize[e_index] + 1 + start_->sizes[1 - s_index]);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: interior loop and bulge loop, this has to be loopMove and stackLoopMove;

		left = loopMove;
		right = stackLoopMove;

		return RateArr(tempRate / 2.0, left, right);

	}

	if ((start->identity == 'I' && end->identity == 'H') || (start->identity == 'H' && end->identity == 'I')) {
		InteriorLoop *start_;
		HairpinLoop *end_;
		int s_index = 0;

		if (start->identity == 'I') {
			start_ = (InteriorLoop *) start;
			end_ = (HairpinLoop *) end;
		} else {
			start_ = (InteriorLoop *) end;
			end_ = (HairpinLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
		}
		// end is the hairpin, which has no extra adjacencies.

		// resulting will be a hairpin loop with previous size, plus interior loop's sizes (both) plus 2 (for the pairing that's now unpaired)

		new_energy = energyModel->HairpinEnergy(start_->int_seq[s_index], end_->hairpinsize + 2 + start_->sizes[0] + start_->sizes[1]);
		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: interior loop and hairpin loop, this has to be loopMove and loopMove;

		left = loopMove;
		right = loopMove;

		return RateArr(tempRate / 2.0, left, right);

	}

	if ((start->identity == 'I' && end->identity == 'M') || (start->identity == 'M' && end->identity == 'I')) {
		InteriorLoop *start_;
		MultiLoop *end_;
		int s_index = 0, e_index = 0;

		if (start->identity == 'I') {
			start_ = (InteriorLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (InteriorLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent];
		char **seqs = new char *[end_->numAdjacent];

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				pairtypes[loop] = end_->pairtype[loop];
				if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
					sidelens[loop] = end_->sidelen[loop];
				else
					sidelens[loop] = end_->sidelen[loop] + start_->sizes[1 - s_index] + 1;
				seqs[loop] = end_->seqs[loop];
			} else {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + start_->sizes[s_index] + 1;
				seqs[loop] = start_->int_seq[s_index];
			}
		}

		new_energy = energyModel->MultiloopEnergy(end_->numAdjacent, sidelens, seqs);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: interior loop and multi loop, this has to be loopMove and something else;

		left = energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);
		right = loopMove;

		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;

		return RateArr(tempRate / 2.0, left, right);

	}

	if ((start->identity == 'I' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'I')) {
		InteriorLoop *start_;
		OpenLoop *end_;
		int s_index = 0, e_index = 0;

		if (start->identity == 'I') {
			start_ = (InteriorLoop *) start;
			end_ = (OpenLoop *) end;
		} else {
			start_ = (InteriorLoop *) end;
			end_ = (OpenLoop *) start;
		}

		for (int loop = 0; (loop < end_->numAdjacent) || (loop < 2); loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (loop < end_->numAdjacent)
				if (end_->adjacentLoops[loop] == start_) {
					e_index = loop;
				}
		}
		// note e_index has different meaning now for openloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent + 1];
		char **seqs = new char *[end_->numAdjacent + 1];

		for (int loop = 0; loop < end_->numAdjacent + 1; loop++) {
			if (loop == e_index) {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + 1 + start_->sizes[1 - s_index];
				seqs[loop] = end_->seqs[loop];
			} else if (loop == e_index + 1) {
				if (loop < end_->numAdjacent)
					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop] + 1 + start_->sizes[s_index];
				seqs[loop] = start_->int_seq[s_index];
			} else {
				if (loop < end_->numAdjacent)
					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop];
				seqs[loop] = end_->seqs[loop];
			}
		}

		new_energy = energyModel->OpenloopEnergy(end_->numAdjacent, sidelens, seqs);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: interior loop and open loop, this has to be loopMove and something else;
		if (energyModel->useArrhenius()) {

			left = energyModel->prefactorOpen(e_index, end_->numAdjacent + 1, end_->sidelen);
			right = loopMove;

		}

		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;

		return RateArr(tempRate / 2.0, left, right); //		return tempRate / 2.0;
	}

// start bulge

	if (start->identity == 'B' && end->identity == 'B') {
		BulgeLoop *end_ = (BulgeLoop *) end, *start_ = (BulgeLoop *) start;
		Loop *start_extra, *end_extra;
		int s_index = 0, e_index = 0;

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				start_extra = start_->adjacentLoops[loop];
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				end_extra = end_->adjacentLoops[loop];
				e_index = loop;
			}
		}

		new_energy = energyModel->InteriorEnergy(start_->bulge_seq[s_index], end_->bulge_seq[e_index],
				end_->bulgesize[1 - e_index] + start_->bulgesize[s_index] + 1, end_->bulgesize[e_index] + start_->bulgesize[1 - s_index] + 1);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: bulge loop and bulge loop, this has to be loopMove and loopMove ;

		left = stackLoopMove;
		right = stackLoopMove;

		return RateArr(tempRate / 2.0, left, right);

	}

	if ((start->identity == 'B' && end->identity == 'H') || (start->identity == 'H' && end->identity == 'B')) {
		BulgeLoop *start_;
		HairpinLoop *end_;
		int s_index = 0;

		if (start->identity == 'B') {
			start_ = (BulgeLoop *) start;
			end_ = (HairpinLoop *) end;
		} else {
			start_ = (BulgeLoop *) end;
			end_ = (HairpinLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
		}
		// end is the hairpin, which has no extra adjacencies.

		// resulting will be a hairpin loop with previous size, plus interior loop's sizes (both) plus 2 (for the pairing that's now unpaired)

		new_energy = energyModel->HairpinEnergy(start_->bulge_seq[s_index], end_->hairpinsize + 2 + start_->bulgesize[0] + start_->bulgesize[1]);
		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: bulge loop and hairpin loop, this has to be stackLoopMove and loopMove ;

		left = stackLoopMove;
		right = loopMove;

		return RateArr(tempRate / 2.0, left, right);

	}

	if ((start->identity == 'B' && end->identity == 'M') || (start->identity == 'M' && end->identity == 'B')) {
		BulgeLoop *start_;
		MultiLoop *end_;
		int s_index = 0, e_index = 0;

		if (start->identity == 'B') {
			start_ = (BulgeLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (BulgeLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent];
		char **seqs = new char *[end_->numAdjacent];

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				pairtypes[loop] = end_->pairtype[loop];
				if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
					sidelens[loop] = end_->sidelen[loop];
				else
					sidelens[loop] = end_->sidelen[loop] + start_->bulgesize[1 - s_index] + 1;
				seqs[loop] = end_->seqs[loop];
			} else {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + start_->bulgesize[s_index] + 1;
				seqs[loop] = start_->bulge_seq[s_index];
			}
		}

		new_energy = energyModel->MultiloopEnergy(end_->numAdjacent, sidelens, seqs);
		old_energy = start->getEnergy() + end->getEnergy();

		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: bulge loop and multi loop, this has to be stackLoopMove and something else;
		if (energyModel->useArrhenius()) {

			left = energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);
			right = stackLoopMove;

		}

		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;
		return RateArr(tempRate / 2.0, left, right);

	}

	if ((start->identity == 'B' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'B')) {
		BulgeLoop *start_;
		OpenLoop *end_;
		int s_index = 0, e_index = 0;

		if (start->identity == 'B') {
			start_ = (BulgeLoop *) start;
			end_ = (OpenLoop *) end;
		} else {
			start_ = (BulgeLoop *) end;
			end_ = (OpenLoop *) start;
		}

		for (int loop = 0; (loop < end_->numAdjacent) || (loop < 2); loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (loop < end_->numAdjacent)
				if (end_->adjacentLoops[loop] == start_) {
					e_index = loop;
				}
		}
		// note e_index has different meaning now for openloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent + 1];
		char **seqs = new char *[end_->numAdjacent + 1];

		for (int loop = 0; loop < end_->numAdjacent + 1; loop++) {
			if (loop == e_index) {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + 1 + start_->bulgesize[1 - s_index];
				seqs[loop] = end_->seqs[loop];
			} else if (loop == e_index + 1) {
				if (loop < end_->numAdjacent)
					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop] + 1 + start_->bulgesize[s_index];
				seqs[loop] = start_->bulge_seq[s_index];
			} else {
				if (loop < end_->numAdjacent)
					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop];
				seqs[loop] = end_->seqs[loop];
			}
		}

		new_energy = energyModel->OpenloopEnergy(end_->numAdjacent, sidelens, seqs);
		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// FD: bulge loop and open loop, this has to be stackLoopMove and something else;
		if (energyModel->useArrhenius()) {

			left = energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);
			right = stackLoopMove;

		}

		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;

		return RateArr(tempRate / 2.0, left, right);

	}

// end bulge

// start hairpin

	if (start->identity == 'H' && end->identity == 'H') {
		fprintf(stderr, "Hairpin/Hairpin deletion move encountered - not currently supported.\n");
		assert(0);
		return RateArr(-1.0, left, right);
	}

	// hairpin,
	if ((start->identity == 'H' && end->identity == 'M') || (start->identity == 'M' && end->identity == 'H')) {
		HairpinLoop *start_;
		MultiLoop *end_;
		int e_index = 0;

		if (start->identity == 'H') {
			start_ = (HairpinLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (HairpinLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for multiloops.

		// Several cases we now need to check for:
		// 1. Multiloop had 3 adjacent stacks, and the resulting figure with two adjacent stacks is now a interior loop.
		// 2. Multiloop had 3 adjacent stacks, and the resulting structure with two adjacent stacks is now a bulge loop
		// 3. Multiloop had >3 adjacent stacks, and the resulting structure is a multiloop with one less adjacent stack.

		if (end_->numAdjacent == 3) {
			int sizes[2];
			int positions[2] = { (e_index + 1) % 3, (e_index + 2) % 3 };
			if (e_index == 0)
				sizes[0] = end_->sidelen[0] + end_->sidelen[2] + start_->hairpinsize + 2;
			else
				sizes[0] = end_->sidelen[e_index] + end_->sidelen[e_index - 1] + start_->hairpinsize + 2;

			sizes[1] = end_->sidelen[(e_index + 1) % 3];

			if (sizes[1] > 0) // interior loop case
					{

				new_energy = energyModel->InteriorEnergy(end_->seqs[positions[1]], end_->seqs[positions[0]], sizes[0], sizes[1]);
				old_energy = start->getEnergy() + end->getEnergy();
				tempRate = energyModel->returnRate(old_energy, new_energy, 0);

				// interior loop, so could be loopMove plus loopMove, or plus stackloopMove
				if (energyModel->useArrhenius()) {

					left = energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);
					right = loopMove;

				}

				return RateArr(tempRate / 2.0, left, right);

			} else  // bulge loop case
			{
				// positions 0 and 1 are the two remaining stacks in the bulge loop now, 0 being the one directly after the removed stack, and 1 being the one before. Thus, i is the 0'th base at position 1, j is the sizes[1]+1 base at position 0, p is the sizes[0]+1 base at position 1, and q is the 0th base at position 0.
				new_energy = energyModel->BulgeEnergy(end_->seqs[positions[1]][0], end_->seqs[positions[0]][sizes[1] + 1],
						end_->seqs[positions[1]][sizes[0] + 1], end_->seqs[positions[0]][0], sizes[0]);
				old_energy = start->getEnergy() + end->getEnergy();
				tempRate = energyModel->returnRate(old_energy, new_energy, 0);

				// bulge loop, so loopMove plus stackStackMove

				left = loopMove;
				right = stackStackMove;

				return RateArr(tempRate / 2.0, left, right);

			}
		}

		else if (end_->numAdjacent > 3)  // multiloop case
				{
			int *pairtypes = new int[end_->numAdjacent - 1];
			int *sidelens = new int[end_->numAdjacent - 1];
			char **seqs = new char *[end_->numAdjacent - 1];

			for (int loop = 0; loop < end_->numAdjacent; loop++) {
				if (loop != e_index) {
					if (loop < e_index) {
						pairtypes[loop] = end_->pairtype[loop];
						if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
							sidelens[loop] = end_->sidelen[loop];
						else
							sidelens[loop] = end_->sidelen[loop] + end_->sidelen[e_index] + 2 + start_->hairpinsize;
						seqs[loop] = end_->seqs[loop];
					}
					if (loop > e_index) {
						pairtypes[loop - 1] = end_->pairtype[loop];
						if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
							sidelens[loop - 1] = end_->sidelen[loop];
						else
							sidelens[loop - 1] = end_->sidelen[loop] + end_->sidelen[e_index] + 2 + start_->hairpinsize;
						seqs[loop - 1] = end_->seqs[loop];

					}
				}
			}

			new_energy = energyModel->MultiloopEnergy(end_->numAdjacent - 1, sidelens, seqs);

			old_energy = start->getEnergy() + end->getEnergy();
			tempRate = energyModel->returnRate(old_energy, new_energy, 0);

			// multi loop, so loopMove plus something else
			if (energyModel->useArrhenius()) {

				left = loopMove; //energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);
				right = stackStackMove;

			}

			delete[] pairtypes;
			delete[] sidelens;
			delete[] seqs;

			return RateArr(tempRate / 2.0, left, right);

		}

	}

	if ((start->identity == 'H' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'H')) {
		HairpinLoop *start_ = (HairpinLoop *) end;
		OpenLoop *end_ = (OpenLoop *) start;
		int e_index = 0;

		if (start->identity == 'H') {
			start_ = (HairpinLoop *) start;
			end_ = (OpenLoop *) end;
		}

		for (int loop = 0; (loop < end_->numAdjacent); loop++) {
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for openloops.

		int *pairtypes = new int[end_->numAdjacent - 1];
		int *sidelens = new int[end_->numAdjacent];
		char **seqs = new char *[end_->numAdjacent];

		for (int loop = 0; loop < end_->numAdjacent + 1; loop++) {
			if (loop < e_index) {
				pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop];
				seqs[loop] = end_->seqs[loop];
			} else if (loop == e_index) {
				if (loop < end_->numAdjacent - 1)
					pairtypes[loop] = end_->pairtype[loop + 1];
				sidelens[loop] = end_->sidelen[loop] + 2 + end_->sidelen[loop + 1] + start_->hairpinsize;
				seqs[loop] = end_->seqs[loop];
			} else if (loop > e_index + 1) {
				if (loop < end_->numAdjacent)
					pairtypes[loop - 1] = end_->pairtype[loop];
				sidelens[loop - 1] = end_->sidelen[loop];
				seqs[loop - 1] = end_->seqs[loop];
			}
		}

		new_energy = energyModel->OpenloopEnergy(end_->numAdjacent - 1, sidelens, seqs);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// hairpin loop and openLoop, so loopMove plus something else
		if (energyModel->useArrhenius()) {

			left = loopMove; //energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);
			right = energyModel->prefactorOpen(e_index, end_->numAdjacent + 1, end_->sidelen);
			;

		}

		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;
		return RateArr(tempRate / 2.0, left, right);

	}

// end hairpin

// start multiloop

	if (start->identity == 'M' && end->identity == 'M') {
		MultiLoop *start_;
		MultiLoop *end_;
		int s_index = 0, e_index = 0, index = 0;

		start_ = (MultiLoop *) start;
		end_ = (MultiLoop *) end;

		for (int loop = 0; loop < end_->numAdjacent || loop < start_->numAdjacent; loop++) {
			if (loop < start_->numAdjacent)
				if (start_->adjacentLoops[loop] == end_) {
					s_index = loop;
				}
			if (loop < end_->numAdjacent)
				if (end_->adjacentLoops[loop] == start_) {
					e_index = loop;
				}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent + start_->numAdjacent - 2];
		int *sidelens = new int[end_->numAdjacent + start_->numAdjacent - 2];
		char **seqs = new char *[end_->numAdjacent + start_->numAdjacent - 2];

		index = 0;
		for (int loop = 0; loop < start_->numAdjacent; loop++) {
			if (loop != s_index) {
				pairtypes[index] = start_->pairtype[loop];
				if ((loop != s_index - 1 && s_index != 0) || (loop != start_->numAdjacent - 1 && s_index == 0))
					sidelens[index] = start_->sidelen[loop];
				else
					sidelens[index] = start_->sidelen[loop] + end_->sidelen[e_index] + 1;
				seqs[index] = start_->seqs[loop];
				index++;
			} else {
				for (int loop2 = 1; loop2 < end_->numAdjacent; loop2++) {
					int temp = (e_index + loop2) % end_->numAdjacent;

					pairtypes[index] = end_->pairtype[temp];
					if ((temp != e_index - 1 && e_index != 0) || (temp != end_->numAdjacent - 1 && e_index == 0))
						sidelens[index] = end_->sidelen[temp];
					else
						sidelens[index] = end_->sidelen[temp] + start_->sidelen[s_index] + 1;
					seqs[index] = end_->seqs[temp];
					index++;
				}
			}
		}
		assert(index == end_->numAdjacent + start_->numAdjacent - 2);

		new_energy = energyModel->MultiloopEnergy(start_->numAdjacent + end_->numAdjacent - 2, sidelens, seqs);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// multiloop and multiLoop, so something plus something else
		if (energyModel->useArrhenius()) {

			left = energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);
			right = energyModel->getPrefactorsMulti(s_index, start_->numAdjacent, start_->sidelen);

		}

		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;

		return RateArr(tempRate / 2.0, left, right);

	}

	if ((start->identity == 'M' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'M')) {
		OpenLoop *start_;
		MultiLoop *end_;
		int s_index = 0, e_index = 0, index = 0;

		start_ = (OpenLoop *) start;
		end_ = (MultiLoop *) end;

		if (start->identity == 'O') {
			start_ = (OpenLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (OpenLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent || loop < start_->numAdjacent; loop++) {
			if (loop < start_->numAdjacent)
				if (start_->adjacentLoops[loop] == end_) {
					s_index = loop;
				}
			if (loop < end_->numAdjacent)
				if (end_->adjacentLoops[loop] == start_) {
					e_index = loop;
				}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent + start_->numAdjacent - 2];
		int *sidelens = new int[end_->numAdjacent + start_->numAdjacent - 1];
		char **seqs = new char *[end_->numAdjacent + start_->numAdjacent - 1];

		index = 0;
		for (int loop = 0; loop <= start_->numAdjacent; loop++) {
			if (loop == s_index + 1) {
				int temp = (e_index + end_->numAdjacent - 1) % end_->numAdjacent;

				if (loop < start_->numAdjacent)
					pairtypes[index] = start_->pairtype[loop];
				sidelens[index] = start_->sidelen[loop] + end_->sidelen[temp] + 1;
				seqs[index] = end_->seqs[temp];
				if (loop < start_->numAdjacent)
					index++;
			} else if (loop != s_index) {
				if (loop < start_->numAdjacent)
					pairtypes[index] = start_->pairtype[loop];
				sidelens[index] = start_->sidelen[loop];
				seqs[index] = start_->seqs[loop];
				if (loop < start_->numAdjacent)
					index++;
			} else {
				for (int loop2 = 1; loop2 < end_->numAdjacent; loop2++) {
					int temp = (e_index + loop2) % end_->numAdjacent;
					int temp2 = (e_index + loop2 - 1) % end_->numAdjacent;

					pairtypes[index] = end_->pairtype[temp];
					if (loop2 == 1) {
						sidelens[index] = end_->sidelen[e_index] + start_->sidelen[s_index] + 1;
						seqs[index] = start_->seqs[s_index];
					} else {
						sidelens[index] = end_->sidelen[temp2];
						seqs[index] = end_->seqs[temp2];
					}
					index++;
				}
			}
		}
		assert(index == end_->numAdjacent + start_->numAdjacent - 2);

		new_energy = energyModel->OpenloopEnergy(start_->numAdjacent + end_->numAdjacent - 2, sidelens, seqs);

		old_energy = start->getEnergy() + end->getEnergy();
		tempRate = energyModel->returnRate(old_energy, new_energy, 0);

		// multiloop and openLoop, so something plus something else
		if (energyModel->useArrhenius()) {

			left = energyModel->getPrefactorsMulti(e_index, end_->numAdjacent, end_->sidelen);
			right = energyModel->getPrefactorsMulti(s_index, start_->numAdjacent, start_->sidelen);

		}

		delete[] pairtypes;
		delete[] sidelens;
		delete[] seqs;

		return RateArr(tempRate / 2.0, left, right);

	}

// end multiloop
// this is the DELETE MOVE function.
// start openloop

	if (start->identity == 'O' && end->identity == 'O') {

		OpenLoop *tempLoop[2];
		double new_energies[2] = { 0.0, 0.0 };
		int index[2], sizes[2], loop, flipflop = 0, temp;
		int *pairtypes = NULL;
		int *sidelen = NULL;
		char **seqs = NULL;

		tempLoop[0] = (OpenLoop *) start;
		tempLoop[1] = (OpenLoop *) end;

		for (loop = 0; loop < tempLoop[1]->numAdjacent || loop < tempLoop[0]->numAdjacent; loop++) {
			if (loop < tempLoop[0]->numAdjacent)
				if (tempLoop[0]->adjacentLoops[loop] == tempLoop[1]) {
					index[0] = loop;
				}
			if (loop < tempLoop[1]->numAdjacent)
				if (tempLoop[1]->adjacentLoops[loop] == tempLoop[0]) {
					index[1] = loop;
				}
		}

		sizes[0] = index[0] + tempLoop[1]->numAdjacent - index[1] - 1;
		sizes[1] = index[1] + tempLoop[0]->numAdjacent - index[0] - 1;

		for (flipflop = 0; flipflop < 2; flipflop++) {

			// pairtypes = new int[sizes[flipflop]];
			sidelen = new int[sizes[flipflop] + 1];
			seqs = new char *[sizes[flipflop] + 1];

			for (loop = 0; loop < sizes[flipflop] + 1; loop++) {
				if (loop < index[flipflop]) {
					sidelen[loop] = tempLoop[flipflop]->sidelen[loop];
					seqs[loop] = tempLoop[flipflop]->seqs[loop];

				}
				if (loop == index[flipflop]) {

					sidelen[loop] = tempLoop[flipflop]->sidelen[loop] + tempLoop[1 - flipflop]->sidelen[index[1 - flipflop] + 1] + 1;
					seqs[loop] = tempLoop[flipflop]->seqs[loop];
				}
				if (loop > index[flipflop]) {
					sidelen[loop] = tempLoop[1 - flipflop]->sidelen[loop - index[flipflop] + index[1 - flipflop] + 1];
					seqs[loop] = tempLoop[1 - flipflop]->seqs[loop - index[flipflop] + index[1 - flipflop] + 1];
				}
			}
			//initialize the new openloops, and connect them correctly, then initialize their moves, etc.
			new_energies[flipflop] = energyModel->OpenloopEnergy(sizes[flipflop], sidelen, seqs);

			delete[] sidelen;
			delete[] seqs;
		}

		old_energy = start->getEnergy() + end->getEnergy();
		new_energy = new_energies[0] + new_energies[1];
		//printf("O/O delete rate calculate: old: %lf %lf, new: %lf %lf\n", start->getEnergy(), end->getEnergy(), new_energies[0], new_energies[1] );

		double tempRate = energyModel->returnRate(old_energy, new_energy, 3);

		// openLoop and openLoop, so something plus something else
		if (energyModel->useArrhenius()) {

			left = energyModel->prefactorOpen(index[0], tempLoop[0]->numAdjacent + 1, tempLoop[0]->sidelen);
			right = energyModel->prefactorOpen(index[1], tempLoop[1]->numAdjacent + 1, tempLoop[1]->sidelen);

		}

		return RateArr(tempRate / 2.0, left, right); //		return tempRate;
	}

	return RateArr(-1.0, left, right);
}

Loop *Loop::performDeleteMove(Move *move) {


	Loop* start = move->affected[0];
	Loop* end = move->affected[1];


	if (start->identity == 'S' && end->identity == 'S') {
		StackLoop *start_, *end_;

		Loop *start_extra, *end_extra, *newLoop;
		int s_index = 0, e_index = 0;
		if (((StackLoop *) start)->seqs[0] > ((StackLoop *) end)->seqs[0]) {
			start_ = (StackLoop *) end;
			end_ = (StackLoop *) start;
		} else {
			start_ = (StackLoop *) start;
			end_ = (StackLoop *) end;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				start_extra = start_->adjacentLoops[loop];
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				end_extra = end_->adjacentLoops[loop];
				e_index = loop;
			}
		}

		// resulting will be an interior loop with sidelengths 1 and 1.
		newLoop = new InteriorLoop(start_->pairtype[s_index], end_->pairtype[e_index], 1, 1, start_->seqs[s_index], end_->seqs[e_index]);

		newLoop->addAdjacent(start_->adjacentLoops[s_index]);
		// TODO: fix this! asserts generate no code when NDEBUG is set!
		assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);
		newLoop->addAdjacent(end_->adjacentLoops[e_index]);
		// TODO: fix this too! see above comment.
		assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);

		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();
		end_->adjacentLoops[e_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'S' && end->identity == 'I') || (start->identity == 'I' && end->identity == 'S')) {

		StackLoop *start_;
		InteriorLoop *end_;
		Loop* newLoop;
		int s_index = 0, e_index = 0;
		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (InteriorLoop *) end;
		} else {
			start_ = (StackLoop *) end;
			end_ = (InteriorLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				e_index = loop;
			}
		}

		// resulting will be an interior loop side lengths equal to the length
		// of the 'input' interior loop+1.
		assert(!(start_->seqs[s_index] + 1 == end_->int_seq[e_index]));
		assert(!(start_->seqs[s_index] == end_->int_seq[e_index] + 1));

		if (start_->seqs[s_index] < end_->int_seq[e_index]) {
			newLoop = new InteriorLoop(start_->pairtype[s_index], end_->pairtype[e_index], end_->sizes[1 - e_index] + 1, end_->sizes[e_index] + 1,
					start_->seqs[s_index], end_->int_seq[e_index]);

			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);
			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
		} else {
			newLoop = new InteriorLoop(end_->pairtype[e_index], start_->pairtype[s_index], end_->sizes[e_index] + 1, end_->sizes[1 - e_index] + 1,
					end_->int_seq[e_index], start_->seqs[s_index]);

			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);

		}
		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();
		end_->adjacentLoops[e_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'S' && end->identity == 'B') || (start->identity == 'B' && end->identity == 'S')) {

		StackLoop *start_;
		BulgeLoop *end_;
		Loop* newLoop;
		int s_index = 0, e_index = 0;
		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (BulgeLoop *) end;
		} else {
			start_ = (StackLoop *) end;
			end_ = (BulgeLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				e_index = loop;
			}
		}

		// resulting will be an interior loop side lengths equal to the length
		// of the 'input' interior loop+1.
		if (start_->seqs[s_index] < end_->bulge_seq[e_index]) {
			newLoop = new InteriorLoop(start_->pairtype[s_index], end_->pairtype[e_index], end_->bulgesize[1 - e_index] + 1, end_->bulgesize[e_index] + 1,
					start_->seqs[s_index], end_->bulge_seq[e_index]);

			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);
			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
		} else {
			newLoop = new InteriorLoop(end_->pairtype[e_index], start_->pairtype[s_index], end_->bulgesize[e_index] + 1, end_->bulgesize[1 - e_index] + 1,
					end_->bulge_seq[e_index], start_->seqs[s_index]);

			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);

		}
		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();
		end_->adjacentLoops[e_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'S' && end->identity == 'H') || (start->identity == 'H' && end->identity == 'S')) {

		StackLoop *start_;
		HairpinLoop *end_;
		Loop* newLoop;

		int s_index = 0, e_index = 0;
		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (HairpinLoop *) end;
		} else {
			start_ = (StackLoop *) end;
			end_ = (HairpinLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
		}
		// end has only a single adjacent.

		// resulting will be an hairpin loop, size+2.

		newLoop = new HairpinLoop(start_->pairtype[s_index], end_->hairpinsize + 2, start_->seqs[s_index]);

		//      printf("index: %d size: %d seqs: %s\n",s_index, end_->hairpinsize+2, start_->seqs[s_index]);
		newLoop->addAdjacent(start_->adjacentLoops[s_index]);
		// TODO: fix this! asserts generate no code when NDEBUG is set!
		assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);

		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'S' && end->identity == 'M') || (start->identity == 'M' && end->identity == 'S')) {

		StackLoop *start_;
		MultiLoop *end_;
		Loop *newLoop;
		int s_index = 0, e_index = 0, temp = 0;

		if (start->identity == 'S') {
			start_ = (StackLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (StackLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent];
		char **seqs = new char *[end_->numAdjacent];

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				pairtypes[loop] = end_->pairtype[loop];
				if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
					sidelens[loop] = end_->sidelen[loop];
				else
					sidelens[loop] = end_->sidelen[loop] + 1;
				seqs[loop] = end_->seqs[loop];
			} else {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + 1;
				seqs[loop] = start_->seqs[s_index];
			}
		}

		//      new_energy = energyModel_Primary->MultiloopEnergy( end_->numAdjacent, pairtypes, sidelens, seqs );

		//      old_energy = start->getEnergy() + end->getEnergy();
		//      temprate = energyModel_Primary->returnRate( old_energy, new_energy, 0);
		//      delete[] pairtypes;
		//      delete[] sidelens;
		//      delete[] seqs;

		// resulting will be an multiloop, same# of adjacent helices, two sides longer by one base, and one pairtype possibly changed.
		newLoop = new MultiLoop(end_->numAdjacent, pairtypes, sidelens, seqs);

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				newLoop->addAdjacent(end_->adjacentLoops[loop]);
				temp = end_->adjacentLoops[loop]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);
			} else {
				newLoop->addAdjacent(start_->adjacentLoops[s_index]);
				temp = start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop);
				assert(temp > 0);

			}
		}
		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < end_->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'S' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'S')) {

		StackLoop *start_;
		OpenLoop *end_;
		Loop *newLoop;
		int s_index = 0;
		int e_index = 0;
		int temp = 0;

		if (start->identity == 'S') {

			start_ = (StackLoop *) start;
			end_ = (OpenLoop *) end;

		} else {

			start_ = (StackLoop *) end;
			end_ = (OpenLoop *) start;

		}

		for (int loop = 0; (loop < end_->numAdjacent) || (loop < 2); loop++) {

			if (loop <= 1 && start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}

			if (loop < end_->numAdjacent && end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}

		}

//		for (int loop = 0; (loop < end_->numAdjacent) || (loop < 2); loop++) {
//
//			if (loop <= 1 && start_->adjacentLoops[loop] != end_) {
//				s_index = loop;
//			}
//
//			if (loop < end_->numAdjacent && end_->adjacentLoops[loop] == start_) {
//				e_index = loop;
//			}
//
//		}

		// note e_index has different meaning now for openloops.

		cout << "s_index, e_index = " << s_index << "  " << e_index << endl;

		cout << "Start is " << endl;
		cout << (start_)->typeInternalsToString() << endl;

		cout << "End is " << endl;
		cout << end_->typeInternalsToString() << endl;

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent + 1];
		char **seqs = new char *[end_->numAdjacent + 1];

//		bool settedPair = false;

		for (int loop = 0; loop < end_->numAdjacent + 1; loop++) {
			if (loop == e_index) {
				pairtypes[loop] = start_->pairtype[s_index];
//				settedPair = true;
				sidelens[loop] = end_->sidelen[loop] + 1;
				seqs[loop] = end_->seqs[loop];
			} else if (loop == e_index + 1) {
				if (loop < end_->numAdjacent){
					pairtypes[loop] = end_->pairtype[loop];
//					settedPair = true;
				}
				sidelens[loop] = end_->sidelen[loop] + 1;
				seqs[loop] = start_->seqs[s_index];
			} else {
				if (loop < end_->numAdjacent){
					pairtypes[loop] = end_->pairtype[loop];
//					settedPair = true;
				}
				sidelens[loop] = end_->sidelen[loop];
				seqs[loop] = end_->seqs[loop];
			}
		}

//		cout << "SettedPair is: " << settedPair << endl;

		cout << "new Pairtype is CUSTOM: " << pairtypes[0] << endl;

		// resulting will be an open loop, same# of adjacent helices, two sides longer by one base, and one pairtype possibly changed.
		newLoop = new OpenLoop(end_->numAdjacent, pairtypes, sidelens, seqs);

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				newLoop->addAdjacent(end_->adjacentLoops[loop]);
				temp = end_->adjacentLoops[loop]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);
			} else {
				newLoop->addAdjacent(start_->adjacentLoops[s_index]);
				temp = start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop);
				assert(temp > 0);

			}
		}
		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < end_->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if (start->identity == 'I' && end->identity == 'I') {
		InteriorLoop *start_, *end_;
		Loop *start_extra, *end_extra, *newLoop;
		int s_index = 0, e_index = 0;
		start_ = (InteriorLoop *) start;
		end_ = (InteriorLoop *) end;

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				start_extra = start_->adjacentLoops[loop];
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				end_extra = end_->adjacentLoops[loop];
				e_index = loop;
			}
		}

		// resulting will be an interior loop side lengths equal to the length
		// of the 'input' interior loop+1.
		if (start_->int_seq[s_index] < end_->int_seq[e_index]) {
			newLoop = new InteriorLoop(start_->pairtype[s_index], end_->pairtype[e_index], end_->sizes[1 - e_index] + start_->sizes[s_index] + 1,
					end_->sizes[e_index] + start_->sizes[1 - s_index] + 1, start_->int_seq[s_index], end_->int_seq[e_index]);

			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);
			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
		} else {
			newLoop = new InteriorLoop(end_->pairtype[e_index], start_->pairtype[s_index], end_->sizes[e_index] + start_->sizes[1 - s_index] + 1,
					end_->sizes[1 - e_index] + start_->sizes[s_index] + 1, end_->int_seq[e_index], start_->int_seq[s_index]);

			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);

		}
		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();
		end_->adjacentLoops[e_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'I' && end->identity == 'B') || (start->identity == 'B' && end->identity == 'I')) {
		InteriorLoop *start_;
		BulgeLoop *end_;
		Loop *newLoop;
		int s_index = 0, e_index = 0;

		if (start->identity == 'S') {
			start_ = (InteriorLoop *) start;
			end_ = (BulgeLoop *) end;
		} else {
			start_ = (InteriorLoop *) end;
			end_ = (BulgeLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				e_index = loop;
			}
		}
		// resulting will be an interior loop side lengths equal to the length
		// of the 'input' bulge loop, plus one on each side (ie, one side will be B+1, the other 1.

		if (start_->int_seq[s_index] < end_->bulge_seq[e_index]) {
			newLoop = new InteriorLoop(start_->pairtype[s_index], end_->pairtype[e_index], end_->bulgesize[1 - e_index] + 1 + start_->sizes[s_index],
					end_->bulgesize[e_index] + 1 + start_->sizes[1 - s_index], start_->int_seq[s_index], end_->bulge_seq[e_index]);

			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);
			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
		} else {
			newLoop = new InteriorLoop(end_->pairtype[e_index], start_->pairtype[s_index], end_->bulgesize[e_index] + 1 + start_->sizes[1 - s_index],
					end_->bulgesize[1 - e_index] + 1 + start_->sizes[s_index], end_->bulge_seq[e_index], start_->int_seq[s_index]);

			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);

		}
		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();
		end_->adjacentLoops[e_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'I' && end->identity == 'H') || (start->identity == 'H' && end->identity == 'I')) {
		InteriorLoop *start_;
		HairpinLoop *end_;
		int s_index = 0;
		Loop *newLoop;

		if (start->identity == 'I') {
			start_ = (InteriorLoop *) start;
			end_ = (HairpinLoop *) end;
		} else {
			start_ = (InteriorLoop *) end;
			end_ = (HairpinLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
		}
		// end is the hairpin, which has no extra adjacencies.

		// resulting will be a hairpin loop with previous size, plus interior loop's sizes (both) plus 2 (for the pairing that's now unpaired)
		newLoop = new HairpinLoop(start_->pairtype[s_index], end_->hairpinsize + 2 + start_->sizes[0] + start_->sizes[1], start_->int_seq[s_index]);

		newLoop->addAdjacent(start_->adjacentLoops[s_index]);
		// TODO: fix this! asserts generate no code when NDEBUG is set!
		assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);

		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	else if ((start->identity == 'I' && end->identity == 'M') || (start->identity == 'M' && end->identity == 'I')) {
		InteriorLoop *start_;
		MultiLoop *end_;
		Loop *newLoop;
		int s_index = 0, e_index = 0, temp = 0;

		if (start->identity == 'I') {
			start_ = (InteriorLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (InteriorLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent];
		char **seqs = new char *[end_->numAdjacent];

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				pairtypes[loop] = end_->pairtype[loop];
				if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
					sidelens[loop] = end_->sidelen[loop];
				else
					sidelens[loop] = end_->sidelen[loop] + start_->sizes[1 - s_index] + 1;
				seqs[loop] = end_->seqs[loop];
			} else {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + start_->sizes[s_index] + 1;
				seqs[loop] = start_->int_seq[s_index];
			}
		}

		// resulting will be an multiloop, same# of adjacent helices, two sides longer by one plus interior loop size, and one pairtype possibly changed.
		newLoop = new MultiLoop(end_->numAdjacent, pairtypes, sidelens, seqs);

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				newLoop->addAdjacent(end_->adjacentLoops[loop]);
				temp = end_->adjacentLoops[loop]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);
			} else {
				newLoop->addAdjacent(start_->adjacentLoops[s_index]);
				temp = start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop);
				assert(temp > 0);

			}
		}
		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < end_->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'I' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'I')) {
		InteriorLoop *start_;
		OpenLoop *end_;
		Loop *newLoop;
		int s_index = 0, e_index = 0, temp = 0;

		if (start->identity == 'I') {
			start_ = (InteriorLoop *) start;
			end_ = (OpenLoop *) end;
		} else {
			start_ = (InteriorLoop *) end;
			end_ = (OpenLoop *) start;
		}

		for (int loop = 0; (loop < end_->numAdjacent) || (loop < 2); loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (loop < end_->numAdjacent)
				if (end_->adjacentLoops[loop] == start_) {
					e_index = loop;
				}
		}
		// note e_index has different meaning now for openloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent + 1];
		char **seqs = new char *[end_->numAdjacent + 1];

		for (int loop = 0; loop < end_->numAdjacent + 1; loop++) {
			if (loop == e_index) {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + 1 + start_->sizes[1 - s_index];
				seqs[loop] = end_->seqs[loop];
			} else if (loop == e_index + 1) {
				if (loop < end_->numAdjacent)
					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop] + 1 + start_->sizes[s_index];
				seqs[loop] = start_->int_seq[s_index];
			} else {
				if (loop < end_->numAdjacent)
					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop];
				seqs[loop] = end_->seqs[loop];
			}
		}

		// resulting will be an open loop, same# of adjacent helices, two sides longer by one base, and one pairtype possibly changed.
		newLoop = new OpenLoop(end_->numAdjacent, pairtypes, sidelens, seqs);

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				newLoop->addAdjacent(end_->adjacentLoops[loop]);
				temp = end_->adjacentLoops[loop]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);
			} else {
				newLoop->addAdjacent(start_->adjacentLoops[s_index]);
				temp = start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop);
				assert(temp > 0);

			}
		}
		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < end_->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

// end interior
// begin bulge

	if (start->identity == 'B' && end->identity == 'B') {
		BulgeLoop *start_, *end_;
		Loop *start_extra, *end_extra, *newLoop;
		int s_index = 0, e_index = 0;
		start_ = (BulgeLoop *) start;
		end_ = (BulgeLoop *) end;

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				start_extra = start_->adjacentLoops[loop];
				s_index = loop;
			}
			if (end_->adjacentLoops[loop] != start_) {
				end_extra = end_->adjacentLoops[loop];
				e_index = loop;
			}
		}

		// resulting will be an interior loop side lengths equal to the length
		// of the 'input' interior loop+1.
		if (start_->bulge_seq[s_index] < end_->bulge_seq[e_index]) {
			newLoop = new InteriorLoop(start_->pairtype[s_index], end_->pairtype[e_index], end_->bulgesize[1 - e_index] + start_->bulgesize[s_index] + 1,
					end_->bulgesize[e_index] + start_->bulgesize[1 - s_index] + 1, start_->bulge_seq[s_index], end_->bulge_seq[e_index]);

			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);
			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
		} else {
			newLoop = new InteriorLoop(end_->pairtype[e_index], start_->pairtype[s_index], end_->bulgesize[e_index] + start_->bulgesize[1 - s_index] + 1,
					end_->bulgesize[1 - e_index] + start_->bulgesize[s_index] + 1, end_->bulge_seq[e_index], start_->bulge_seq[s_index]);

			newLoop->addAdjacent(end_->adjacentLoops[e_index]);
			// TODO: fix this too! see above comment.
			assert(end_->adjacentLoops[e_index]->replaceAdjacent(end_, newLoop) > 0);
			newLoop->addAdjacent(start_->adjacentLoops[s_index]);
			// TODO: fix this! asserts generate no code when NDEBUG is set!
			assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);

		}
		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();
		end_->adjacentLoops[e_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'B' && end->identity == 'H') || (start->identity == 'H' && end->identity == 'B')) {
		BulgeLoop *start_;
		HairpinLoop *end_;
		int s_index = 0;
		Loop *newLoop;

		if (start->identity == 'B') {
			start_ = (BulgeLoop *) start;
			end_ = (HairpinLoop *) end;
		} else {
			start_ = (BulgeLoop *) end;
			end_ = (HairpinLoop *) start;
		}

		for (int loop = 0; loop <= 1; loop++) {
			if (start_->adjacentLoops[loop] != end_) {
				s_index = loop;
			}
		}
		// end is the hairpin, which has no extra adjacencies.

		// resulting will be a hairpin loop with previous size, plus interior loop's sizes (both) plus 2 (for the pairing that's now unpaired)
		newLoop = new HairpinLoop(start_->pairtype[s_index], end_->hairpinsize + 2 + start_->bulgesize[0] + start_->bulgesize[1], start_->bulge_seq[s_index]);

		newLoop->addAdjacent(start_->adjacentLoops[s_index]);
		// TODO: fix this! asserts generate no code when NDEBUG is set!
		assert(start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop) > 0);

		newLoop->generateMoves();

		// need to re-generate the moves for the two adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		start_->adjacentLoops[s_index]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	else if ((start->identity == 'B' && end->identity == 'M') || (start->identity == 'M' && end->identity == 'B')) {
		BulgeLoop *start_;
		MultiLoop *end_;
		Loop *newLoop;
		int s_index = 0, e_index = 0, temp = 0;

		if (start->identity == 'B') {
			start_ = (BulgeLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (BulgeLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent];
		char **seqs = new char *[end_->numAdjacent];

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				pairtypes[loop] = end_->pairtype[loop];
				if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
					sidelens[loop] = end_->sidelen[loop];
				else
					sidelens[loop] = end_->sidelen[loop] + start_->bulgesize[1 - s_index] + 1;
				seqs[loop] = end_->seqs[loop];
			} else {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + start_->bulgesize[s_index] + 1;
				seqs[loop] = start_->bulge_seq[s_index];
			}
		}

		// resulting will be an multiloop, same# of adjacent helices, two sides longer by one plus interior loop size, and one pairtype possibly changed.
		newLoop = new MultiLoop(end_->numAdjacent, pairtypes, sidelens, seqs);

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				newLoop->addAdjacent(end_->adjacentLoops[loop]);
				temp = end_->adjacentLoops[loop]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);
			} else {
				newLoop->addAdjacent(start_->adjacentLoops[s_index]);
				temp = start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop);
				assert(temp > 0);

			}
		}
		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < end_->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'B' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'B')) {
		BulgeLoop *start_;
		OpenLoop *end_;
		Loop *newLoop;
		int s_index = 0, e_index = 0, temp = 0;

		if (start->identity == 'B') {
			start_ = (BulgeLoop *) start;
			end_ = (OpenLoop *) end;
		} else {
			start_ = (BulgeLoop *) end;
			end_ = (OpenLoop *) start;
		}

		for (int loop = 0; (loop < end_->numAdjacent) || (loop < 2); loop++) {
			if (loop <= 1)
				if (start_->adjacentLoops[loop] != end_) {
					s_index = loop;
				}
			if (loop < end_->numAdjacent)
				if (end_->adjacentLoops[loop] == start_) {
					e_index = loop;
				}
		}
		// note e_index has different meaning now for openloops.

		int *pairtypes = new int[end_->numAdjacent];
		int *sidelens = new int[end_->numAdjacent + 1];
		char **seqs = new char *[end_->numAdjacent + 1];

		for (int loop = 0; loop < end_->numAdjacent + 1; loop++) {
			if (loop == e_index) {
				pairtypes[loop] = start_->pairtype[s_index];
				sidelens[loop] = end_->sidelen[loop] + 1 + start_->bulgesize[1 - s_index];
				seqs[loop] = end_->seqs[loop];
			} else if (loop == e_index + 1) {
				if (loop < end_->numAdjacent)
					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop] + 1 + start_->bulgesize[s_index];
				seqs[loop] = start_->bulge_seq[s_index];
			} else {
				if (loop < end_->numAdjacent)
					pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop];
				seqs[loop] = end_->seqs[loop];
			}
		}

		// resulting will be an open loop, same# of adjacent helices, two sides longer by one base, and one pairtype possibly changed.
		newLoop = new OpenLoop(end_->numAdjacent, pairtypes, sidelens, seqs);

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				newLoop->addAdjacent(end_->adjacentLoops[loop]);
				temp = end_->adjacentLoops[loop]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);
			} else {
				newLoop->addAdjacent(start_->adjacentLoops[s_index]);
				temp = start_->adjacentLoops[s_index]->replaceAdjacent(start_, newLoop);
				assert(temp > 0);

			}
		}
		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < end_->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

// end bulge

// start hairpin

	if (start->identity == 'H' && end->identity == 'H') {
		fprintf(stderr, "Hairpin/Hairpin deletion move encountered - not currently supported.\n");
		assert(0);
		return NULL;
	}

	if ((start->identity == 'H' && end->identity == 'M') || (start->identity == 'M' && end->identity == 'H')) {
		HairpinLoop *start_;
		MultiLoop *end_;
		Loop *newLoop;
		int e_index = 0, temp;

		if (start->identity == 'H') {
			start_ = (HairpinLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (HairpinLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for multiloops.

		// Several cases we now need to check for:
		// 1. Multiloop had 3 adjacent stacks, and the resulting figure with two adjacent stacks is now a interior loop.
		// 2. Multiloop had 3 adjacent stacks, and the resulting structure with two adjacent stacks is now a bulge loop
		// 3. Multiloop had >3 adjacent stacks, and the resulting structure is a multiloop with one less adjacent stack.

		if (end_->numAdjacent == 3) {
			int sizes[2];
			int positions[2] = { (e_index + 1) % 3, (e_index + 2) % 3 };
			if (e_index == 0)
				sizes[0] = end_->sidelen[0] + end_->sidelen[2] + start_->hairpinsize + 2;
			else
				sizes[0] = end_->sidelen[e_index] + end_->sidelen[e_index - 1] + start_->hairpinsize + 2;

			sizes[1] = end_->sidelen[(e_index + 1) % 3];

			if (sizes[1] > 0) // interior loop case
					{
				// TODO: check for other idiot uses of sequence compares
				// and see if they caused as much problems as this one
				// (false case was putting the adjacencies back in the wrong order)
				/*
				 if( end_->seqs[positions[0]] > end_->seqs[positions[1]] )
				 {*/
				newLoop = new InteriorLoop(end_->pairtype[positions[1]], end_->pairtype[positions[0]], sizes[0], sizes[1], end_->seqs[positions[1]],
						end_->seqs[positions[0]]);

				newLoop->addAdjacent(end_->adjacentLoops[positions[1]]);
				temp = end_->adjacentLoops[positions[1]]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);

				newLoop->addAdjacent(end_->adjacentLoops[positions[0]]);
				temp = end_->adjacentLoops[positions[0]]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);

				newLoop->generateMoves();

				// need to re-generate the moves for all adjacent loops.
				// TODO: change this to only re-generate the deletion moves.
				newLoop->adjacentLoops[0]->generateMoves();
				newLoop->adjacentLoops[1]->generateMoves();

				start_->cleanupAdjacent();
				delete start_;
				end_->cleanupAdjacent();
				delete end_;
				return newLoop;

			} else  // bulge loop case
			{
				newLoop = new BulgeLoop(end_->pairtype[positions[0]], end_->pairtype[positions[1]], sizes[1], sizes[0], end_->seqs[positions[0]],
						end_->seqs[positions[1]]);

				newLoop->addAdjacent(end_->adjacentLoops[positions[0]]);
				temp = end_->adjacentLoops[positions[0]]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);

				newLoop->addAdjacent(end_->adjacentLoops[positions[1]]);
				temp = end_->adjacentLoops[positions[1]]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);

				newLoop->generateMoves();

				// need to re-generate the moves for all adjacent loops.
				// TODO: change this to only re-generate the deletion moves.
				newLoop->adjacentLoops[0]->generateMoves();
				newLoop->adjacentLoops[1]->generateMoves();
				start_->cleanupAdjacent();
				delete start_;
				end_->cleanupAdjacent();
				delete end_;
				return newLoop;
			}
		}

		else if (end_->numAdjacent > 3)              // multiloop case
				{
			int *pairtypes = new int[end_->numAdjacent - 1];
			int *sidelens = new int[end_->numAdjacent - 1];
			char **seqs = new char *[end_->numAdjacent - 1];

			for (int loop = 0; loop < end_->numAdjacent; loop++) {
				if (loop != e_index) {
					if (loop < e_index) {
						pairtypes[loop] = end_->pairtype[loop];
						if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
							sidelens[loop] = end_->sidelen[loop];
						else
							sidelens[loop] = end_->sidelen[loop] + end_->sidelen[e_index] + 2 + start_->hairpinsize;
						seqs[loop] = end_->seqs[loop];
					}
					if (loop > e_index) {
						pairtypes[loop - 1] = end_->pairtype[loop];
						if ((loop != e_index - 1 && e_index != 0) || (loop != end_->numAdjacent - 1 && e_index == 0))
							sidelens[loop - 1] = end_->sidelen[loop];
						else
							sidelens[loop - 1] = end_->sidelen[loop] + end_->sidelen[e_index] + 2 + start_->hairpinsize;
						seqs[loop - 1] = end_->seqs[loop];

					}
				}
			}

			newLoop = new MultiLoop(end_->numAdjacent - 1, pairtypes, sidelens, seqs);
			for (int loop = 0; loop < end_->numAdjacent; loop++) {
				if (loop != e_index) {
					newLoop->addAdjacent(end_->adjacentLoops[loop]);
					temp = end_->adjacentLoops[loop]->replaceAdjacent(end_, newLoop);
					assert(temp > 0);
				}
			}

			newLoop->generateMoves();

			// need to re-generate the moves for all adjacent loops.
			// TODO: change this to only re-generate the deletion moves.
			for (int loop = 0; loop < newLoop->numAdjacent; loop++)
				newLoop->adjacentLoops[loop]->generateMoves();

			start_->cleanupAdjacent();
			delete start_;
			end_->cleanupAdjacent();
			delete end_;
			return newLoop;
		}

	}

	if ((start->identity == 'H' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'H')) {
		HairpinLoop *start_;
		OpenLoop *end_;
		Loop *newLoop;
		int e_index = 0, temp;

		if (start->identity == 'H') {
			start_ = (HairpinLoop *) start;
			end_ = (OpenLoop *) end;
		} else {
			start_ = (HairpinLoop *) end;
			end_ = (OpenLoop *) start;
		}

		for (int loop = 0; (loop < end_->numAdjacent); loop++) {
			if (end_->adjacentLoops[loop] == start_) {
				e_index = loop;
			}
		}
		// note e_index has different meaning now for openloops.

		int *pairtypes = new int[end_->numAdjacent - 1];
		int *sidelens = new int[end_->numAdjacent];
		char **seqs = new char *[end_->numAdjacent];

		for (int loop = 0; loop < end_->numAdjacent + 1; loop++) {
			if (loop < e_index) {
				pairtypes[loop] = end_->pairtype[loop];
				sidelens[loop] = end_->sidelen[loop];
				seqs[loop] = end_->seqs[loop];
			} else if (loop == e_index) {
				if (loop < end_->numAdjacent - 1)
					pairtypes[loop] = end_->pairtype[loop + 1];
				sidelens[loop] = end_->sidelen[loop] + 2 + end_->sidelen[loop + 1] + start_->hairpinsize;
				seqs[loop] = end_->seqs[loop];
			} else if (loop > e_index + 1) {
				if (loop < end_->numAdjacent)
					pairtypes[loop - 1] = end_->pairtype[loop];
				sidelens[loop - 1] = end_->sidelen[loop];
				seqs[loop - 1] = end_->seqs[loop];
			}
		}

		newLoop = new OpenLoop(end_->numAdjacent - 1, pairtypes, sidelens, seqs);
		for (int loop = 0; loop < end_->numAdjacent; loop++) {
			if (loop != e_index) {
				newLoop->addAdjacent(end_->adjacentLoops[loop]);
				temp = end_->adjacentLoops[loop]->replaceAdjacent(end_, newLoop);
				assert(temp > 0);
			}
		}

		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < newLoop->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

// end hairpin

// start multiloop

	if (start->identity == 'M' && end->identity == 'M') {
		MultiLoop *start_;
		MultiLoop *end_;
		Loop *newLoop;
		int s_index = 0, e_index = 0, index = 0, temp;

		start_ = (MultiLoop *) start;
		end_ = (MultiLoop *) end;

		for (int loop = 0; loop < end_->numAdjacent || loop < start_->numAdjacent; loop++) {
			if (loop < start_->numAdjacent)
				if (start_->adjacentLoops[loop] == end_) {
					s_index = loop;
				}
			if (loop < end_->numAdjacent)
				if (end_->adjacentLoops[loop] == start_) {
					e_index = loop;
				}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent + start_->numAdjacent - 2];
		int *sidelens = new int[end_->numAdjacent + start_->numAdjacent - 2];
		char **seqs = new char *[end_->numAdjacent + start_->numAdjacent - 2];

		index = 0;
		for (int loop = 0; loop < start_->numAdjacent; loop++) {
			if (loop != s_index) {
				pairtypes[index] = start_->pairtype[loop];
				if ((loop != s_index - 1 && s_index != 0) || (loop != start_->numAdjacent - 1 && s_index == 0))
					sidelens[index] = start_->sidelen[loop];
				else
					sidelens[index] = start_->sidelen[loop] + end_->sidelen[e_index] + 1;
				seqs[index] = start_->seqs[loop];
				index++;
			} else {
				for (int loop2 = 1; loop2 < end_->numAdjacent; loop2++) {
					temp = (e_index + loop2) % end_->numAdjacent;

					pairtypes[index] = end_->pairtype[temp];
					if ((temp != e_index - 1 && e_index != 0) || (temp != end_->numAdjacent - 1 && e_index == 0))
						sidelens[index] = end_->sidelen[temp];
					else
						sidelens[index] = end_->sidelen[temp] + start_->sidelen[s_index] + 1;
					seqs[index] = end_->seqs[temp];
					index++;
				}
			}
		}
		assert(index == end_->numAdjacent + start_->numAdjacent - 2);

		newLoop = new MultiLoop(start_->numAdjacent + end_->numAdjacent - 2, pairtypes, sidelens, seqs);
		for (int loop = 0; loop < start_->numAdjacent; loop++) {
			if (loop != s_index) {
				newLoop->addAdjacent(start_->adjacentLoops[loop]);
				temp = start_->adjacentLoops[loop]->replaceAdjacent(start_, newLoop);
				assert(temp > 0);
			} else {
				for (int loop2 = 1; loop2 < end_->numAdjacent; loop2++) {
					int temp2 = (e_index + loop2) % end_->numAdjacent;
					newLoop->addAdjacent(end_->adjacentLoops[temp2]);
					temp = end_->adjacentLoops[temp2]->replaceAdjacent(end_, newLoop);
					assert(temp > 0);
				}
			}
		}

		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < newLoop->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

	if ((start->identity == 'M' && end->identity == 'O') || (start->identity == 'O' && end->identity == 'M')) {
		Loop *newLoop;
		OpenLoop *start_;
		MultiLoop *end_;
		int s_index = 0, e_index = 0, index = 0, temp;

		start_ = (OpenLoop *) start;
		end_ = (MultiLoop *) end;

		if (start->identity == 'O') {
			start_ = (OpenLoop *) start;
			end_ = (MultiLoop *) end;
		} else {
			start_ = (OpenLoop *) end;
			end_ = (MultiLoop *) start;
		}

		for (int loop = 0; loop < end_->numAdjacent || loop < start_->numAdjacent; loop++) {
			if (loop < start_->numAdjacent)
				if (start_->adjacentLoops[loop] == end_) {
					s_index = loop;
				}
			if (loop < end_->numAdjacent)
				if (end_->adjacentLoops[loop] == start_) {
					e_index = loop;
				}
		}
		// note e_index has different meaning now for multiloops.

		int *pairtypes = new int[end_->numAdjacent + start_->numAdjacent - 2];
		int *sidelens = new int[end_->numAdjacent + start_->numAdjacent - 1];
		char **seqs = new char *[end_->numAdjacent + start_->numAdjacent - 1];

		index = 0;
		for (int loop = 0; loop <= start_->numAdjacent; loop++) {
			if (loop == s_index + 1) {
				temp = (e_index + end_->numAdjacent - 1) % end_->numAdjacent;

				if (loop < start_->numAdjacent)
					pairtypes[index] = start_->pairtype[loop];
				sidelens[index] = start_->sidelen[loop] + end_->sidelen[temp] + 1;
				seqs[index] = end_->seqs[temp];
				if (loop < start_->numAdjacent)
					index++;
			} else if (loop != s_index) {
				if (loop < start_->numAdjacent)
					pairtypes[index] = start_->pairtype[loop];
				sidelens[index] = start_->sidelen[loop];
				seqs[index] = start_->seqs[loop];
				if (loop < start_->numAdjacent)
					index++;
			} else {
				for (int loop2 = 1; loop2 < end_->numAdjacent; loop2++) {
					temp = (e_index + loop2) % end_->numAdjacent;
					int temp2 = (e_index + loop2 - 1) % end_->numAdjacent;

					pairtypes[index] = end_->pairtype[temp];
					if (loop2 == 1) {
						sidelens[index] = end_->sidelen[e_index] + start_->sidelen[s_index] + 1;
						seqs[index] = start_->seqs[s_index];
					} else {
						sidelens[index] = end_->sidelen[temp2];
						seqs[index] = end_->seqs[temp2];
					}
					index++;
				}
			}
		}
		assert(index == end_->numAdjacent + start_->numAdjacent - 2);

		newLoop = new OpenLoop(start_->numAdjacent + end_->numAdjacent - 2, pairtypes, sidelens, seqs);
		for (int loop = 0; loop < start_->numAdjacent; loop++) {
			if (loop != s_index) {
				newLoop->addAdjacent(start_->adjacentLoops[loop]);
				temp = start_->adjacentLoops[loop]->replaceAdjacent(start_, newLoop);
				assert(temp > 0);
			} else {
				for (int loop2 = 1; loop2 < end_->numAdjacent; loop2++) {
					int temp2 = (e_index + loop2) % end_->numAdjacent;
					newLoop->addAdjacent(end_->adjacentLoops[temp2]);
					temp = end_->adjacentLoops[temp2]->replaceAdjacent(end_, newLoop);
					assert(temp > 0);
				}
			}
		}

		newLoop->generateMoves();

		// need to re-generate the moves for all adjacent loops.
		// TODO: change this to only re-generate the deletion moves.
		for (int loop = 0; loop < newLoop->numAdjacent; loop++)
			newLoop->adjacentLoops[loop]->generateMoves();

		start_->cleanupAdjacent();
		delete start_;
		end_->cleanupAdjacent();
		delete end_;
		return newLoop;
	}

// end multiloop

// start openloop

// Control flow should never reach here, as Scomplex shortcuts O/O deletion moves (complex breaks) to a different function - performComplexSplit
	if (start->identity == 'O' && end->identity == 'O') {
		fprintf(stderr, "Openloop/Openloop deletion reached via performDeleteMove, bad control flow\n");
		assert(0);
		return NULL;
	}

// end openloop

	else
		return NULL;
}

/* StackLoop */

void StackLoop::calculateEnergy(void) {

	assert(Loop::energyModel != NULL);

//	if (Loop::energyModel == NULL){
//		return; // we can't handle this error. I'm trying to work out a way around it, but generally if the loops try to get used before the energy model initializes, it's all over.
//	}

	energy = Loop::energyModel->StackEnergy(seqs[0][0], seqs[1][1], seqs[0][1], seqs[1][0]);

	return;
}

void StackLoop::generateMoves(void) {
	generateDeleteMoves();
}

void StackLoop::generateDeleteMoves(void) {
	double temprate;
	if (moves != NULL)
		delete moves;
	moves = new MoveList(0); // always have 2 delete moves, no shift moves and no creation moves.

	generateAndSaveDeleteMove(adjacentLoops[0], 0);
	generateAndSaveDeleteMove(adjacentLoops[1], 1);

	totalRate = moves->getRate();
}

void StackLoop::printMove(Loop *comefrom, char *structure_p, char *seq_p) {
	int loop;
	structure_p[seqs[0] - seq_p] = (seqs[0] < seqs[1]) ? '(' : ')';
	structure_p[seqs[1] - seq_p + 1] = (seqs[0] < seqs[1]) ? ')' : '(';
	structure_p[seqs[0] - seq_p + 1] = (seqs[0] < seqs[1]) ? '(' : ')';
	structure_p[seqs[1] - seq_p] = (seqs[0] < seqs[1]) ? ')' : '(';

	for (loop = 0; loop < curAdjacent; loop++) {
		if (adjacentLoops[loop] != comefrom && adjacentLoops[loop] != NULL)
			// shouldn't happen, being careful.
			adjacentLoops[loop]->printMove(this, structure_p, seq_p);
		assert(adjacentLoops[loop] != NULL);
	}
// for now, Stack loops also don't have deletion moves.
}

Move *StackLoop::getChoice(double *randomchoice, Loop *from) {
	Move *stor;
	assert(randomchoice != NULL);
	assert(*randomchoice >= 0.0); // never should see a negative choice value.

	if (*randomchoice < totalRate) // something was chosen, do this
		return moves->getChoice(randomchoice);
	else {
		*randomchoice -= totalRate;
		if (adjacentLoops[0] != from) {
			stor = adjacentLoops[0]->getChoice(randomchoice, this);
			if (stor != NULL)
				return stor;
		}
		if (adjacentLoops[1] != from) {
			stor = adjacentLoops[1]->getChoice(randomchoice, this);
			if (stor != NULL)
				return stor;
		}
	}
	return NULL;
}

double StackLoop::doChoice(Move *move, Loop **returnLoop) {
	;
// currently no moves in stackloop, in the future will include deletion moves
}

char *StackLoop::getLocation(Move *move, int index) {
	if (move->getType() & MOVE_CREATE) // then something's wrong!
		assert(0);

	else if (move->getType() & MOVE_DELETE) {
		if (adjacentLoops[0] == move->affected[0] || adjacentLoops[0] == move->affected[1])
			return seqs[0];
		else if (adjacentLoops[1] == move->affected[0] || adjacentLoops[1] == move->affected[1])
			return seqs[1];
	}
	assert(0);
}

char *StackLoop::verifyLoop(char *incoming_sequence, int incoming_pairtype, Loop *from) {
	char *ret_seq;
	if (adjacentLoops[0] == from) {
		if (incoming_sequence != seqs[0]) {
			fprintf(stderr, "Verification Failed\n");
			assert(incoming_sequence == seqs[0]);
		}
		ret_seq = adjacentLoops[1]->verifyLoop(&seqs[0][1], pairtype[1], this);
		assert(ret_seq == seqs[1]);
		return seqs[1] + 1;
	} else if (adjacentLoops[1] == from) {
		if (incoming_sequence != seqs[1]) {
			fprintf(stderr, "Verification Failed\n");
			assert(incoming_sequence == seqs[1]);
		}
		ret_seq = adjacentLoops[0]->verifyLoop(seqs[1] + 1, pairtype[0], this);
		assert(ret_seq == seqs[0]);
		return seqs[0] + 1;
	} else
		assert(0);
	return NULL;
}

StackLoop::StackLoop(void) {
	numAdjacent = 2;
	adjacentLoops = new Loop *[2];
	identity = 'S';
}



StackLoop::StackLoop(int type1, int type2, char *seq1, char *seq2, Loop *left, Loop *right) // left and right default to NULL, see header.
		{

//	pairtype[0] = type1;
//	pairtype[1] = type2;
	numAdjacent = 2;
	adjacentLoops = new Loop *[2];
	adjacentLoops[0] = left;
	adjacentLoops[1] = right;
	curAdjacent = (left == NULL ? 0 : 1) + (right == NULL ? 0 : 1);
	identity = 'S';
	seqs[0] = seq1;
	seqs[1] = seq2;
	pairtype[0] = seqs[0][0];
	pairtype[1] = seqs[0][1];

}

string StackLoop::typeInternalsToString(void) {

	std::stringstream ss;

	ss << "(" << baseTypeString[seqs[0][0]] << "/";
	ss << baseTypeString[seqs[1][1]] << ", ";

	ss << baseTypeString[seqs[0][1]] << "/";
	ss << baseTypeString[seqs[1][0]];

	ss << ") \n";

	ss << "  --   ";
	ss << basepairString[pairtype[0]] << ",   ";
	ss << basepairString[pairtype[1]] << endl;


	return ss.str();

}

/*
 HairpinLoop Functions
 */

HairpinLoop::HairpinLoop(void) {
	numAdjacent = 1;
	adjacentLoops = new Loop *[1];

	hairpinsize = 0;
	hairpin_seq = NULL;
	identity = 'H';
}

HairpinLoop::HairpinLoop(int type, int size, char *hairpin_sequence, Loop *previous) {
	numAdjacent = 1;
	adjacentLoops = new Loop *[1];
	adjacentLoops[0] = previous;
	if (previous != NULL)
		curAdjacent = 1;

	pairtype = type;
	hairpinsize = size;
	hairpin_seq = hairpin_sequence;
	identity = 'H';
}

string HairpinLoop::typeInternalsToString(void) {

	std::stringstream ss;

	ss << "pairType =" << pairtype << ", ";
	ss << "seq= " << utility::sequenceToString(hairpin_seq, hairpinsize) << "\n";

	return ss.str();

}

void HairpinLoop::calculateEnergy(void) {
	if (energyModel == NULL)
		return; // we can't handle this error. I'm trying to work out a way around it, but generally if the loops try to get used before the energy model initializes, it's all over.

	energy = energyModel->HairpinEnergy(hairpin_seq, hairpinsize);
	return;
}

Move *HairpinLoop::getChoice(double *randomchoice, Loop *from) {
	Move *stor;
	assert(randomchoice != NULL);
	assert(*randomchoice >= 0.0); // never should see a negative choice value.

	if (*randomchoice < totalRate) // something was chosen, do this
		return moves->getChoice(randomchoice);
	else {
		*randomchoice -= totalRate;
		if (adjacentLoops[0] != from) // should not occur, can remove later. FIXME
				{
			stor = adjacentLoops[0]->getChoice(randomchoice, this);
			if (stor != NULL)
				return stor;
		}
	}
	return NULL;
}

double HairpinLoop::doChoice(Move *move, Loop **returnLoop) {
	Loop *newLoop[2];
	int pt, loop, loop2;
	if (move->type & MOVE_CREATE) {
		loop = move->index[0];
		loop2 = move->index[1];
		pt = pairtypes[hairpin_seq[loop]][hairpin_seq[loop2]];
		if (move->type & MOVE_1) // stack and hairpin
				{
			newLoop[0] = new StackLoop(pairtype, pt, hairpin_seq, &hairpin_seq[loop2]);
			newLoop[1] = new HairpinLoop(pt, hairpinsize - 2, &hairpin_seq[1]);
			newLoop[0]->addAdjacent(adjacentLoops[0]);
			adjacentLoops[0]->replaceAdjacent(this, newLoop[0]);
			newLoop[0]->addAdjacent(newLoop[1]);
			newLoop[1]->addAdjacent(newLoop[0]);
			adjacentLoops[0]->generateMoves();
			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}
		if (move->type & MOVE_2) // bulge and hairpin
				{
			if (loop == 1)
				newLoop[0] = new BulgeLoop(pairtype, pt, 0, (hairpinsize - loop2), &hairpin_seq[0], &hairpin_seq[loop2]);
			else
				newLoop[0] = new BulgeLoop(pairtype, pt, loop - 1, 0, &hairpin_seq[0], &hairpin_seq[loop2]);
			newLoop[1] = new HairpinLoop(pt, loop2 - loop - 1, &hairpin_seq[loop]);
			newLoop[0]->addAdjacent(adjacentLoops[0]);
			adjacentLoops[0]->replaceAdjacent(this, newLoop[0]);
			newLoop[0]->addAdjacent(newLoop[1]);
			newLoop[1]->addAdjacent(newLoop[0]);
			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			adjacentLoops[0]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}
		if (move->type & MOVE_3) // interior and hairpin
				{
			newLoop[0] = new InteriorLoop(pairtype, pt, loop - 1, hairpinsize - loop2, &hairpin_seq[0], &hairpin_seq[loop2]);
			newLoop[1] = new HairpinLoop(pt, loop2 - loop - 1, &hairpin_seq[loop]);
			newLoop[0]->addAdjacent(adjacentLoops[0]);
			adjacentLoops[0]->replaceAdjacent(this, newLoop[0]);
			newLoop[0]->addAdjacent(newLoop[1]);
			newLoop[1]->addAdjacent(newLoop[0]);
			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			adjacentLoops[0]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}

	}

	return -totalRate;
}

void HairpinLoop::generateMoves(void) {
	double energies[2];
	int pt = 0;
	int loop, loop2;
	double tempRate = 0;
	RateEnv rateEnv;

// Creation moves
	if (hairpinsize <= 4) {
		// We cannot form any creation moves in the hairpin unless it has at least 5 bases.
		if (moves != NULL)
			delete moves;

		moves = new MoveList(0);
		totalRate = 0.0;
		generateDeleteMoves();
		return;
	} else {
		if (moves != NULL)
			delete moves;
		moves = new MoveList(1);

		// Indice 0 is the starting hairpin base. hairpinsize+1 is the ending hairpin base. Thus we want to start at hairpin indice 1, and go to hairpinsize - 3. (which could pair to indice hairpinsize)
		for (loop = 1; loop <= hairpinsize - 4; loop++)
			for (loop2 = loop + 4; loop2 <= hairpinsize; loop2++) {

				pt = pairtypes[hairpin_seq[loop]][hairpin_seq[loop2]];

				if (pt != 0) { // the two could pair. Work out energies of the resulting pair of loops.

					// Case handling time.
					// possibilities: new stack + hairpin.
					//                bulge + hairpin
					//                interior + hairpin

					// new stack + hairpin
					if (loop == 1 && loop2 == hairpinsize) {

						energies[0] = energyModel->StackEnergy(hairpin_seq[0], hairpin_seq[hairpinsize + 1], hairpin_seq[loop], hairpin_seq[loop2]);
						energies[1] = energyModel->HairpinEnergy(&hairpin_seq[1], hairpinsize - 2);
						tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

						// stack and hairpin, so this is loop and stack
						rateEnv = RateEnv(tempRate, energyModel, loopMove, stackMove);
						moves->addMove(new Move(MOVE_CREATE | MOVE_1, rateEnv, this, loop, loop2));

					}
					// bulge + hairpin
					else if (loop == 1 || loop2 == hairpinsize) {
						// total bulge size is the difference between either loop and 1 (if that's the edge that has the bulge), or hairpinsize and loop2 (if that's the edge that moved). One must have moved with the other stationary to have a bulge, so the sum will always total whichever moved.
						energies[0] = energyModel->BulgeEnergy(hairpin_seq[0], hairpin_seq[hairpinsize + 1], hairpin_seq[loop], hairpin_seq[loop2],
								(loop - 1) + (hairpinsize - loop2));

						// loop2 - loop - 1 is the new hairpin size.

						energies[1] = energyModel->HairpinEnergy(&hairpin_seq[loop], loop2 - loop - 1);

						tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

						// new bulgeloop + hairpin: this is openMove and stackLoopMove

						rateEnv = RateEnv(tempRate, energyModel, loopMove, stackLoopMove);
						moves->addMove(new Move(MOVE_CREATE | MOVE_2, rateEnv, this, loop, loop2));

					} else // interior loop + hairpin case.
					{
						//		  char mismatches[4] = { hairpin_seq[1], hairpin_seq[hairpinsize], hairpin_seq[loop-1], hairpin_seq[loop2+1]};
						energies[0] = energyModel->InteriorEnergy(hairpin_seq, &hairpin_seq[loop2], loop - 1, hairpinsize - loop2);

						//pairtype, pt, loop - 1, hairpinsize - loop2, mismatches );

						// loop2 - loop - 1 is the new hairpin size.
						energies[1] = energyModel->HairpinEnergy(&hairpin_seq[loop], loop2 - loop - 1);
						tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

						// interiorLoop + hairpin, so this is open + open

						rateEnv = RateEnv(tempRate, energyModel, loopMove, loopMove);

						moves->addMove(new Move(MOVE_CREATE | MOVE_3, rateEnv, this, loop, loop2));
					}
				}
			}
		totalRate = moves->getRate();
	}

// Shift moves

// Delete moves
	generateDeleteMoves();
}

void HairpinLoop::generateDeleteMoves(void) {
	double temprate;

	generateAndSaveDeleteMove(adjacentLoops[0], 0);

	totalRate = moves->getRate();
}

void HairpinLoop::printMove(Loop *comefrom, char *structure_p, char *seq_p) {
	int loop;
	structure_p[hairpin_seq - seq_p] = '(';
	structure_p[hairpin_seq - seq_p + 1 + hairpinsize] = ')';

	for (loop = 0; loop < curAdjacent; loop++) {
		if (adjacentLoops[loop] != comefrom && adjacentLoops[loop] != NULL)
			// shouldn't happen, being careful.
			adjacentLoops[loop]->printMove(this, structure_p, seq_p);
		assert(adjacentLoops[loop] != NULL);
	}
}

char *HairpinLoop::getLocation(Move *move, int index) {
	if (move->getType() & MOVE_CREATE) {
		if (index == 0)
			return &hairpin_seq[move->index[0]];
		else
			return &hairpin_seq[move->index[1]];

	} else if (move->getType() & MOVE_DELETE)
		return hairpin_seq;

	assert(0);
}

char *HairpinLoop::verifyLoop(char *incoming_sequence, int incoming_pairtype, Loop *from) {
	char *ret_seq;
	if (adjacentLoops[0] == from) {
		if (incoming_sequence != hairpin_seq) {
			fprintf(stderr, "Verification Failed\n");
			assert(incoming_sequence == hairpin_seq);
		}
		return hairpin_seq + hairpinsize + 1;
	} else
		assert(0);
	return NULL;
}

/*
 BulgeLoop Functions
 */

BulgeLoop::BulgeLoop(void) {
	numAdjacent = 2;
	adjacentLoops = new Loop *[2];

	bulgesize[0] = 0;
	bulgesize[1] = 0;
	bulge_seq[0] = NULL;
	bulge_seq[1] = NULL;
	identity = 'B';
}

BulgeLoop::BulgeLoop(int type1, int type2, int size1, int size2, char *bulge_sequence1, char *bulge_sequence2, Loop *left, Loop *right) {
	numAdjacent = 2;
	curAdjacent = 0;
	adjacentLoops = new Loop *[2];
	adjacentLoops[0] = left;
	adjacentLoops[1] = right;
	if (left != NULL)
		curAdjacent++;
	if (right != NULL)
		curAdjacent++;

	bulgesize[0] = size1;
	bulgesize[1] = size2;
	bulge_seq[0] = bulge_sequence1;
	bulge_seq[1] = bulge_sequence2;
	pairtype[0] = type1;
	pairtype[1] = type2;
	identity = 'B';
}

string BulgeLoop::typeInternalsToString(void) {

	std::stringstream ss;

	ss << " seq0=" << utility::sequenceToString(bulge_seq[0], bulgesize[0]) << " \n";
	ss << " seq1= " << utility::sequenceToString(bulge_seq[1], bulgesize[1]) << "\n";

	return ss.str();

}

Move *BulgeLoop::getChoice(double *randomchoice, Loop *from) {
	Move *stor;
	assert(randomchoice != NULL);
	assert(*randomchoice >= 0.0); // never should see a negative choice value.

	if (*randomchoice < totalRate) // something was chosen, do this
		return moves->getChoice(randomchoice);
	else {
		*randomchoice -= totalRate;
		if (adjacentLoops[0] != from) {
			stor = adjacentLoops[0]->getChoice(randomchoice, this);
			if (stor != NULL)
				return stor;
		}
		if (adjacentLoops[1] != from) {
			stor = adjacentLoops[1]->getChoice(randomchoice, this);
			if (stor != NULL)
				return stor;
		}
	}
	return NULL;
}

void BulgeLoop::calculateEnergy(void) {
	if (energyModel == NULL)
		return; // we can't handle this error. I'm trying to work out a way around it, but generally if the loops try to get used before the energy model initializes, it's all over.

	energy = energyModel->BulgeEnergy(bulge_seq[0][0], bulge_seq[1][bulgesize[1] + 1], bulge_seq[0][bulgesize[0] + 1], bulge_seq[1][0],
			bulgesize[0] + bulgesize[1]);
	return;
}

double BulgeLoop::doChoice(Move *move, Loop **returnLoop) {
	Loop *newLoop[2];
	int pt, loop, loop2;
	int *ptypes = new int[3];
	int *sidelen = new int[3];
	char **seqs = new char *[3];
	int bsize = bulgesize[0] + bulgesize[1];
	int bside = (bulgesize[0] == 0) ? 1 : 0;

	if (move->type & MOVE_CREATE) {
		loop = move->index[0];
		loop2 = move->index[1];
		pt = pairtypes[bulge_seq[bside][loop]][bulge_seq[bside][loop2]];
		if (bside == 1) {
			ptypes[0] = pairtype[0];
			ptypes[1] = pairtype[1];
			ptypes[2] = pt;
			sidelen[0] = 0;
			sidelen[1] = loop - 1;
			sidelen[2] = bsize - loop2;
			seqs[0] = bulge_seq[0];
			seqs[1] = bulge_seq[1];
			seqs[2] = &bulge_seq[1][loop2];
		} else {
			ptypes[0] = pairtype[0];
			ptypes[1] = pt;
			ptypes[2] = pairtype[1];
			sidelen[0] = loop - 1;
			sidelen[1] = bsize - loop2;
			sidelen[2] = 0;
			seqs[0] = bulge_seq[0];
			seqs[1] = &bulge_seq[0][loop2];
			seqs[2] = bulge_seq[1];
		}
		// creation moves in a bulge loop are always multi+hairpin
		newLoop[0] = new MultiLoop(3, ptypes, sidelen, seqs);
		newLoop[1] = new HairpinLoop(pt, loop2 - loop - 1, &bulge_seq[bside][loop]);
		newLoop[0]->addAdjacent(adjacentLoops[0]);
		adjacentLoops[0]->replaceAdjacent(this, newLoop[0]);
		if (bside == 0) {
			newLoop[0]->addAdjacent(newLoop[1]);
			newLoop[0]->addAdjacent(adjacentLoops[1]);
		} else {
			newLoop[0]->addAdjacent(adjacentLoops[1]);
			newLoop[0]->addAdjacent(newLoop[1]);
		}
		adjacentLoops[1]->replaceAdjacent(this, newLoop[0]);
		newLoop[1]->addAdjacent(newLoop[0]);
		newLoop[0]->generateMoves();
		newLoop[1]->generateMoves();
		adjacentLoops[0]->generateMoves();
		adjacentLoops[1]->generateMoves();
		*returnLoop = newLoop[0];
		return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
	}
	return -totalRate;
}

void BulgeLoop::generateMoves(void) {
	double energies[2];
	int loop, loop2, pt;
	double tempRate;
	RateEnv rateEnv;
	int bsize = bulgesize[0] + bulgesize[1];
	int bside = (bulgesize[0] == 0) ? 1 : 0;

//	std::cout << "Trying to generate moves";
//	std: cout.flush();

// Creation moves
	if (bsize <= 3) {
		if (moves != NULL)
			delete moves;

		moves = new MoveList(0);
		totalRate = 0.0;
		generateDeleteMoves();
		return;
	} else {
		if (moves != NULL)
			delete moves;
		moves = new MoveList(bsize); // what's the optimal #?

		// Indice 0 is the starting bulge base. bulgesize+1 is the ending hairpin base. Thus we want to start at hairpin indice 1, and go to hairpinsize - 4. (which could pair to indice hairpinsize)
		for (loop = 1; loop <= bsize - 4; loop++)
			for (loop2 = loop + 4; loop2 <= bsize; loop2++) {

				pt = pairtypes[bulge_seq[bside][loop]][bulge_seq[bside][loop2]];

				if (pt != 0) { // the two could pair. Work out energies of the resulting pair of loops.

					// Case handling time.
					// it will always be a multiloop and hairpin.

					// Multiloop energy - CHECK THIS/FIXME
					int ptypes[3];
					int sidelen[3];
					char *sequences[3];

					if (bside == 0) {
						ptypes[0] = pairtype[0];
						sidelen[0] = loop - 1;
						sequences[0] = bulge_seq[0];

						ptypes[1] = pt;
						sidelen[1] = bsize - loop2;
						sequences[1] = &bulge_seq[0][loop2];

						ptypes[2] = pairtype[1];
						sidelen[2] = 0;
						sequences[2] = bulge_seq[1];

					} else {
						ptypes[0] = pairtype[0];
						sidelen[0] = 0;
						sequences[0] = bulge_seq[0];

						ptypes[1] = pairtype[1];
						sidelen[1] = loop - 1;
						sequences[1] = bulge_seq[1];

						ptypes[2] = pt;
						sidelen[2] = bsize - loop2;
						sequences[2] = &bulge_seq[1][loop2];

					}
					// need to add sequence info -definate FIXME for dangles != 0
					energies[0] = energyModel->MultiloopEnergy(3, sidelen, sequences);
					// loop2 - loop + 1 is the new hairpin size.
					energies[1] = energyModel->HairpinEnergy(&bulge_seq[bside][loop], loop2 - loop - 1);

					tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

					// hairpin and multiloop, so this is loopMove and something

					MoveType multiMove = stackMove; // default init value;

					if (bside == 0) {
						multiMove = energyModel->prefactorInternal(sidelen[0], sidelen[1]);
					} else {
						multiMove = energyModel->prefactorInternal(sidelen[1], sidelen[2]);
					}

					rateEnv = RateEnv(tempRate, energyModel, loopMove, multiMove);

					moves->addMove(new Move(MOVE_CREATE, rateEnv, this, loop, loop2));
				}
			}
	}
	totalRate = moves->getRate();

	generateDeleteMoves();
}

void BulgeLoop::generateDeleteMoves(void) {

	double temprate;

//	std::cout << "Trying to generate delete moves";
//	std: cout.flush();

	assert(moves != NULL);

	generateAndSaveDeleteMove(adjacentLoops[0], 0);
	generateAndSaveDeleteMove(adjacentLoops[1], 1);

	totalRate = moves->getRate();
}

void BulgeLoop::printMove(Loop *comefrom, char *structure_p, char *seq_p) {
	int loop;
	int item = (bulge_seq[0] < bulge_seq[1]);

	structure_p[bulge_seq[0] - seq_p] = item ? '(' : ')';
	structure_p[bulge_seq[1] - seq_p + 1 + bulgesize[1]] = item ? ')' : '(';
	structure_p[bulge_seq[0] - seq_p + 1 + bulgesize[0]] = item ? '(' : ')';
	structure_p[bulge_seq[1] - seq_p] = item ? ')' : '(';
	/**/
	for (loop = 0; loop < curAdjacent; loop++) {
		if (adjacentLoops[loop] != comefrom && adjacentLoops[loop] != NULL)
			// shouldn't happen, being careful.
			adjacentLoops[loop]->printMove(this, structure_p, seq_p);
		assert(adjacentLoops[loop] != NULL);
	}
}

char *BulgeLoop::getLocation(Move *move, int index) {
	int bside = (bulgesize[0] == 0) ? 1 : 0;
	if (move->getType() & MOVE_CREATE) {
		return &bulge_seq[bside][move->index[index]];
	} else if (move->getType() & MOVE_DELETE) {
		if (adjacentLoops[0] == move->affected[0] || adjacentLoops[0] == move->affected[1])
			return bulge_seq[0];
		else if (adjacentLoops[1] == move->affected[0] || adjacentLoops[1] == move->affected[1])
			return bulge_seq[1];
	}
	assert(0);
}

char *BulgeLoop::verifyLoop(char *incoming_sequence, int incoming_pairtype, Loop *from) {
	char *ret_seq;
	if (adjacentLoops[0] == from) {
		if (incoming_sequence != bulge_seq[0]) {
			fprintf(stderr, "Verification Failed\n");
			assert(incoming_sequence == bulge_seq[0]);
		}
		ret_seq = adjacentLoops[1]->verifyLoop(bulge_seq[0] + 1 + bulgesize[0], pairtype[1], this);
		assert(ret_seq == bulge_seq[1]);
		return bulge_seq[1] + 1 + bulgesize[1];
	} else if (adjacentLoops[1] == from) {
		if (incoming_sequence != bulge_seq[1]) {
			fprintf(stderr, "Verification Failed\n");
			assert(incoming_sequence == bulge_seq[1]);
		}
		ret_seq = adjacentLoops[0]->verifyLoop(bulge_seq[1] + 1 + bulgesize[1], pairtype[0], this);
		assert(ret_seq == bulge_seq[0]);
		return bulge_seq[0] + 1 + bulgesize[0];
	} else
		assert(0);
	return NULL;
}

/*
 InteriorLoop functions
 */

InteriorLoop::InteriorLoop(void) {
	numAdjacent = 2;
	adjacentLoops = new Loop *[2];

	sizes[0] = sizes[1] = 0;
	int_seq[0] = int_seq[1] = NULL;
	identity = 'I';
}

InteriorLoop::InteriorLoop(int type1, int type2, int size1, int size2, char *int_seq1, char *int_seq2, Loop *left, Loop *right) {
	numAdjacent = 2;
	adjacentLoops = new Loop *[2];

	adjacentLoops[0] = left;
	adjacentLoops[1] = right;
	if (left != NULL)
		curAdjacent++;
	if (right != NULL)
		curAdjacent++;

	pairtype[0] = type1;
	pairtype[1] = type2;
	sizes[0] = size1;
	sizes[1] = size2;
	int_seq[0] = int_seq1;
	int_seq[1] = int_seq2;
	identity = 'I';

}

string InteriorLoop::typeInternalsToString(void) {

	std::stringstream ss;

	ss << " seq0=" << utility::sequenceToString(int_seq[0], sizes[0] + 1) << " \n";
	ss << " seq1= " << utility::sequenceToString(int_seq[1], sizes[1] + 1) << "\n";

	return ss.str();

}

Move *InteriorLoop::getChoice(double *randomchoice, Loop *from) {
	Move *stor;
	assert(randomchoice != NULL);
	assert(*randomchoice >= 0.0); // never should see a negative choice value.

	if (*randomchoice < totalRate) // something was chosen, do this
		return moves->getChoice(randomchoice);
	else {
		*randomchoice -= totalRate;
		if (adjacentLoops[0] != from) {
			stor = adjacentLoops[0]->getChoice(randomchoice, this);
			if (stor != NULL)
				return stor;
		}
		if (adjacentLoops[1] != from) {
			stor = adjacentLoops[1]->getChoice(randomchoice, this);
			if (stor != NULL)
				return stor;
		}
	}
	return NULL;
}

void InteriorLoop::calculateEnergy(void) {
	char mismatches[4];
	if (int_seq[0] == NULL || int_seq[1] == NULL)
		return;

// CHECK THESE
	mismatches[0] = int_seq[0][1];
	mismatches[1] = int_seq[1][sizes[1]];
	mismatches[2] = int_seq[0][sizes[0]];
	mismatches[3] = int_seq[1][1];

	if (energyModel == NULL)
		return; // we can't handle this error. I'm trying to work out a way around it, but generally if the loops try to get used before the energy model initializes, it's all over.

	energy = energyModel->InteriorEnergy(int_seq[0], int_seq[1], sizes[0], sizes[1]);
	return;
}

double InteriorLoop::doChoice(Move *move, Loop **returnLoop) {
	Loop *newLoop[2];
	int pt, loop, loop2;
	int *ptypes = new int[3];
	int *sidelen = new int[3];
	char **seqs = new char *[3];

	if (move->type & MOVE_CREATE) {
		if (move->type & MOVE_1) {
			loop = move->index[0];
			loop2 = move->index[1];
			pt = pairtypes[int_seq[0][loop]][int_seq[0][loop2]];

			ptypes[0] = pairtype[0];
			ptypes[1] = pt;
			ptypes[2] = pairtype[1];
			sidelen[0] = loop - 1;
			sidelen[1] = sizes[0] - loop2;
			sidelen[2] = sizes[1];
			seqs[0] = int_seq[0];
			seqs[1] = &int_seq[0][loop2];
			seqs[2] = int_seq[1];

			// creation moves in a interior loop's sides are always multi+hairpin
			newLoop[0] = new MultiLoop(3, ptypes, sidelen, seqs);
			newLoop[1] = new HairpinLoop(pt, loop2 - loop - 1, &int_seq[0][loop]);

			newLoop[0]->addAdjacent(adjacentLoops[0]);
			adjacentLoops[0]->replaceAdjacent(this, newLoop[0]);
			newLoop[0]->addAdjacent(newLoop[1]);
			newLoop[0]->addAdjacent(adjacentLoops[1]);
			adjacentLoops[1]->replaceAdjacent(this, newLoop[0]);

			newLoop[1]->addAdjacent(newLoop[0]);

			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			adjacentLoops[0]->generateMoves();
			adjacentLoops[1]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}

		else if (move->type & MOVE_2) {
			loop = move->index[0];
			loop2 = move->index[1];
			pt = pairtypes[int_seq[1][loop]][int_seq[1][loop2]];

			ptypes[0] = pairtype[0];
			ptypes[1] = pairtype[1];
			ptypes[2] = pt;
			sidelen[0] = sizes[0];
			sidelen[1] = loop - 1;
			sidelen[2] = sizes[1] - loop2;
			seqs[0] = int_seq[0];
			seqs[1] = int_seq[1];
			seqs[2] = &int_seq[1][loop2];

			// creation moves in a bulge loop are always multi+hairpin
			newLoop[0] = new MultiLoop(3, ptypes, sidelen, seqs);
			newLoop[1] = new HairpinLoop(pt, loop2 - loop - 1, &int_seq[1][loop]);

			newLoop[0]->addAdjacent(adjacentLoops[0]);
			adjacentLoops[0]->replaceAdjacent(this, newLoop[0]);
			newLoop[0]->addAdjacent(adjacentLoops[1]);
			adjacentLoops[1]->replaceAdjacent(this, newLoop[0]);
			newLoop[0]->addAdjacent(newLoop[1]);

			newLoop[1]->addAdjacent(newLoop[0]);

			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			adjacentLoops[0]->generateMoves();
			adjacentLoops[1]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		} else if (move->type & MOVE_3) {
			loop = move->index[0];
			loop2 = move->index[1];
			pt = pairtypes[int_seq[0][loop]][int_seq[1][loop2]];
			// Need to check conditions for each side in order to determine what the two new loops types would be.
			// adjacent to first pair side:
			if (loop == 1 && loop2 == sizes[1]) // stack
				newLoop[0] = new StackLoop(pairtype[0], pt, int_seq[0], &int_seq[1][loop2]);
			else if (loop == 1 || loop2 == sizes[1]) // bulge
				newLoop[0] = new BulgeLoop(pairtype[0], pt, (loop - 1), (sizes[1] - loop2), int_seq[0], &int_seq[1][loop2]);
			else
				// interior
				newLoop[0] = new InteriorLoop(pairtype[0], pt, loop - 1, sizes[1] - loop2, int_seq[0], &int_seq[1][loop2]);

			// other side
			if (loop == sizes[0] && loop2 == 1) {		// stack

				newLoop[1] = new StackLoop(pt, pairtype[1], &int_seq[0][loop], int_seq[1]);

			} else {

				if (loop == sizes[0] || loop2 == 1) 	// bulge

					newLoop[1] = new BulgeLoop(pt, pairtype[1], (sizes[0] - loop), (loop2 - 1), &int_seq[0][loop], int_seq[1]);
				else
					// interior
					newLoop[1] = new InteriorLoop(pt, pairtype[1], sizes[0] - loop, loop2 - 1, &int_seq[0][loop], int_seq[1]);

			}
			newLoop[0]->addAdjacent(adjacentLoops[0]);
			adjacentLoops[0]->replaceAdjacent(this, newLoop[0]);
			newLoop[0]->addAdjacent(newLoop[1]);
			newLoop[1]->addAdjacent(newLoop[0]);
			newLoop[1]->addAdjacent(adjacentLoops[1]);
			adjacentLoops[1]->replaceAdjacent(this, newLoop[1]);

			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();

			adjacentLoops[0]->generateMoves();
			adjacentLoops[1]->generateMoves();

			*returnLoop = newLoop[0];

			delete[] seqs;
			delete[] sidelen;
			delete[] ptypes;
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		} else {
			delete[] seqs;
			delete[] sidelen;
			delete[] ptypes;
		}
	} else {
		delete[] seqs;
		delete[] sidelen;
		delete[] ptypes;
	}
	return -totalRate;
}

void InteriorLoop::generateMoves(void) {
	double energies[2];
	int pt = 0;
	int loop, loop2;
	double tempRate = 0;
	RateEnv rateEnv;

	int nummoves = sizes[0] * sizes[1];
	if (sizes[0] > 4)
		nummoves += sizes[0] - 4;
	if (sizes[1] > 4)
		nummoves += sizes[1] - 4;

// this number is wrong... sigh...
	nummoves = (int) ((sizes[0] * sizes[1]) / 16 + 1);

// Creation moves
	if (moves != NULL)
		delete moves;
	moves = new MoveList(nummoves);

// three loops here, the first is only side 0's possible creation moves
//                   the second is only side 1's possible creation moves
//                   the third is only creation moves that cross 0-1.

// Loop #1: Side 0 only Creation Moves
	for (loop = 1; loop <= sizes[0] - 4; loop++) {

		for (loop2 = loop + 4; loop2 <= sizes[0]; loop2++) { // each possibility will always result in a new hairpin + multiloop.

			pt = pairtypes[int_seq[0][loop]][int_seq[0][loop2]];

			if (pt != 0) {
				energies[0] = energyModel->HairpinEnergy(&int_seq[0][loop], loop2 - loop - 1);

				// Multiloop energy

				int ptypes[3] = { pairtype[0], pt, pairtype[1] };
				int sidelen[3] = { loop - 1, sizes[0] - loop2, sizes[1] };
				char *sequences[3] = { &int_seq[0][0], &int_seq[0][loop2], &int_seq[1][0] };

				energies[1] = energyModel->MultiloopEnergy(3, sidelen, sequences);
				tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

//				// hairpin and multiloop, so this is loopMove and something

				MoveType multiMove = energyModel->prefactorInternal(sidelen[0], sidelen[1]);
				rateEnv = RateEnv(tempRate, energyModel, multiMove, loopMove);

				moves->addMove(new Move(MOVE_CREATE | MOVE_1, rateEnv, this, loop, loop2));
			}
		}
	}

// Loop #2: Side 1 only Creation Moves
	for (loop = 1; loop <= sizes[1] - 4; loop++)
		for (loop2 = loop + 4; loop2 <= sizes[1]; loop2++) { // each possibility will always result in a new hairpin + multiloop.
			pt = pairtypes[int_seq[1][loop]][int_seq[1][loop2]];
			if (pt != 0) {
				energies[0] = energyModel->HairpinEnergy(&int_seq[1][loop], loop2 - loop - 1);

				// Multiloop energy - CHECK THIS

				int ptypes[3] = { pairtype[0], pairtype[1], pt };
				int sidelen[3] = { sizes[0], loop - 1, sizes[1] - loop2 };
				char *sequences[3] = { &int_seq[0][0], &int_seq[1][0], &int_seq[1][loop2] };

				energies[1] = energyModel->MultiloopEnergy(3, sidelen, sequences);
				tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

				// hairpin and multiloop, so this is loopMove and something

				MoveType multiMove = energyModel->prefactorInternal(sidelen[1], sidelen[2]);
				rateEnv = RateEnv(tempRate, energyModel, loopMove, multiMove);

				moves->addMove(new Move(MOVE_CREATE | MOVE_2, rateEnv, this, loop, loop2));
			}
		}

// Loop #3: Side 0 to Side 1 crossing moves ONLY

	for (loop = 1; loop <= sizes[0]; loop++)
		for (loop2 = 1; loop2 <= sizes[1]; loop2++) {

			pt = pairtypes[int_seq[0][loop]][int_seq[1][loop2]];
			if (pt != 0) {
				// Need to check conditions for each side in order to determine what the two new loops types would be.
				// adjacent to first pair side:

				MoveType leftMove = stackMove;
				MoveType rightMove = stackMove;

				if (loop == 1 && loop2 == sizes[1]) {			// stack
					energies[0] = energyModel->StackEnergy(int_seq[0][0], int_seq[1][sizes[1] + 1], int_seq[0][loop], int_seq[1][loop2]);
				} else if (loop == 1 || loop2 == sizes[1]) {		// bulge
					energies[0] = energyModel->BulgeEnergy(int_seq[0][0], int_seq[1][sizes[1] + 1], int_seq[0][loop], int_seq[1][loop2],
							(loop - 1) + (sizes[1] - loop2));
					leftMove = stackLoopMove;
				} else {  // interior
					energies[0] = energyModel->InteriorEnergy(int_seq[0], &int_seq[1][loop2], loop - 1, sizes[1] - loop2);
					leftMove = loopMove;
				}

				// other side
				if (loop == sizes[0] && loop2 == 1) { // stack
					energies[1] = energyModel->StackEnergy(int_seq[0][loop], int_seq[1][loop2], int_seq[0][sizes[0] + 1], int_seq[1][0]);
				} else if (loop == sizes[0] || loop2 == 1) { // bulge
					energies[1] = energyModel->BulgeEnergy(int_seq[0][loop], int_seq[1][loop2], int_seq[0][sizes[0] + 1], int_seq[1][0],
							(loop2 - 1) + (sizes[0] - loop));
					rightMove = stackLoopMove;
				} else { // interior
					energies[1] = energyModel->InteriorEnergy(&int_seq[0][loop], int_seq[1], sizes[0] - loop, loop2 - 1);
					rightMove = loopMove;
				}

				tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

				// interior loop is closing, so this could be anything.
				rateEnv = RateEnv(tempRate, energyModel, leftMove, rightMove);

				moves->addMove(new Move(MOVE_CREATE | MOVE_3, rateEnv, this, loop, loop2));
			}
		}

// totaling the rate
	totalRate = moves->getRate();

// Shift moves

// Delete moves
	generateDeleteMoves();
}

void InteriorLoop::generateDeleteMoves(void) {
	double temprate;

	assert(moves != NULL);

	generateAndSaveDeleteMove(adjacentLoops[0], 0);
	generateAndSaveDeleteMove(adjacentLoops[1], 1);

	totalRate = moves->getRate();

}

char *InteriorLoop::getLocation(Move *move, int index) {
	if (move->getType() & MOVE_CREATE) {
		if (move->getType() & MOVE_1)
			return &int_seq[0][move->index[index]];
		if (move->getType() & MOVE_2)
			return &int_seq[1][move->index[index]];
		if (move->getType() & MOVE_3)
			return &int_seq[index][move->index[index]];
	} else if (move->getType() & MOVE_DELETE) {
		if (adjacentLoops[0] == move->affected[0] || adjacentLoops[0] == move->affected[1])
			return int_seq[0];
		else if (adjacentLoops[1] == move->affected[0] || adjacentLoops[1] == move->affected[1])
			return int_seq[1];
	}
	assert(0);
}

char *InteriorLoop::verifyLoop(char *incoming_sequence, int incoming_pairtype, Loop *from) {
	char *ret_seq;
	if (adjacentLoops[0] == from) {
		if (incoming_sequence != int_seq[0]) {
			fprintf(stderr, "Verification Failed\n");
			assert(incoming_sequence == int_seq[0]);
		}
		ret_seq = adjacentLoops[1]->verifyLoop(int_seq[0] + 1 + sizes[0], pairtype[1], this);
		assert(ret_seq == int_seq[1]);
		return int_seq[1] + 1 + sizes[1];
	} else if (adjacentLoops[1] == from) {
		if (incoming_sequence != int_seq[1]) {
			fprintf(stderr, "Verification Failed\n");
			assert(incoming_sequence == int_seq[1]);
		}
		ret_seq = adjacentLoops[0]->verifyLoop(int_seq[1] + 1 + sizes[1], pairtype[0], this);
		assert(ret_seq == int_seq[0]);
		return int_seq[0] + 1 + sizes[0];
	} else
		assert(0);
	return NULL;
}

void InteriorLoop::printMove(Loop *comefrom, char *structure_p, char *seq_p) {
	int loop;
	int item = (int_seq[0] < int_seq[1]);
	structure_p[int_seq[0] - seq_p] = item ? '(' : ')';
	structure_p[int_seq[1] - seq_p + 1 + sizes[1]] = item ? ')' : '(';
	structure_p[int_seq[0] - seq_p + 1 + sizes[0]] = item ? '(' : ')';
	structure_p[int_seq[1] - seq_p] = item ? ')' : '(';

	for (loop = 0; loop < curAdjacent; loop++) {
		if (adjacentLoops[loop] != comefrom && adjacentLoops[loop] != NULL)
			// shouldn't happen, being careful.
			adjacentLoops[loop]->printMove(this, structure_p, seq_p);
		assert(adjacentLoops[loop] != NULL);
	}
}

/* MultiLoop functions */

MultiLoop::MultiLoop(void) {
	numAdjacent = 0;
	adjacentLoops = NULL;

	pairtype = NULL;
	sidelen = NULL;
	seqs = NULL;
	identity = 'M';
}

MultiLoop::MultiLoop(int branches, int *pairtypes, int *sidelengths, char **sequences) {
	numAdjacent = branches;
	adjacentLoops = new Loop *[branches];
	for (int loop = 0; loop < branches; loop++)
		adjacentLoops[loop] = NULL;
	pairtype = pairtypes;
	sidelen = sidelengths;
	seqs = sequences;
	identity = 'M';
}

MultiLoop::~MultiLoop(void) {
	delete[] pairtype;
	delete[] sidelen;
	delete[] seqs;
}

string MultiLoop::typeInternalsToString(void) {

	std::stringstream ss;

	ss << "sideLength = ";

	for (int i = 0; i < numAdjacent + 1; i++) {

		ss << sidelen[i];

	}

	ss << "\n";

	return ss.str();

}

Move *MultiLoop::getChoice(double* randomchoice, Loop *from) {
	Move *stor;
	assert(randomchoice != NULL);
	assert(*randomchoice >= 0.0); // never should see a negative choice value.

	if (*randomchoice < totalRate) // something was chosen, do this
		return moves->getChoice(randomchoice);
	else {
		*randomchoice -= totalRate;
		for (int loop = 0; loop < curAdjacent; loop++)
			if (adjacentLoops[loop] != from) {
				stor = adjacentLoops[loop]->getChoice(randomchoice, this);
				if (stor != NULL)
					return stor;
			}
	}
	return NULL;
}

void MultiLoop::calculateEnergy(void) {
	if (energyModel == NULL)
		return; // we can't handle this error. I'm trying to work out a way around it, but generally if the loops try to get used before the energy model initializes, it's all over.

	energy = energyModel->MultiloopEnergy(numAdjacent, sidelen, seqs);
	return;
}

double MultiLoop::doChoice(Move *move, Loop **returnLoop) {
	Loop *newLoop[2];
	int pt, loop, loop2, loop3, loop4, temploop, tempindex;
	int *ptypes;
	int *sidelengths;
	char **sequences;

	if (move->type & MOVE_CREATE) {
		loop = move->index[0];
		loop2 = move->index[1];
		loop3 = move->index[2];
		loop4 = move->index[3];

		if (move->type & MOVE_1) {
			//single side, hairpin + multi with 1 higher mag.
			ptypes = new int[numAdjacent + 1];
			sidelengths = new int[numAdjacent + 1];
			sequences = new char *[numAdjacent + 1];

			pt = pairtypes[seqs[loop3][loop]][seqs[loop3][loop2]];

			for (temploop = 0, tempindex = 0; temploop < numAdjacent + 1; temploop++, tempindex++) {
				if (temploop == loop3) {
					ptypes[temploop] = pairtype[temploop];
					sidelengths[temploop] = loop - 1;
					sequences[temploop] = seqs[temploop];
					ptypes[temploop + 1] = pt;
					sidelengths[temploop + 1] = sidelen[temploop] - loop2;
					sequences[temploop + 1] = &seqs[temploop][loop2];
					temploop = temploop + 1;
				} else {
					ptypes[temploop] = pairtype[tempindex];
					sidelengths[temploop] = sidelen[tempindex];
					sequences[temploop] = seqs[tempindex];
				}
			}

			newLoop[0] = new MultiLoop(numAdjacent + 1, ptypes, sidelengths, sequences);

			newLoop[1] = new HairpinLoop(pt, loop2 - loop - 1, &seqs[loop3][loop]);

			for (temploop = 0; temploop < numAdjacent; temploop++) {
				if (temploop == loop3) {
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
					newLoop[0]->addAdjacent(newLoop[1]);
				} else {
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
				}
			}

			newLoop[1]->addAdjacent(newLoop[0]);
			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}

		// MOVE_2
		if (move->type & MOVE_2) {
			//adjacent sides, one of: stack    + multi with same mag
			//                        bulge    + multi with same mag
			//                        interior + multi with same mag
			ptypes = new int[numAdjacent];
			sidelengths = new int[numAdjacent];
			sequences = new char *[numAdjacent];

			loop4 = (loop3 + 1) % numAdjacent;
			pt = pairtypes[seqs[loop3][loop]][seqs[loop4][loop2]];
			for (temploop = 0; temploop < numAdjacent; temploop++) {
				if (temploop == loop3) {
					ptypes[temploop] = pairtype[temploop];
					sidelengths[temploop] = loop - 1;
					sequences[temploop] = seqs[temploop];
				} else {
					ptypes[temploop] = pairtype[temploop];
					if (temploop == loop4) {
						ptypes[temploop] = pt;
						sidelengths[temploop] = sidelen[temploop] - loop2;
						sequences[temploop] = &seqs[temploop][loop2];
					} else {
						sidelengths[temploop] = sidelen[temploop];
						sequences[temploop] = seqs[temploop];
					}
				}
			}
			newLoop[0] = new MultiLoop(numAdjacent, ptypes, sidelengths, sequences);

			// three cases for which type of move:
			// #2a: stack
			if (loop == sidelen[loop3] && loop2 == 1) {
				if (loop3 < loop4)
					newLoop[1] = new StackLoop(pt, pairtype[loop4], &seqs[loop3][loop], &seqs[loop4][loop2 - 1]);
				else
					newLoop[1] = new StackLoop(pairtype[loop4], pt, &seqs[loop4][loop2 - 1], &seqs[loop3][loop]);
			}
			// #2b: bulge
			else if (loop == sidelen[loop3] || loop2 == 1) {
				if (loop3 < loop4)
					newLoop[1] = new BulgeLoop(pt, pairtype[loop4], sidelen[loop3] - loop, loop2 - 1, &seqs[loop3][loop], &seqs[loop4][0]);
				else
					newLoop[1] = new BulgeLoop(pairtype[loop4], pt, loop2 - 1, sidelen[loop3] - loop, &seqs[loop4][0], &seqs[loop3][loop]);
			}

			// #2c: interior
			else {
				if (loop3 < loop4)
					newLoop[1] = new InteriorLoop(pt, pairtype[loop4], sidelen[loop3] - loop, loop2 - 1, &seqs[loop3][loop], &seqs[loop4][0]);
				else
					newLoop[1] = new InteriorLoop(pairtype[loop4], pt, loop2 - 1, sidelen[loop3] - loop, &seqs[loop4][0], &seqs[loop3][loop]);
			}

			for (temploop = 0; temploop < numAdjacent; temploop++) {
				if (temploop == loop4) {
					newLoop[0]->addAdjacent(newLoop[1]);
				} else {
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
				}
			}

			if (loop3 < loop4) {
				newLoop[1]->addAdjacent(newLoop[0]);
				newLoop[1]->addAdjacent(adjacentLoops[loop4]);
			} else {
				newLoop[1]->addAdjacent(adjacentLoops[loop4]);
				newLoop[1]->addAdjacent(newLoop[0]);
			}
			adjacentLoops[loop4]->replaceAdjacent(this, newLoop[1]);
			adjacentLoops[loop4]->generateMoves();

			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}

		// MOVE_3
		if (move->type & MOVE_3) {
			//non-adjacent sides, multi + open loop

			ptypes = new int[loop4 - loop3 + 1];
			sidelengths = new int[loop4 - loop3 + 1];
			sequences = new char *[loop4 - loop3 + 1];

			pt = pairtypes[seqs[loop3][loop]][seqs[loop4][loop2]];

			for (temploop = 0, tempindex = 0; temploop < (loop4 - loop3 + 1); tempindex++) // note that loop4 - loop3 is the number of pairings that got included in the multiloop. The extra closing pair makes the +1.
					{
				if (tempindex == loop3) {
					ptypes[temploop] = pt;
					sidelengths[temploop] = sidelen[tempindex] - loop;
					sequences[temploop] = &seqs[tempindex][loop];
					temploop++;
				}
				if (tempindex > loop3 && tempindex < loop4) {
					ptypes[temploop] = pairtype[tempindex];
					sidelengths[temploop] = sidelen[tempindex];
					sequences[temploop] = seqs[tempindex];
					temploop++;
				}
				if (tempindex == loop4) {
					ptypes[temploop] = pairtype[tempindex];
					sidelengths[temploop] = loop2 - 1;
					sequences[temploop] = seqs[tempindex];
					temploop++;
				}
			}

			newLoop[0] = new MultiLoop(loop4 - loop3 + 1, ptypes, sidelengths, sequences);

			ptypes = new int[numAdjacent - (loop4 - loop3 - 1)];
			sidelengths = new int[numAdjacent - (loop4 - loop3 - 1)];
			sequences = new char *[numAdjacent - (loop4 - loop3 - 1)];

			for (temploop = 0, tempindex = 0; temploop < numAdjacent - (loop4 - loop3 - 1); tempindex++) {
				if (tempindex == loop3) {
					ptypes[temploop] = pairtype[tempindex];
					sidelengths[temploop] = loop - 1;
					sequences[temploop] = seqs[tempindex];
					temploop++;
				} else if (tempindex == loop4) {
					ptypes[temploop] = pt;
					sidelengths[temploop] = sidelen[tempindex] - loop2;
					sequences[temploop] = &seqs[tempindex][loop2];
					temploop++;
				} else if (!((tempindex > loop3) && (tempindex < loop4))) {
					ptypes[temploop] = pairtype[tempindex];
					sidelengths[temploop] = sidelen[tempindex];
					sequences[temploop] = seqs[tempindex];
					temploop++;
				}
			}
			newLoop[1] = new MultiLoop(numAdjacent - (loop4 - loop3 - 1), ptypes, sidelengths, sequences);

			// fix all the connections
			//
			for (temploop = 0; temploop < numAdjacent; temploop++) {
				if (temploop < loop3) {
					newLoop[1]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[1]);
					adjacentLoops[temploop]->generateMoves();
				} else if (temploop == loop3) {
					newLoop[1]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[1]);
					adjacentLoops[temploop]->generateMoves();
					newLoop[1]->addAdjacent(newLoop[0]);
					newLoop[0]->addAdjacent(newLoop[1]);

				} else if (temploop <= loop4) {
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
				} else {
					newLoop[1]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[1]);
					adjacentLoops[temploop]->generateMoves();
				}
			}

			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}

	}
	return -totalRate;
}

void MultiLoop::generateMoves(void) {
	int loop, loop2, loop3, loop4, temploop, tempindex, loops[4];
	int pt;
	double tempRate;
	RateEnv rateEnv;
	double energies[2];

	if (moves != NULL)
		delete moves;
	moves = new MoveList(sidelen[0] + 1);
// This is almost identical to OpenLoop::generateMoves, which was written first.
//  Several options here:
//     #1: creation move within a side this results in a hairpin and a multi loop with 1 greater magnitude.
//     #2a: creation move between sides resulting in a stack and multi loop
//     #2b: creation move between sides resulting in a bulge and multi loop
//     #2c: creation move between sides resulting in a interior loop and multi loop
//     #3: creation move between sides resulting in two multiloops.
// #2a-#2c can only happen for adjacent sides, #3 only happens for non-adjacent sides (and is always the case for such). We separate these into cases #2a-#2c and #3 .

// these pointers are needed to set up all of the multiloop energy calls.
	int *ptypes = NULL;
	int *sideLengths = NULL;
	char **sequences = NULL;

// the most storage we'll need is for case #1, which will have a multiloop of 1 greater magnitude.
	ptypes = new int[numAdjacent + 1];
	sideLengths = new int[numAdjacent + 1];
	sequences = new char *[numAdjacent + 1];

// Case #1: Single Side only Creation Moves
	for (loop3 = 0; loop3 < numAdjacent; loop3++) {
		for (loop = 1; loop <= sidelen[loop3] - 4; loop++) {
			for (loop2 = loop + 4; loop2 <= sidelen[loop3]; loop2++) { // each possibility is a hairpin and multiloop, see above.

				//FD: loop3 is the strand that will split.
				//FD: loop and loop2 are the nucleotide indices.
				//FD: Loop2 - loop is at least 4, e.g. this is the hairpin length.
				//FD: the length of the right-side remaining loop is sidelen[loop3]-loop2;

				pt = pairtypes[seqs[loop3][loop]][seqs[loop3][loop2]];

				if (pt != 0) {

					energies[0] = energyModel->HairpinEnergy(&seqs[loop3][loop], loop2 - loop - 1);

					for (temploop = 0, tempindex = 0; temploop < numAdjacent + 1; temploop++, tempindex++) {
						if (temploop == loop3) {

							ptypes[temploop] = pairtype[temploop];
							sideLengths[temploop] = loop - 1;
							sequences[temploop] = seqs[temploop];
							ptypes[temploop + 1] = pt;
							sideLengths[temploop + 1] = sidelen[temploop] - loop2;
							sequences[temploop + 1] = &seqs[temploop][loop2];

							// FD: This places an additional side to the multiloop.
							// FD: The left-side retains loop3 location, the right-side is now indexed at loop3+1.

							temploop = temploop + 1;

						} else {

							ptypes[temploop] = pairtype[tempindex];
							sideLengths[temploop] = sidelen[tempindex];
							sequences[temploop] = seqs[tempindex];

						}
					}
					energies[1] = energyModel->MultiloopEnergy(numAdjacent + 1, sideLengths, sequences);

					tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

					// multiLoop is closing, so this an loopMove and something else
//					if (energyModel->useArrhenius()) {

					MoveType rightMove = energyModel->prefactorInternal(sideLengths[loop3], sideLengths[loop3]);

//					}
					rateEnv = RateEnv(tempRate, energyModel, loopMove, rightMove);

					moves->addMove(new Move(MOVE_CREATE | MOVE_1, rateEnv, this, loop, loop2, loop3));
				}
			}
		}
	}

// Case #2a-c: adjacent loop creation moves
	for (loop3 = 0; loop3 <= numAdjacent - 1; loop3++) { // CHECK: is numAdjacent really correct? it could be numAdjacent+1
		for (loop = 1; loop <= sidelen[loop3]; loop++) {
			for (loop2 = 1; loop2 <= sidelen[(loop3 + 1) % numAdjacent]; loop2++) { // each possibility is a hairpin and open loop, see above.
				loop4 = (loop3 + 1) % numAdjacent;

				pt = pairtypes[seqs[loop3][loop]][seqs[loop4][loop2]];

				if (pt != 0) {

					MoveType leftMove = stackMove;

					// three cases for which type of move:
					// #2a: stack
					if (loop == sidelen[loop3] && loop2 == 1) {

						energies[0] = energyModel->StackEnergy(seqs[loop3][loop], seqs[loop4][loop2], seqs[loop3][sidelen[loop3] + 1], seqs[loop4][0]);

					} else if (loop == sidelen[loop3] || loop2 == 1) { 			// #2b: bulge

						if (loop2 == 1) {
							energies[0] = energyModel->BulgeEnergy(seqs[loop3][loop], seqs[loop4][loop2], seqs[loop3][sidelen[loop3] + 1], seqs[loop4][0],
									sidelen[loop3] - loop);
						} else {
							energies[0] = energyModel->BulgeEnergy(seqs[loop3][loop], seqs[loop4][loop2], seqs[loop3][sidelen[loop3] + 1], seqs[loop4][0],
									loop2 - 1);
						}

						leftMove = stackLoopMove;

					} else {				 					// #2c: interior

						energies[0] = energyModel->InteriorEnergy(&seqs[loop3][loop], seqs[loop4], sidelen[loop3] - loop, loop2 - 1);
						leftMove = loopMove;

					}

					for (temploop = 0; temploop < numAdjacent; temploop++) {
						if (temploop == loop3) {
							//FD: This is computing the sideLengths for the remaining loops.

							ptypes[temploop] = pairtype[temploop];
							sideLengths[temploop] = loop - 1;
							sequences[temploop] = seqs[temploop];

						} else {

							ptypes[temploop] = pairtype[temploop];
							if (temploop == loop4) {
								ptypes[temploop] = pt;
								sideLengths[temploop] = sidelen[temploop] - loop2;
								sequences[temploop] = &seqs[temploop][loop2];
							} else {
								sideLengths[temploop] = sidelen[temploop];
								sequences[temploop] = seqs[temploop];
							}
						}
					}
					energies[1] = energyModel->MultiloopEnergy(numAdjacent, sideLengths, sequences);
					tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

					// multiLoop is forming an stack/bulge/interior, which is something and something else
					MoveType rightMove = energyModel->prefactorInternal(sideLengths[loop3], sideLengths[loop4]);

					rateEnv = RateEnv(tempRate, energyModel, leftMove, rightMove);
					moves->addMove(new Move(MOVE_CREATE | MOVE_2, rateEnv, this, loop, loop2, loop3));
				}
			}
		}
	}

// FIXME 01/17/05 JS - This appears to be directly copied from the openloop section - the second half always refers to openloops, not multi as should be expected... Needs checking. RESOLVED 01/17/05 This is now updated to be exactly appropriate to the multiloop case.

// Case #3: non-adjacent loop creation moves (2d)
// Revamped so it actually works. Algorithm follows:
// This is all connections between non-adjacent sides. Thus we must exclude adjacent sides, and must try all possible combinations which match. This means we have to cover ~n^2 side combinations, where n is the total number of sides. Note that in this data structure, n is numAdjacent+1, and they are labelled 0,1,...,numAdjacent
// Loop over all sides. Within this loop, cover all sides that are labelled higher that the first, and are non adjacent. For each pair of bases in these two sides, check whether they can pair. For each pair, compute energies and add move to list.
	for (loop3 = 0; loop3 < numAdjacent - 2; loop3++) { // The last 2 entries are not needed as neither have higher numbered non-adjacent sections.

		for (loop4 = loop3 + 2; (loop4 < numAdjacent) && (loop3 != 0 || loop4 < numAdjacent - 1); loop4++) {

			// FD: Loop3, loop4 are the selected strands that will form a new base-pair.

			for (loop = 1; loop <= sidelen[loop3]; loop++) {

				for (loop2 = 1; loop2 <= sidelen[loop4]; loop2++) {

					pt = pairtypes[seqs[loop3][loop]][seqs[loop4][loop2]];

					if (pt != 0) { // result is a multiloop and multi loop.
								   // Multiloop

						for (temploop = 0, tempindex = 0; temploop < (loop4 - loop3 + 1); tempindex++) // note that loop4 - loop3 is the number of pairings that got included in the multiloop. The extra closing pair makes the +1.
								{
							if (tempindex == loop3) {
								ptypes[temploop] = pt;
								sideLengths[temploop] = sidelen[tempindex] - loop;
								sequences[temploop] = &seqs[tempindex][loop];
								temploop++;
							}

							if (tempindex > loop3 && tempindex < loop4) {
								ptypes[temploop] = pairtype[tempindex];
								sideLengths[temploop] = sidelen[tempindex];
								sequences[temploop] = seqs[tempindex];
								temploop++;
							}

							if (tempindex == loop4) {
								ptypes[temploop] = pairtype[tempindex];
								sideLengths[temploop] = loop2 - 1;
								sequences[temploop] = seqs[tempindex];
								temploop++;
							}
						}

						energies[0] = energyModel->MultiloopEnergy(loop4 - loop3 + 1, sideLengths, sequences);
						MoveType leftMove = energyModel->prefactorInternal(sideLengths[loop3], sideLengths[loop4]);

						// Multi loop
						for (temploop = 0, tempindex = 0; temploop < numAdjacent - (loop4 - loop3 - 1); tempindex++) {
							if (tempindex == loop3) {
								ptypes[temploop] = pairtype[tempindex];
								sideLengths[temploop] = loop - 1;
								sequences[temploop] = seqs[tempindex];
								temploop++;
							} else if (tempindex == loop4) {
								ptypes[temploop] = pt;
								sideLengths[temploop] = sidelen[tempindex] - loop2;
								sequences[temploop] = &seqs[tempindex][loop2];
								temploop++;
							} else if (!((tempindex > loop3) && (tempindex < loop4))) {
								ptypes[temploop] = pairtype[tempindex];
								sideLengths[temploop] = sidelen[tempindex];
								sequences[temploop] = seqs[tempindex];
								temploop++;
							}

						}
						energies[1] = energyModel->MultiloopEnergy(numAdjacent - (loop4 - loop3 - 1), sideLengths, sequences);
						tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);
						loops[0] = loop;
						loops[1] = loop2;
						loops[2] = loop3;
						loops[3] = loop4;

						// multiLoop is splitting into two multiLoops. Which is something, and something else

						MoveType rightMove = energyModel->prefactorInternal(sideLengths[loop3], sideLengths[loop4]);

						rateEnv = RateEnv(tempRate, energyModel, leftMove, rightMove);
						moves->addMove(new Move(MOVE_CREATE | MOVE_3, rateEnv, this, loops));
					}

				}

			}

		}

	}

	totalRate = moves->getRate();
	if (ptypes != NULL)
		delete[] ptypes;
	if (sideLengths != NULL)
		delete[] sideLengths;
	if (sequences != NULL)
		delete[] sequences;

	generateDeleteMoves();
}

void MultiLoop::generateDeleteMoves(void) {
	double temprate;

	assert(moves != NULL);

	for (int loop = 0; loop < numAdjacent; loop++) {

		generateAndSaveDeleteMove(adjacentLoops[loop], loop);

	}

	totalRate = moves->getRate();
}

void MultiLoop::printMove(Loop *comefrom, char *structure_p, char *seq_p) {
	int loop, item, loop2;
	for (loop = 0; loop < numAdjacent; loop++) {
		loop2 = (loop + 1) % numAdjacent;
		item = (seqs[loop] < seqs[loop2]);
		structure_p[seqs[loop] - seq_p + 1 + sidelen[loop]] = item ? '(' : ')';
		structure_p[seqs[loop2] - seq_p] = item ? ')' : '(';
	}

	for (loop = 0; loop < curAdjacent; loop++) {
		if (adjacentLoops[loop] != comefrom && adjacentLoops[loop] != NULL)
			// shouldn't happen, being careful.
			adjacentLoops[loop]->printMove(this, structure_p, seq_p);
		assert(adjacentLoops[loop] != NULL);
	}
}

char *MultiLoop::getLocation(Move *move, int index) {
	if (move->getType() & MOVE_CREATE) {
		if (move->getType() & MOVE_1)
			return &seqs[move->index[2]][move->index[index]];
		if (move->getType() & MOVE_2)
			return &seqs[(move->index[2] + index) % numAdjacent][move->index[index]];
		if (move->getType() & MOVE_3)
			return &seqs[move->index[2 + index]][move->index[index]];
	} else if (move->getType() & MOVE_DELETE) {
		for (int loop = 0; loop < numAdjacent; loop++)
			if (adjacentLoops[loop] == move->affected[0] || adjacentLoops[loop] == move->affected[1])
				return seqs[loop];
	}
	assert(0);
}

char *MultiLoop::verifyLoop(char *incoming_sequence, int incoming_pairtype, Loop *from) {
	char *ret_seq;
	int call_index, call_adjacent;
	int adjacent;
	int loop;
	for (loop = 0, adjacent = numAdjacent - 1; loop < numAdjacent; loop++, adjacent++) {
		if (adjacent == numAdjacent)
			adjacent = 0;

		if (adjacentLoops[loop] == from) {
			call_index = loop;
			call_adjacent = adjacent;
		} else {
			ret_seq = adjacentLoops[loop]->verifyLoop(&seqs[adjacent][sidelen[adjacent] + 1], pairtype[loop], this);
			assert(ret_seq == seqs[loop]);
		}
	}

	if (incoming_sequence != seqs[call_index]) {
		fprintf(stderr, "Verification Failed\n");
		assert(incoming_sequence == seqs[call_index]);
	}
	return seqs[call_adjacent] + sidelen[call_adjacent] + 1;
}

/* OpenLoop functions */
OpenLoop::OpenLoop(void) {
	numAdjacent = 0;
	adjacentLoops = NULL;

	pairtype = NULL;
	sidelen = NULL;
	seqs = NULL;
	identity = 'O';

}

OpenLoop::~OpenLoop(void) {
	delete[] pairtype;
	delete[] sidelen;
	delete[] seqs;
}

OpenLoop::OpenLoop(int branches, int *pairtypes, int *sidelengths, char **sequences) {
	numAdjacent = branches;
	if (branches > 0) {
		adjacentLoops = new Loop *[branches];
		for (int loop = 0; loop < branches; loop++)
			adjacentLoops[loop] = NULL;
	} else
		adjacentLoops = NULL;

	pairtype = pairtypes;
	sidelen = sidelengths;
	seqs = sequences;
	identity = 'O';
}

string OpenLoop::typeInternalsToString(void) {

	std::stringstream ss;

	for (int i = 0; i < numAdjacent + 1; i++) {

		ss << utility::sequenceToString(seqs[i], sidelen[i]) << "  -- ";

	}

	for (int i = 0; i < numAdjacent; i++) {

		ss << "        pairTypes " << basepairString[pairtype[i]] << ", ";

	}

	ss << " \n";

	ss << openInfo;

	return ss.str();

}

//void OpenLoop::halfContextToString(std::stringstream& ss) {
//
//	ss << context;
//
//}

// if the half contexts are updated, print them

//	if (updatedContext) {
//
//		ss << "hContext: \n";
//
//		for (int i = 0; i < numAdjacent + 1; i++) {
//
//			ss << "   c" << i << " - ";
//
//			for (int j = 0; j < context[i].size(); j++) {
//
//				ss << context[i][j];
//
//			}
//
//			ss << "  \n";
//
//		}
//
//	}

//	ss << " \n";

//	return ss.str();

void OpenLoop::calculateEnergy(void) {
	if (energyModel == NULL)
		return; // if the loops try to get used before the energy model initializes, it's all over.

	energy = energyModel->OpenloopEnergy(numAdjacent, sidelen, seqs);
	return;
}

Move *OpenLoop::getChoice(double *randomchoice, Loop *from) {
	Move *stor;
	assert(randomchoice != NULL);
	assert(*randomchoice >= 0.0); // never should see a negative choice value.

	if (*randomchoice < totalRate) { // something was chosen, do this
		return moves->getChoice(randomchoice);
	} else {
		*randomchoice -= totalRate;
		for (int loop = 0; loop < curAdjacent; loop++) {
			if (adjacentLoops[loop] != from) {
				stor = adjacentLoops[loop]->getChoice(randomchoice, this);
				if (stor != NULL)
					return stor;
			}
		}
	}
	return NULL;
}

double OpenLoop::doChoice(Move *move, Loop **returnLoop) {
	Loop *newLoop[2];
	int pt, loop, loop2, loop3, loop4, temploop, tempindex;
	int *ptypes;
	int *sidelengths;
	char **sequences;

	if (move->type & MOVE_CREATE) {
		loop = move->index[0];
		loop2 = move->index[1];
		loop3 = move->index[2];
		loop4 = move->index[3];

		if (move->type & MOVE_1) {
			//single side, hairpin + open with 1 higher mag.
			ptypes = new int[numAdjacent + 1];
			sidelengths = new int[numAdjacent + 2];
			sequences = new char *[numAdjacent + 2];

			pt = pairtypes[seqs[loop3][loop]][seqs[loop3][loop2]];

			for (temploop = 0, tempindex = 0; temploop <= numAdjacent + 1; temploop++, tempindex++) {
				if (temploop == loop3) {
					ptypes[temploop] = pt;
					sidelengths[temploop] = loop - 1;
					sequences[temploop] = seqs[temploop];
					if (temploop != numAdjacent)
						ptypes[temploop + 1] = pairtype[temploop];
					sidelengths[temploop + 1] = sidelen[temploop] - loop2;
					sequences[temploop + 1] = seqs[temploop] + loop2;
					temploop = temploop + 1;
				} else {
					if (temploop != numAdjacent + 1)
						ptypes[temploop] = pairtype[tempindex];
					sidelengths[temploop] = sidelen[tempindex];
					sequences[temploop] = seqs[tempindex];
				}
			}

			newLoop[0] = new OpenLoop(numAdjacent + 1, ptypes, sidelengths, sequences);

			newLoop[1] = new HairpinLoop(pt, loop2 - loop - 1, &seqs[loop3][loop]);

			for (temploop = 0; temploop < numAdjacent; temploop++) {
				if (temploop == loop3) {
					newLoop[0]->addAdjacent(newLoop[1]);
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
				} else {
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
				}
			}
			if (temploop == loop3)
				newLoop[0]->addAdjacent(newLoop[1]);

			newLoop[1]->addAdjacent(newLoop[0]);
			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}

		// MOVE_2
		if (move->type & MOVE_2) {
			//adjacent sides, one of: stack    + open with same mag
			//                        bulge    + open with same mag
			//                        interior + open with same mag
			ptypes = new int[numAdjacent];
			sidelengths = new int[numAdjacent + 1];
			sequences = new char *[numAdjacent + 1];

			pt = pairtypes[seqs[loop3][loop]][seqs[loop3 + 1][loop2]];

			for (temploop = 0; temploop <= numAdjacent; temploop++) {
				if (temploop == loop3) {
					ptypes[temploop] = pt;
					sidelengths[temploop] = loop - 1;
					sequences[temploop] = seqs[temploop];
				} else {
					if (temploop != numAdjacent)
						ptypes[temploop] = pairtype[temploop];
					if (temploop == loop3 + 1) {
						sidelengths[temploop] = sidelen[temploop] - loop2;
						sequences[temploop] = &seqs[temploop][loop2];
					} else {
						sidelengths[temploop] = sidelen[temploop];
						sequences[temploop] = seqs[temploop];
					}
				}
			}

			newLoop[0] = new OpenLoop(numAdjacent, ptypes, sidelengths, sequences);

			// three cases for which type of move:
			// #2a: stack
			if (loop == sidelen[loop3] && loop2 == 1) {
				newLoop[1] = new StackLoop(pt, pairtype[loop3], &seqs[loop3][loop], &seqs[loop3 + 1][0]);
			}
			// #2b: bulge
			else if (loop == sidelen[loop3] || loop2 == 1) {
				newLoop[1] = new BulgeLoop(pt, pairtype[loop3], sidelen[loop3] - loop, loop2 - 1, &seqs[loop3][loop], &seqs[loop3 + 1][0]);
			}

			//FIXME 01/17/05: not necessarily in this location: need to make sure that we have a consistent case: sequences in multiloops/openloops are always /before/ (as is the case with open loops) or after, the pairing. 01/17/05 - this definately is the case, open loops have the sequence /before/ the pair with the same index - multi and all others have the squence /after/ the pair with the same index. It appears to be consistently used in most cases, but perhaps i should add an assert into the code to ensure this is the case.

			// #2c: interior
			else {
				newLoop[1] = new InteriorLoop(pt, pairtype[loop3], sidelen[loop3] - loop, loop2 - 1, &seqs[loop3][loop], &seqs[loop3 + 1][0]);
			}

			for (temploop = 0; temploop < numAdjacent; temploop++) {
				if (temploop == loop3) {
					newLoop[0]->addAdjacent(newLoop[1]);
				} else {
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
				}
			}

			newLoop[1]->addAdjacent(newLoop[0]);
			newLoop[1]->addAdjacent(adjacentLoops[loop3]);
			adjacentLoops[loop3]->replaceAdjacent(this, newLoop[1]);
			adjacentLoops[loop3]->generateMoves();

			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			*returnLoop = newLoop[0];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}

		// MOVE_3
		if (move->type & MOVE_3) {
			//non-adjacent sides, multi + open loop

			ptypes = new int[loop4 - loop3 + 1];
			sidelengths = new int[loop4 - loop3 + 1];
			sequences = new char *[loop4 - loop3 + 1];

			pt = pairtypes[seqs[loop3][loop]][seqs[loop4][loop2]];

			for (temploop = 0, tempindex = 0; temploop < (loop4 - loop3 + 1); tempindex++) // note that loop4 - loop3 is the number of pairings that got included in the multiloop. The extra closing pair makes the +1.
					{
				if (tempindex == loop3) {
					ptypes[temploop] = pt;
					sidelengths[temploop] = sidelen[tempindex] - loop;
					sequences[temploop] = &seqs[tempindex][loop];
					temploop++;
				}
				if (tempindex > loop3 && tempindex < loop4) {
					ptypes[temploop] = pairtype[tempindex - 1];
					sidelengths[temploop] = sidelen[tempindex];
					sequences[temploop] = seqs[tempindex];
					temploop++;
				}
				if (tempindex == loop4) {
					ptypes[temploop] = pairtype[tempindex - 1];
					sidelengths[temploop] = loop2 - 1;
					sequences[temploop] = seqs[tempindex];
					temploop++;
				}
			}

			newLoop[0] = new MultiLoop(loop4 - loop3 + 1, ptypes, sidelengths, sequences);

			ptypes = new int[numAdjacent - (loop4 - loop3 - 1)];
			sidelengths = new int[numAdjacent - (loop4 - loop3 - 1) + 1];
			sequences = new char *[numAdjacent - (loop4 - loop3 - 1) + 1];

			for (temploop = 0, tempindex = 0; temploop <= numAdjacent - (loop4 - loop3 - 1); tempindex++) {
				if (tempindex == loop3) {
					ptypes[temploop] = pt;
					sidelengths[temploop] = loop - 1;
					sequences[temploop] = seqs[tempindex];
					temploop++;
				} else if (tempindex == loop4) {
					if (temploop < numAdjacent - (loop4 - loop3 - 1))
						ptypes[temploop] = pairtype[tempindex];
					sidelengths[temploop] = sidelen[tempindex] - loop2;
					sequences[temploop] = &seqs[tempindex][loop2];
					temploop++;
				} else if (!((tempindex > loop3) && (tempindex < loop4))) {
					if (temploop != numAdjacent - (loop4 - loop3 - 1))
						ptypes[temploop] = pairtype[tempindex];
					sidelengths[temploop] = sidelen[tempindex];
					sequences[temploop] = seqs[tempindex];
					temploop++;
				}
			}
			newLoop[1] = new OpenLoop(numAdjacent - (loop4 - loop3 - 1), ptypes, sidelengths, sequences);

			// fix all the connections

			for (temploop = 0; temploop < numAdjacent; temploop++) {
				if (temploop < loop3) {
					newLoop[1]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[1]);
					adjacentLoops[temploop]->generateMoves();
				} else if (temploop == loop3) {
					newLoop[1]->addAdjacent(newLoop[0]);
					newLoop[0]->addAdjacent(newLoop[1]);
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
				} else if (temploop < loop4) {
					newLoop[0]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[0]);
					adjacentLoops[temploop]->generateMoves();
				} else {
					newLoop[1]->addAdjacent(adjacentLoops[temploop]);
					adjacentLoops[temploop]->replaceAdjacent(this, newLoop[1]);
					adjacentLoops[temploop]->generateMoves();
				}
			}

			newLoop[0]->generateMoves();
			newLoop[1]->generateMoves();
			*returnLoop = newLoop[1];
			return ((newLoop[0]->getTotalRate() + newLoop[1]->getTotalRate()) - totalRate);
		}

	}
	return -totalRate;
}

void OpenLoop::generateMoves(void) {

	int loop, loop2, loop3, loop4, temploop, tempindex, loops[4];
	int pairType;
	double tempRate;
	RateEnv rateEnv;
	double energies[2];

	if (moves != NULL)
		delete moves;
	moves = new MoveList(1);
//  Several options here:
//     #1: creation move within a side this results in a hairpin and a open loop with 1 greater magnitude.
//     #2a: creation move between sides resulting in a stack and open loop
//     #2b: creation move between sides resulting in a bulge and open loop
//     #2c: creation move between sides resulting in a interior loop and open loop
//     #2d: creation move between sides resulting in a multiloop and open loop. (ICK).
// #2a-#2c can only happen for adjacent sides, #2d only happens for non-adjacent sides (and is always the case for such). We separate these into cases #2 (#2a-#2c) and #3 (#2d).

// these three pointers are needed to set up the open loop's energy calls.
// i'd like to optimize so they don't need to be created/deleted very often
// but i'm  not sure of a good way of handling that yet.

	int *ptypes = NULL;
	int *sideLengths = NULL;
	char **sequences = NULL;

	ptypes = new int[numAdjacent + 1];
	sideLengths = new int[numAdjacent + 2];
	sequences = new char *[numAdjacent + 2];
// Case #1: Single Side only Creation Moves
	for (loop3 = 0; loop3 < numAdjacent + 1; loop3++) {

		char* mySequence = seqs[loop3]; // this is the sequence of the strand that we use

		for (loop = 1; loop < sidelen[loop3] - 3; loop++) {

			for (loop2 = loop + 4; loop2 <= sidelen[loop3]; loop2++) { // each possibility is a hairpin and open loop, see above.

				pairType = pairtypes[mySequence[loop]][mySequence[loop2]];

				if (pairType != 0) { // FD: the NUPACK model puts terms here to be non-zero.    in NUPACK, G-T stacking is a thing. Hairpin loops are size 3 or more.

					energies[0] = energyModel->HairpinEnergy(&mySequence[loop], loop2 - loop - 1);

					for (temploop = 0, tempindex = 0; temploop < numAdjacent + 2; temploop++, tempindex++) {
						if (temploop == loop3) {
							ptypes[temploop] = pairType;
							sideLengths[temploop] = loop - 1;
							sequences[temploop] = seqs[temploop];
							if (temploop != numAdjacent)
								ptypes[temploop + 1] = pairtype[temploop];
							sideLengths[temploop + 1] = sidelen[temploop] - loop2;
							sequences[temploop + 1] = seqs[temploop] + loop2;
							temploop = temploop + 1;
						} else {
							if (temploop != numAdjacent + 1)
								ptypes[temploop] = pairtype[tempindex];
							sideLengths[temploop] = sidelen[tempindex];
							sequences[temploop] = seqs[tempindex];
						}
					}
					energies[1] = energyModel->OpenloopEnergy(numAdjacent + 1, sideLengths, sequences);
					tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

					// if the new Arrhenius model is used, modify the existing rate based on the local context.
					// to start, we need to learn what the local context is, AFTER the nucleotide is put in place.

					// OpenLoop is splitting off an hairpin. Which is loopMove, and something else

					MoveType rightMove = energyModel->prefactorOpen(loop3, numAdjacent + 2, sideLengths);
					rateEnv = RateEnv(tempRate, energyModel, loopMove, rightMove);

					Move *tmove = new Move(MOVE_CREATE | MOVE_1, rateEnv, this, loop, loop2, loop3);
					moves->addMove(tmove);
				}
			}
		}
	}

// Case #2a-c: adjacent loop creation moves
	for (loop3 = 0; loop3 < numAdjacent; loop3++) // CHECK: is numAdjacent really correct? it could be numAdjacent+1
		for (loop = 1; loop <= sidelen[loop3]; loop++)
			for (loop2 = 1; loop2 <= sidelen[loop3 + 1]; loop2++) { // each possibility is a hairpin and open loop, see above.

				pairType = pairtypes[seqs[loop3][loop]][seqs[loop3 + 1][loop2]];

				if (pairType != 0) {

					// three cases for which type of move:
					MoveType leftMove = stackMove;

					if (loop == sidelen[loop3] && loop2 == 1) { 					// #2a: stack

						energies[0] = energyModel->StackEnergy(seqs[loop3][loop], seqs[loop3 + 1][loop2], seqs[loop3][sidelen[loop3] + 1], seqs[loop3 + 1][0]);

					} else if (loop == sidelen[loop3] || loop2 == 1) { 			// #2b: bulge

						if (loop2 == 1) {

							energies[0] = energyModel->BulgeEnergy(seqs[loop3][loop], seqs[loop3 + 1][loop2], seqs[loop3][sidelen[loop3] + 1],
									seqs[loop3 + 1][0], sidelen[loop3] - loop);

						} else {

							energies[0] = energyModel->BulgeEnergy(seqs[loop3][loop], seqs[loop3 + 1][loop2], seqs[loop3][sidelen[loop3] + 1],
									seqs[loop3 + 1][0], loop2 - 1);

						}

						leftMove = stackLoopMove;

					} else { 					// #2c: interior

						energies[0] = energyModel->InteriorEnergy(&seqs[loop3][loop], seqs[loop3 + 1], sidelen[loop3] - loop, loop2 - 1);

						leftMove = loopMove;

					}

					for (temploop = 0; temploop < numAdjacent + 1; temploop++) {
						if (temploop == loop3) {
							ptypes[temploop] = pairType;
							sideLengths[temploop] = loop - 1;
							sequences[temploop] = seqs[temploop];
						} else {
							if (temploop != numAdjacent)
								ptypes[temploop] = pairtype[temploop];
							if (temploop == loop3 + 1) {
								sideLengths[temploop] = sidelen[temploop] - loop2;
								sequences[temploop] = &seqs[temploop][loop2];
							} else {
								sideLengths[temploop] = sidelen[temploop];
								sequences[temploop] = seqs[temploop];
							}
						}
					}
					energies[1] = energyModel->OpenloopEnergy(numAdjacent, sideLengths, sequences);

					tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

					// openLoop is splitting off an stack/bulge/interior, and another openloop.
					// Which is something, and something else

					// the new stack/bulge/interior is the LeftMove (see above);
					// the new Openloop:
					MoveType rightMove = energyModel->prefactorOpen(loop3, numAdjacent + 1, sideLengths);

					rateEnv = RateEnv(tempRate, energyModel, leftMove, rightMove);
					moves->addMove(new Move(MOVE_CREATE | MOVE_2, rateEnv, this, loop, loop2, loop3));
				}
			}

// Case #3: non-adjacent loop creation moves (2d)
// Revamped so it actually works. Algorithm follows:
// This is all connections between non-adjacent sides. Thus we must exclude adjacent sides, and must try all possible combinations which match. This means we have to cover ~n^2 side combinations, where n is the total number of sides. Note that in this data structure, n is numAdjacent+1, and they are labelled 0,1,...,numAdjacent
// Loop over all sides. Within this loop, cover all sides that are labelled higher that the first, and are non adjacent. For each pair of bases in these two sides, check whether they can pair. For each pair, compute energies and add move to list.
	for (loop3 = 0; loop3 <= numAdjacent - 2; loop3++) // The last 2 entries are not needed as neither have higher numbered non-adjacent sections.
		for (loop4 = loop3 + 2; loop4 <= numAdjacent; loop4++) {

			for (loop = 1; loop <= sidelen[loop3]; loop++) { // new version with all sequences in openloop starting at 1.

				for (loop2 = 1; loop2 <= sidelen[loop4]; loop2++) {

					pairType = pairtypes[seqs[loop3][loop]][seqs[loop4][loop2]];

					if (pairType != 0) { // result is a multiloop and open loop.

						for (temploop = 0, tempindex = 0; temploop < (loop4 - loop3 + 1); tempindex++) { // note that loop4 - loop3 is the number of pairings that got included in the multiloop. The extra closing pair makes the +1.

							if (tempindex == loop3) {
								ptypes[temploop] = pairType;
								sideLengths[temploop] = sidelen[tempindex] - loop;
								sequences[temploop] = &seqs[tempindex][loop];
								temploop++;
							}
							if (tempindex > loop3 && tempindex < loop4) {
								ptypes[temploop] = pairtype[tempindex - 1];
								sideLengths[temploop] = sidelen[tempindex];
								sequences[temploop] = seqs[tempindex];
								temploop++;
							}
							if (tempindex == loop4) {
								ptypes[temploop] = pairtype[tempindex - 1];
								sideLengths[temploop] = loop2 - 1;
								sequences[temploop] = seqs[tempindex];
								temploop++;
							}
						}

						energies[0] = energyModel->MultiloopEnergy(loop4 - loop3 + 1, sideLengths, sequences);
						MoveType leftMove = energyModel->prefactorInternal(sideLengths[loop3], sideLengths[loop4]);

						// Open loop
						for (temploop = 0, tempindex = 0; temploop <= numAdjacent - (loop4 - loop3 - 1); tempindex++) {
							if (tempindex == loop3) {
								ptypes[temploop] = pairType;
								sideLengths[temploop] = loop - 1;
								sequences[temploop] = seqs[tempindex];
								temploop++;
							} else if (tempindex == loop4) {
								if (temploop < numAdjacent - (loop4 - loop3 - 1))
									ptypes[temploop] = pairtype[tempindex];
								sideLengths[temploop] = sidelen[tempindex] - loop2;
								sequences[temploop] = &seqs[tempindex][loop2];
								temploop++;
							} else if (!((tempindex > loop3) && (tempindex < loop4))) {
								if (temploop != numAdjacent - (loop4 - loop3 - 1))
									ptypes[temploop] = pairtype[tempindex];
								sideLengths[temploop] = sidelen[tempindex];
								sequences[temploop] = seqs[tempindex];
								temploop++;
							}
						}
						energies[1] = energyModel->OpenloopEnergy(numAdjacent - (loop4 - loop3 - 1), sideLengths, sequences);
						tempRate = energyModel->returnRate(getEnergy(), (energies[0] + energies[1]), 0);

						// openLoop is splitting off . Which is something, and something else

						MoveType rightMove = energyModel->prefactorOpen(loop3, numAdjacent - (loop4 - loop3) + 2, sideLengths);

						loops[0] = loop;
						loops[1] = loop2;
						loops[2] = loop3;
						loops[3] = loop4;

						rateEnv = RateEnv(tempRate, energyModel, leftMove, rightMove);

						moves->addMove(new Move(MOVE_CREATE | MOVE_3, rateEnv, this, loops));
					}

				}
			}
		}

	totalRate = moves->getRate();

	if (ptypes != NULL)
		delete[] ptypes;
	if (sideLengths != NULL)
		delete[] sideLengths;
	if (sequences != NULL)
		delete[] sequences;

	generateDeleteMoves();
}

void OpenLoop::generateDeleteMoves(void) {
	double temprate;

	assert(moves != NULL);

	for (int loop = 0; loop < numAdjacent; loop++) {
		generateAndSaveDeleteMove(adjacentLoops[loop], loop);

	}

	totalRate = moves->getRate();
}

void OpenLoop::printMove(Loop *comefrom, char *structure_p, char *seq_p) {

	int loop, item, loop2;

	for (loop = 0; loop < numAdjacent; loop++) {
		loop2 = (loop + 1) % (numAdjacent + 1);
		item = (seqs[loop] < seqs[loop2]);
		structure_p[seqs[loop] - seq_p + 1 + sidelen[loop]] = item ? '(' : ')';
		structure_p[seqs[loop2] - seq_p] = item ? ')' : '(';
	}

	for (loop = 0; loop < curAdjacent; loop++) {
		if (adjacentLoops[loop] != comefrom && adjacentLoops[loop] != NULL)
			// shouldn't happen, being careful.
			adjacentLoops[loop]->printMove(this, structure_p, seq_p);
		assert(adjacentLoops[loop] != NULL);
	}

}

char *OpenLoop::getLocation(Move *move, int index) {

	if (move->getType() & MOVE_CREATE) {
		if (move->getType() & MOVE_1)
			return &seqs[move->index[2]][move->index[index]];
		if (move->getType() & MOVE_2)
			return &seqs[move->index[2] + index][move->index[index]];
		if (move->getType() & MOVE_3)
			return &seqs[move->index[2 + index]][move->index[index]];
	} else if (move->getType() & MOVE_DELETE) {
		for (int loop = 0; loop < numAdjacent; loop++)
			if (adjacentLoops[loop] == move->affected[0] || adjacentLoops[loop] == move->affected[1])
				return seqs[loop + 1];
	}

	assert(0);
	return NULL;

}

//char *OpenLoop::getBase(char type, int index) {
//	int loop, loop2;
//	int newindex = index;
//	for (loop = 0; loop <= numAdjacent; loop++) {
//		for (loop2 = 1; loop2 <= sidelen[loop]; loop2++)
//			if (seqs[loop][loop2] == type) {
//				if (newindex == 0)
//					return &seqs[loop][loop2];
//				else
//					newindex--;
//			}
//	}
//
//	assert(0);
//	return NULL;
//}

char* OpenLoop::getBase(char type, int index, HalfContext half) {

	for (int loop = 0; loop <= numAdjacent; loop++) {

		int end = sidelen[loop] + 1;

		for (int loop2 = 1; loop2 < end; loop2++)

			if (seqs[loop][loop2] == type) {

				// potential match, if the halfContext matches.

				HalfContext thisHalf = getHalfContext(loop, loop2);

				if (half == thisHalf) {
					// it's a match.

					if (index == 0) {

						return &seqs[loop][loop2];

					} else {

						index--;

					}

				}

			}
	}

	if (utility::debugTraces) {

		cout << "Failing with  \n";

		cout << "type: " << (int) type << "\n";
		cout << "index: " << index << "\n";
		cout << "HalfContext: " << half << "\n";
		cout << "This OpenLoop Info: \n" << openInfo << endl;

	}

	assert(0);
	return NULL;
}

char* OpenLoop::getBase(char type, int index, bool useArr) {

	// FD 2016-11-14
	// adjusting this to work with arrhenius rates.

	for (int loop = 0; loop <= numAdjacent; loop++) {

		int loop2 = 1;
		int end = sidelen[loop] + 1;

		if (useArr) {

			loop2++;
			end--;

		}

		for (; loop2 < end; loop2++)

			if (seqs[loop][loop2] == type) {

				if (index == 0) {

					return &seqs[loop][loop2];

				} else {
					index--;
				}

			}
	}

	assert(0);
	return NULL;
}

//int *OpenLoop::getFreeBases(void) {
//	int *results;
//	int loop, loop2;
//	results = new int[5];
//	for (loop = 0; loop < 5; loop++)
//		results[loop] = 0;
//
//	for (loop = 0; loop <= numAdjacent; loop++) {
//		for (loop2 = 1; loop2 <= sidelen[loop]; loop2++) {
//			if (seqs[loop][loop2] < 5)
//				results[seqs[loop][loop2]]++;
//			else
//				results[0]++;
//		}
//	}
//	return results;
//}

// if using Arr, do not count the external bases.
BaseCount& OpenLoop::getFreeBases() {

	// do nothing if already computed last time
	if (updatedContext2) {

		return exposedBases;

	} else {

		exposedBases.clear();

		for (int loop = 0; loop < (numAdjacent + 1); loop++) {

			int loop2 = 1;
			int end = sidelen[loop] + 1;

//			if (false) {
//				// avoid counting external nucleotides
//				loop2 = loop2 + 1;
//				end = end - 1;
//			}

			for (; loop2 < end; loop2++) {

				int base = seqs[loop][loop2];

				if (base < 5) {

					exposedBases.count[base]++;

				} else {

					// I'd like the software to fail if
					// errors in the sequence exist.
					exposedBases.count[0]++;
					cout << "Encountered nonsense base in sequence.";
					assert(0);

				}
			}

		}

		updatedContext2 = true;

	}

	return exposedBases;
}

/*
 OpenLoop::performComplexJoin

 static void OpenLoop::performComplexJoin( OpenLoop **oldLoops, OpenLoop **newLoops, char *types, int *index);

 Joins the two open loops given in the array oldLoops (size 2) at the locations given by the the types/index arrays (size 2).
 Resulting open loops are placed into the newLoops array (as pointers) for the calling function to use.

 */

void OpenLoop::performComplexJoin(OpenLoop **oldLoops, OpenLoop **newLoops, char *types, int *index, HalfContext* halfs, bool useArr) {

	// FD: nov 23 2016. Need to adapt this to take into account the local context.

	int seqnum[2] = { -1, -1 };
	int seqindex[2] = { -1, -1 };
	int sizes[2];
	int loop, loop2;
	int newindex;
	int toggle;
	OpenLoop *newLoop;
	int *pairtype;
	int *sidelen;
	char **seqs;

	// FD: After this loop, SEQNUM, SEQINDEX are properly set

	for (toggle = 0; toggle <= 1; toggle++) {

		newindex = index[toggle];

		for (loop = 0; loop <= oldLoops[toggle]->numAdjacent && newindex >= 0; loop++) {

			int end = oldLoops[toggle]->sidelen[loop] + 1;

//			if (useArr) {
//
//				loop2++;
//				end--;
//
//			}

			for (int loop2 = 1; loop2 < end; loop2++) {

				HalfContext thisHalf = oldLoops[toggle]->getHalfContext(loop, loop2);

				if (!useArr || (thisHalf == halfs[toggle])) {

					if (oldLoops[toggle]->seqs[loop][loop2] == types[toggle]) {

						if (newindex == 0) {

							seqnum[toggle] = loop;
							seqindex[toggle] = loop2;
							loop2 = oldLoops[toggle]->sidelen[loop] + 1;
							newindex--;

						} else {

							newindex--;

						}
					}
				}

			}
		}

	}

// seqnum and seqindex now have the appropriate locations within each openloop. Time to compute the new #'s of adjacent loops.

	sizes[0] = seqnum[0] + 1 + (oldLoops[1]->numAdjacent - seqnum[1]);
	sizes[1] = seqnum[1] + 1 + (oldLoops[0]->numAdjacent - seqnum[0]);

	for (toggle = 0; toggle <= 1; toggle++) {
		pairtype = new int[sizes[toggle]];
		sidelen = new int[sizes[toggle] + 1];
		seqs = new char *[sizes[toggle] + 1];

		for (loop = 0; loop < sizes[toggle] + 1; loop++) {
			if (loop < seqnum[toggle]) {
				if (loop < sizes[toggle]) // should never be false?
					pairtype[loop] = oldLoops[toggle]->pairtype[loop];
				sidelen[loop] = oldLoops[toggle]->sidelen[loop];
				seqs[loop] = oldLoops[toggle]->seqs[loop];

			} else if (loop == seqnum[toggle]) {
				if (loop < sizes[toggle]) // again, never false
					pairtype[loop] = pairtypes[types[toggle]][types[1 - toggle]];

				sidelen[loop] = seqindex[toggle] - 1;
				seqs[loop] = oldLoops[toggle]->seqs[loop];
			} else if (loop > seqnum[toggle]) {
				if (loop < sizes[toggle])
					pairtype[loop] = oldLoops[1 - toggle]->pairtype[loop - seqnum[toggle] - 1 + seqnum[1 - toggle]];

				if (loop == seqnum[toggle] + 1) {
					sidelen[loop] = oldLoops[1 - toggle]->sidelen[loop - seqnum[toggle] - 1 + seqnum[1 - toggle]] - seqindex[1 - toggle];
					seqs[loop] = &(oldLoops[1 - toggle]->seqs[loop - seqnum[toggle] + seqnum[1 - toggle] - 1][seqindex[1 - toggle]]);
				} else {
					sidelen[loop] = oldLoops[1 - toggle]->sidelen[loop - seqnum[toggle] + seqnum[1 - toggle] - 1];
					seqs[loop] = oldLoops[1 - toggle]->seqs[loop - seqnum[toggle] + seqnum[1 - toggle] - 1];
				}
			}
		}
		//initialize the new openloops, and connect them correctly, then initialize their moves, etc.
		newLoop = new OpenLoop(sizes[toggle], pairtype, sidelen, seqs);

		for (loop = 0; loop < sizes[toggle]; loop++) {
			if (loop < seqnum[toggle]) {
				newLoop->addAdjacent(oldLoops[toggle]->adjacentLoops[loop]);
				loop2 = oldLoops[toggle]->adjacentLoops[loop]->replaceAdjacent(oldLoops[toggle], newLoop);
				assert(loop2 > 0);
			} else if (loop == seqnum[toggle]) {
				newLoop->addAdjacent( NULL); // This is a Trick. We'll later replace the null value with the first/second loop generated. Note that this doesn't require a replacement as in the other cases.
			} else {
				newLoop->addAdjacent(oldLoops[1 - toggle]->adjacentLoops[loop - seqnum[toggle] + seqnum[1 - toggle] - 1]);

				loop2 = oldLoops[1 - toggle]->adjacentLoops[loop - seqnum[toggle] + seqnum[1 - toggle] - 1]->replaceAdjacent(oldLoops[1 - toggle], newLoop);
				assert(loop2 > 0);
			}
		}

		newLoops[toggle] = newLoop;
	}

// must do these after, otherwise the links between the two loops don't exist
	newLoops[0]->replaceAdjacent( NULL, newLoops[1]);
	newLoops[1]->replaceAdjacent( NULL, newLoops[0]);

	newLoops[0]->generateMoves();
	newLoops[1]->generateMoves();

// need to re-generate the moves for all adjacent loops.
// TODO: change this to only re-generate the deletion moves.
	for (toggle = 0; toggle <= 1; toggle++)
		for (loop = 0; loop < sizes[toggle]; loop++)
			if (loop != seqnum[toggle])
				newLoops[toggle]->adjacentLoops[loop]->generateMoves();

}

char *OpenLoop::verifyLoop(char *incoming_sequence, int incoming_pairtype, Loop *from) {
	char *ret_seq;
	int call_index = -1, call_adjacent;
	int adjacent;
	int loop;
	for (loop = 0; loop < numAdjacent; loop++) {
		if (adjacentLoops[loop] == from) {
			call_index = loop;
		} else {
			ret_seq = adjacentLoops[loop]->verifyLoop(&seqs[loop][sidelen[loop] + 1], pairtype[loop], this);
			assert(ret_seq == seqs[loop + 1]);
		}
	}

	if (call_index != -1 && incoming_sequence != seqs[call_index + 1]) {
		fprintf(stderr, "Verification Failed\n");
		assert(incoming_sequence == seqs[call_index + 1]);
	}
	if (call_index != -1)
		return seqs[call_index] + sidelen[call_index] + 1;
	else
		return NULL;
}

OpenInfo& OpenLoop::getOpenInfo(void) {

// do nothing if not required
	if (openInfo.upToDate) {

		return openInfo;

	} // else

	openInfo.clear();

	for (int i = 0; i < numAdjacent + 1; i++) {

		parseLocalContext(i);

	}

	openInfo.upToDate = true;

	return openInfo;

}

HalfContext OpenLoop::getHalfContext(int loop, int loop2) {

	char* mySeq = seqs[loop];

	// intialize for strands on both sides.
	HalfContext thisHalf = HalfContext(strandC, strandC);

	if (loop2 == 1) {
		// left could be end or stack
		thisHalf.left = moveutil::getContext(mySeq[0]);
	}

	if (loop2 == sidelen[loop]) {
		// right could be end or stack;
		thisHalf.right = moveutil::getContext(mySeq[sidelen[loop] + 1]);
	}

	return thisHalf;

}

void OpenLoop::parseLocalContext(int index) {

// FD: redoing this to save more information,
// and only record local context for exterior
// zero, one or two nucleotides

	char* mySeq = seqs[index];
	int size = sidelen[index];

	openInfo.numExposed += size;

	if (size > 2) {	 // there are internal nucleotides

		openInfo.numExposedInternal += size - 2;

	}

	if (size > 0) {

		char base = mySeq[1];

		QuartContext left = moveutil::getContext(mySeq[0]);
		QuartContext right;

		if (size == 1) { // exactly one nucleotide

			right = moveutil::getContext(mySeq[2]);

		} else {

			right = strandC;
		}

		openInfo.increment(left, base, right);

	}

// also record the second external nucleotide, if it exists
	if (size > 1) {

		int base = mySeq[size];

		QuartContext left = strandC;
		QuartContext right = moveutil::getContext(mySeq[size + 1]);

		openInfo.increment(left, base, right);

	}

// now update the internal exposed toeholds:
// the rates for these are easier to compute because they are
// loop by loop local contexts.

	if (openInfo.numExposedInternal > 0) {

		BaseCount myCount = BaseCount();

		for (int loop2 = 2; loop2 <= sidelen[index] - 1; loop2++) {

			// removing checks because I'd like to software to fail if
			// errors in the sequence exist.
			myCount.count[seqs[index][loop2]]++;

		}

		openInfo.increment(HalfContext(strandC, strandC), myCount);
	}

}
