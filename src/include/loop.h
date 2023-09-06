/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

/* Loop class header for main loop class and derivatives. */

#ifndef __LOOP_H__
#define __LOOP_H__

#include <string>
#include <vector>
#include "energymodel.h"
#include "move.h"
#include "moveutil.h"
#include "basetype.h"

using std::vector;

class EnergyOptions;

struct RateArr {

	double rate;
	MoveType left;
	MoveType right;

	RateArr() {
		rate = -1.0;
		left = stackMove;
		right = stackMove;
	}

	RateArr(double inRate, MoveType inLeft, MoveType inRight) {
		rate = inRate;
		left = inLeft;
		right = inRight;
	}
};

class Loop {
public:
	inline double getEnergy(EnergyModel *energyModel);
	inline double getEnthalpy(EnergyModel *energyModel);
	inline double getTotalRate(void);
	char getType(void);
	Loop(void);
	virtual ~Loop(void);
	virtual void calculateEnergy(EnergyModel *energyModel) = 0;
	virtual void calculateEnthalpy(EnergyModel *energyModel) = 0;
	virtual void generateMoves(EnergyModel *energyModel) = 0;
	virtual void generateDeleteMoves(EnergyModel*) = 0;
	virtual void regenerateDeleteMoves(EnergyModel *energyModel) = 0;
	virtual Move *getChoice(SimTimer& timer, Loop *from) = 0;
	virtual double doChoice(Move *move, Loop **returnLoop, EnergyModel *energyModel) = 0;
	virtual BaseType *getLocation(Move *move, int index) = 0;
	virtual BaseType *verifyLoop(BaseType *incoming_sequence, Loop *from) = 0;
	virtual string typeInternalsToString(void) = 0;
	virtual void printMove(Loop *comefrom, char *structure_p, BaseType *seq_p) = 0;
	Loop *getAdjacent(int index);
	int getCurAdjacent(void);
	int getNumAdjacent(void);
	void addAdjacent(Loop *loopToAdd);
	void initAdjacency(int index);
	int replaceAdjacent(Loop *loopToReplace, Loop *loopToReplaceWith);
	// sets adjacentLoops up to be deleted.
	void cleanupAdjacent(void);
	// returns the total energy of all loops underneath this one.
	double returnEnergies(Loop *comefrom, EnergyModel *energyModel);
	// sums the enthalpy of the contained loops.
	double returnEnthalpies(Loop *comefrom, EnergyModel *energyModel);
	// returns the total rate of all loops underneath this one.
	double returnFlux(Loop *comefrom);
	// returns the number of transitions in this loop.
	int getMoveCount(Loop *comefrom);

	void firstGen(Loop *comefrom, EnergyModel *energyModel);
	static RateArr generateDeleteMoveRate(Loop *start, Loop *end, EnergyModel *energyModel);
	static Loop *performDeleteMove(Move *move, EnergyModel *energyModel);
	static void performComplexSplit(
		Move *move, Loop **firstOpen, Loop **secondOpen, EnergyModel *energyModel);

	// FD functions
    static bool identify(Loop*, Loop* , char , char );
    static std::tuple<int,int> findExternalAdjacent(Loop*, Loop*);
    static std::pair<Loop*, Loop*> orderMyLoops(Loop*, Loop*, char);

	string toString(void);
	string toStringShort(void);
	void printAllMoves(Loop*, bool useArrhenius);
	void generateAndSaveDeleteMove(Loop*, int, EnergyModel*);

	// FD: moving private to public
	int numAdjacent;

protected:
	Loop** adjacentLoops;
	int curAdjacent;
	double energy;
	double enthalpy = 0;
	bool energyComputed = false;
	bool enthalpyComputed = false;
	double totalRate;
	MoveContainer *moves;
	char identity;
	int add_index;
};

class StackLoop: public Loop {
public:
	void calculateEnergy(EnergyModel *energyModel);
	void calculateEnthalpy(EnergyModel *energyModel);
	void generateMoves(EnergyModel *energyModel);
	void generateDeleteMoves(EnergyModel*);
	void regenerateDeleteMoves(EnergyModel *energyModel);
	Move *getChoice(SimTimer& timer, Loop *from);
	double doChoice(Move *move, Loop **returnLoop, EnergyModel *energyModel);
	void printMove(Loop *comefrom, char *structure_p, BaseType *seq_p);
	BaseType *getLocation(Move *move, int index);
	BaseType *verifyLoop(BaseType *incoming_sequence,  Loop *from);

	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end, EnergyModel *energyModel);
	friend Loop * Loop::performDeleteMove(Move *move, EnergyModel *energyModel);
	friend void Loop::performComplexSplit(
		Move *move, Loop **firstOpen, Loop **secondOpen, EnergyModel *energyModel);

	StackLoop(void);
	StackLoop(BaseType *seq1, BaseType *seq2, Loop *left = NULL, Loop *right = NULL);
	string typeInternalsToString(void);

private:
	BaseType *seqs[2];
};

class HairpinLoop: public Loop {
public:
	void calculateEnergy(EnergyModel *energyModel);
	void calculateEnthalpy(EnergyModel *energyModel);
	void generateMoves(EnergyModel *energyModel);
	void generateDeleteMoves(EnergyModel*);
	void regenerateDeleteMoves(EnergyModel *energyModel);
	Move *getChoice(SimTimer& timer, Loop *from);
	double doChoice(Move *move, Loop **returnLoop, EnergyModel *energyModel);
	void printMove(Loop *comefrom, char *structure_p, BaseType *seq_p);
	BaseType *getLocation(Move *move, int index);
	BaseType *verifyLoop(BaseType *incoming_sequence,  Loop *from);

	HairpinLoop(void);
	HairpinLoop( int size, BaseType *hairpin_sequence, Loop *previous = NULL);

	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end, EnergyModel *energyModel);
	friend Loop * Loop::performDeleteMove(Move *move, EnergyModel *energyModel);
	friend void Loop::performComplexSplit(
		Move *move, Loop **firstOpen, Loop **secondOpen, EnergyModel *energyModel);
	string typeInternalsToString(void);

private:
	int hairpinsize;
	BaseType *hairpin_seq;
};

class BulgeLoop: public Loop {
public:
	void calculateEnergy(EnergyModel *energyModel);
	void calculateEnthalpy(EnergyModel *energyModel);
	void generateMoves(EnergyModel *energyModel);
	void generateDeleteMoves(EnergyModel*);
	void regenerateDeleteMoves(EnergyModel *energyModel);
	Move *getChoice(SimTimer& timer, Loop *from);
	double doChoice(Move *move, Loop **returnLoop, EnergyModel *energyModel);
	void printMove(Loop *comefrom, char *structure_p, BaseType *seq_p);
	BaseType *getLocation(Move *move, int index);
	BaseType *verifyLoop(BaseType *incoming_sequence, Loop *from);

	BulgeLoop(void);
	BulgeLoop(int size1, int size2, BaseType *bulge_sequence1, BaseType *bulge_sequence2, Loop *left = NULL, Loop *right = NULL);

	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end, EnergyModel *energyModel);
	friend Loop * Loop::performDeleteMove(Move *move, EnergyModel *energyModel);
	friend void Loop::performComplexSplit(
		Move *move, Loop **firstOpen, Loop **secondOpen, EnergyModel *energyModel);
	string typeInternalsToString(void);

private:
	int bulgesize[2];
	BaseType *bulge_seq[2];
};

class InteriorLoop: public Loop {
public:
	void calculateEnergy(EnergyModel *energyModel);
	void calculateEnthalpy(EnergyModel *energyModel);
	void generateMoves(EnergyModel *energyModel);
	void generateDeleteMoves(EnergyModel*);
	void regenerateDeleteMoves(EnergyModel *energyModel);
	Move *getChoice(SimTimer& timer, Loop *from);
	double doChoice(Move *move, Loop **returnLoop, EnergyModel *energyModel);
	void printMove(Loop *comefrom, char *structure_p, BaseType *seq_p);
	BaseType *getLocation(Move *move, int index);
	BaseType *verifyLoop(BaseType *incoming_sequence,  Loop *from);

	InteriorLoop(void);
	InteriorLoop(int size1, int size2, BaseType *int_seq1, BaseType *int_seq2, Loop *left = NULL, Loop *right = NULL);

	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end, EnergyModel *energyModel);
	friend Loop * Loop::performDeleteMove(Move *move, EnergyModel *energyModel);
	friend void Loop::performComplexSplit(
		Move *move, Loop **firstOpen, Loop **secondOpen, EnergyModel *energyModel);
	string typeInternalsToString(void);

private:
	int sizes[2];
	BaseType *int_seq[2];
};

class MultiLoop: public Loop {
public:
	void calculateEnergy(EnergyModel *energyModel);
	void calculateEnthalpy(EnergyModel *energyModel);
	void generateMoves(EnergyModel *energyModel);
	void generateDeleteMoves(EnergyModel*);
	void regenerateDeleteMoves(EnergyModel *energyModel);
	Move *getChoice(SimTimer& timer, Loop *from);
	double doChoice(Move *move, Loop **returnLoop, EnergyModel *energyModel);
	void printMove(Loop *comefrom, char *structure_p, BaseType *seq_p);
	BaseType *getLocation(Move *move, int index);
	BaseType *verifyLoop(BaseType *incoming_sequence, Loop *from);

	MultiLoop(void);
	MultiLoop(int branches, int *sidelengths, BaseType **sequences);
	~MultiLoop(void);

	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end, EnergyModel *energyModel);
	friend Loop * Loop::performDeleteMove(Move *move, EnergyModel *energyModel);
	friend void Loop::performComplexSplit(
		Move *move, Loop **firstOpen, Loop **secondOpen, EnergyModel *energyModel);
	string typeInternalsToString(void);

private:
	int *sidelen;
	BaseType **seqs;
};

class OpenLoop: public Loop {
public:
	void calculateEnergy(EnergyModel *energyModel);
	void calculateEnthalpy(EnergyModel *energyModel);
	void generateMoves(EnergyModel *energyModel);
	void generateDeleteMoves(EnergyModel*);
	void regenerateDeleteMoves(EnergyModel *energyModel);
	Move *getChoice(SimTimer& timer, Loop *from);
	double doChoice(Move *move, Loop **returnLoop, EnergyModel *energyModel);
	void printMove(Loop *comefrom, char *structure_p, BaseType *seq_p);
	BaseType *getLocation(Move *move, int index);
	BaseType *verifyLoop(BaseType *incoming_sequence,  Loop *from);

	// OpenLoop::getFreeBases returns the base composition information for the
	//   open loop. Return form is a pointer to an array of size 5, containing
	//   frequence counts of each type of base, plus invalid/out of bounds bases.
	//   Calling function needs to free the array.

	// Possibly change this so it returns the index of a static array inside the function,
	// so there's no memory overhead for the function call. (Data is accessed single-threaded,
	// and information is copied out by the calling function.) Could also make it a private data
	// member that's just returned.  Not completed currently.
	// FD: 2016-11-14 the above is now implemented.
	BaseCount& getFreeBases();
	OpenInfo& getOpenInfo(void);
	HalfContext getHalfContext(int, int);

	// check for cotranscriptional folding.
	bool nucleotideIsActive(
		const BaseType* sequence, const BaseType* initial,
		const int pos, EnergyModel *energyModel);
	bool nucleotideIsActive(
		const BaseType* sequence, const BaseType* initial,
		const int pos1, const int pos2, EnergyModel *energyModel);
	BaseType* getBase(char type, int index);
	BaseType* getBase(char type, int index, HalfContext);
	static void performComplexJoin(
		OpenLoop **oldLoops, OpenLoop **newLoops, BaseType *types, int *index,
		HalfContext[], EnergyModel *energyModel);

	OpenLoop(void);
	OpenLoop(int branches,  int *sidelengths, BaseType **sequences);
	~OpenLoop(void);
	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end, EnergyModel *energyModel);
	friend Loop * Loop::performDeleteMove(Move *move, EnergyModel *energyModel);
	friend void Loop::performComplexSplit(
		Move *move, Loop **firstOpen, Loop **secondOpen, EnergyModel *energyModel);
	string typeInternalsToString(void);
	void parseLocalContext(int);

	// non-private because we trust each other;
	// so: only the loop itself is allowed to set these.
	BaseCount exposedBases = BaseCount();
	bool updatedContext = false; // this is the update toggle for exposedBases.

	OpenInfo openInfo;
	bool initial = false; // FD: if true, then the loop is the initial open loop and seqs[0][0] is out of bounds.

private:
	int *sidelen;
	BaseType **seqs;
};

#endif
