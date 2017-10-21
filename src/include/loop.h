/*
Copyright (c) 2017 California Institute of Technology. All rights reserved.
Multistrand nucleic acid kinetic simulator
help@multistrand.org
*/
/* Loop class header for main loop class and derivatives. */

#ifndef __LOOP_H__
#define __LOOP_H__

#include <string>
#include <vector>
#include "energymodel.h"
#include "move.h"
#include "moveutil.h"

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
	inline double getEnergy(void);
	inline double getTotalRate(void);
	char getType(void);
	Loop(void);
	virtual ~Loop(void);
	virtual void calculateEnergy(void) = 0;
	virtual void generateMoves(void) = 0;
	virtual void generateDeleteMoves(void) = 0;
	virtual Move *getChoice(double *randomchoice, Loop *from) = 0;
	virtual double doChoice(Move *move, Loop **returnLoop) = 0;
	virtual char *getLocation(Move *move, int index) =0;
	virtual char *verifyLoop(char *incoming_sequence, Loop *from) =0;
	virtual string typeInternalsToString(void) = 0;
	virtual void printMove(Loop *comefrom, char *structure_p, char *seq_p) = 0;
	Loop *getAdjacent(int index);
	int getCurAdjacent(void);
	int getNumAdjacent(void);
	void addAdjacent(Loop *loopToAdd);
	void initAdjacency(int index);
	int replaceAdjacent(Loop *loopToReplace, Loop *loopToReplaceWith);
	void cleanupAdjacent(void); // sets adjacentLoops up to be deleted.
	double returnEnergies(Loop *comefrom); // returns the total energy of all loops underneath this one.
	double returnFlux(Loop *comefrom); // returns the total rate of all loops underneath this one.
	void firstGen(Loop *comefrom);
	static void SetEnergyModel(EnergyModel *newEnergyModel);
	static EnergyModel *GetEnergyModel(void);
	static RateArr generateDeleteMoveRate(Loop *start, Loop *end);
	static Loop *performDeleteMove(Move *move);
	static void performComplexSplit(Move *move, Loop **firstOpen, Loop **secondOpen);

	// FD functions
    static bool identify(Loop*, Loop* , char , char );
    static std::tuple<int,int> findExternalAdjacent(Loop*, Loop*);
    static std::pair<Loop*, Loop*> orderMyLoops(Loop*, Loop*, char);

	string toString(void);
	string toStringShort(void);
	void printAllMoves(Loop*);
	void generateAndSaveDeleteMove(Loop*, int);

	// FD: moving private to public
	int numAdjacent;

protected:
	static EnergyModel *energyModel;

	Loop** adjacentLoops;
	int curAdjacent;
	double energy;
	bool energyComputed = false;
	double totalRate;
	MoveContainer *moves;
	char identity;
	int add_index;
};

class StackLoop: public Loop {
public:
	void calculateEnergy(void);
	void generateMoves(void);
	void generateDeleteMoves(void);
	Move *getChoice(double *randnum, Loop *from);
	double doChoice(Move *move, Loop **returnLoop);
	void printMove(Loop *comefrom, char *structure_p, char *seq_p);
	char *getLocation(Move *move, int index);
	char *verifyLoop(char *incoming_sequence,  Loop *from);
	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end);
	friend Loop * Loop::performDeleteMove(Move *move);
	friend void Loop::performComplexSplit(Move *move, Loop **firstOpen, Loop **secondOpen);
	StackLoop(void);
	StackLoop(char *seq1, char *seq2, Loop *left = NULL, Loop *right = NULL);
	string typeInternalsToString(void);

private:
	char *seqs[2];
};

class HairpinLoop: public Loop {
public:
	void calculateEnergy(void);
	void generateMoves(void);
	void generateDeleteMoves(void);
	Move *getChoice(double *randnum, Loop *from);
	double doChoice(Move *move, Loop **returnLoop);
	void printMove(Loop *comefrom, char *structure_p, char *seq_p);
	char *getLocation(Move *move, int index);
	char *verifyLoop(char *incoming_sequence,  Loop *from);

	HairpinLoop(void);
	HairpinLoop( int size, char *hairpin_sequence, Loop *previous = NULL);
	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end);
	friend Loop * Loop::performDeleteMove(Move *move);
	friend void Loop::performComplexSplit(Move *move, Loop **firstOpen, Loop **secondOpen);
	string typeInternalsToString(void);

private:
	int hairpinsize;
	char *hairpin_seq;
};

class BulgeLoop: public Loop {
public:
	void calculateEnergy(void);
	void generateMoves(void);
	void generateDeleteMoves(void);
	Move *getChoice(double *randnum, Loop *from);
	double doChoice(Move *move, Loop **returnLoop);
	void printMove(Loop *comefrom, char *structure_p, char *seq_p);
	char *getLocation(Move *move, int index);
	char *verifyLoop(char *incoming_sequence, Loop *from);
	BulgeLoop(void);
	BulgeLoop(int size1, int size2, char *bulge_sequence1, char *bulge_sequence2, Loop *left = NULL, Loop *right = NULL);
	friend Loop * Loop::performDeleteMove(Move *move);
	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end);
	friend void Loop::performComplexSplit(Move *move, Loop **firstOpen, Loop **secondOpen);
	string typeInternalsToString(void);

private:
	int bulgesize[2];
	char *bulge_seq[2];
};

class InteriorLoop: public Loop {
public:
	void calculateEnergy(void);
	void generateMoves(void);
	void generateDeleteMoves(void);
	Move *getChoice(double *randnum, Loop *from);
	double doChoice(Move *move, Loop **returnLoop);
	void printMove(Loop *comefrom, char *structure_p, char *seq_p);
	char *getLocation(Move *move, int index);
	char *verifyLoop(char *incoming_sequence,  Loop *from);
	InteriorLoop(void);
	InteriorLoop(int size1, int size2, char *int_seq1, char *int_seq2, Loop *left = NULL, Loop *right = NULL);

	friend Loop * Loop::performDeleteMove(Move *move);
	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end);
	friend void Loop::performComplexSplit(Move *move, Loop **firstOpen, Loop **secondOpen);
	string typeInternalsToString(void);

private:
	int sizes[2];
	char *int_seq[2];
};

class MultiLoop: public Loop {
public:
	void calculateEnergy(void);
	void generateMoves(void);
	void generateDeleteMoves(void);
	Move *getChoice(double *randnum, Loop *from);
	double doChoice(Move *move, Loop **returnLoop);
	void printMove(Loop *comefrom, char *structure_p, char *seq_p);
	char *getLocation(Move *move, int index);
	char *verifyLoop(char *incoming_sequence, Loop *from);
	MultiLoop(void);
	MultiLoop(int branches, int *sidelengths, char **sequences);
	~MultiLoop(void);

	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end);
	friend Loop * Loop::performDeleteMove(Move *move);
	friend void Loop::performComplexSplit(Move *move, Loop **firstOpen, Loop **secondOpen);
	string typeInternalsToString(void);

private:
	int *sidelen;
	char **seqs;
};

class OpenLoop: public Loop {
public:
	void calculateEnergy(void);
	void generateMoves(void);
	void generateDeleteMoves(void);
	Move *getChoice(double *randomchoice, Loop *from);
	double doChoice(Move *move, Loop **returnLoop);
	void printMove(Loop *comefrom, char *structure_p, char *seq_p);
	char *getLocation(Move *move, int index);
	char *verifyLoop(char *incoming_sequence,  Loop *from);

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
	bool nucleotideIsActive(const char* sequence, const char* initial, const int pos);
	bool nucleotideIsActive(const char* sequence, const char* initial, const int pos1, const int pos2);

	char* getBase(char type, int index);
	char* getBase(char type, int index, HalfContext);

	OpenLoop(void);
	OpenLoop(int branches,  int *sidelengths, char **sequences);
	~OpenLoop(void);
	friend RateArr Loop::generateDeleteMoveRate(Loop *start, Loop *end);
	friend Loop * Loop::performDeleteMove(Move *move);
	friend void Loop::performComplexSplit(Move *move, Loop **firstOpen, Loop **secondOpen);
	static void performComplexJoin(OpenLoop **oldLoops, OpenLoop **newLoops, char *types, int *index, HalfContext[], bool);
	string typeInternalsToString(void);
	void parseLocalContext(int);

	// non-private because we trust each other;
	// so: only the loop itself is allowed to set these.
	BaseCount exposedBases = BaseCount();
	bool updatedContext2 = false; // this is the update toggle for exposedBases.

	OpenInfo openInfo;
	bool initial = false; // FD: if true, then the loop is the initial open loop and seqs[0][0] is out of bounds.
private:

	int *sidelen;
	char **seqs;



};

#endif
