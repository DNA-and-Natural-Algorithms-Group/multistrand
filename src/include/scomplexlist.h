/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef __SCOMPLEXLIST_H__
#define __SCOMPLEXLIST_H__

#include <stdio.h>
#include <iostream>

#include "scomplex.h"
#include "energymodel.h"
#include "optionlists.h"
#include "simtimer.h"
#include "strandordering.h"

using std::cout;

class SComplexListEntry;
class JoinCriterea;

class SComplexList {
public:

	SComplexList(EnergyModel *energyModel);
	~SComplexList(void);

	SComplexListEntry *addComplex(StrandComplex *newComplex);
	void initializeList(void);
	void regenerateMoves(void);
	double getTotalFlux(void);
	double getJoinFlux(void);
	int getMoveCount(void);

	BaseCount getExposedBases();
	OpenInfo getOpenInfo();

	double getJoinFluxArr(void);
	double computeArrBiRate(SComplexListEntry*);
	double cycleCrossRateArr(StrandOrdering*, StrandOrdering*);
	int getCount(void);
	double *getEnergy(int volume_flag);
	void printComplexList();
	SComplexListEntry *getFirst(void);
	double doBasicChoice(SimTimer& timer);
	JoinCriteria cycleForJoinChoice(SimTimer& timer);
	JoinCriteria cycleForJoinChoiceArr(SimTimer& timer);
	JoinCriteria findJoinNucleotides(BaseType, int, BaseCount&, SComplexListEntry*, HalfContext* = NULL);
	double doJoinChoice(SimTimer& choice);
	void doJoinChoiceArr(double choice);
	bool checkStopComplexList(class complexItem *stoplist);
	string toString(void);
	void updateOpenInfo(void);

private:
	bool checkStopComplexList_Bound(class complexItem *stoplist);
	bool checkStopComplexList_Structure_Disassoc(class complexItem *stoplist);
	bool checkLooseStructure(const char *our_struc, const char *stop_struc, int count);
	bool checkCountStructure(const char *our_struc, const char *stop_struc, int count);

	int numOfComplexes = 0;
	int idcounter = 0;

	SComplexListEntry* first = NULL;
	EnergyModel* eModel = NULL;

	double joinRate = 0.0;	// joinrate is the sum of collision rates in the state.

}
;

class SComplexListEntry {
public:
	SComplexListEntry(StrandComplex *newComplex, int newid);
	~SComplexListEntry(void);
	void initializeComplex(bool debug);
	void regenerateMoves(void);
	void fillData(EnergyModel *em);
	string toString(EnergyModel *em);
	void dumpComplexEntryToPython(ExportData& data);

	int id;
	StrandComplex* thisComplex;
	energyS ee_energy;
	double energy;
	double rate;

	SComplexListEntry *next;
};

#endif

