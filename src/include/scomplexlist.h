/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */
#ifndef __SCOMPLEXLIST_H__
#define __SCOMPLEXLIST_H__

#include "scomplex.h"
#include "energymodel.h"
#include "optionlists.h"
#include <stdio.h>

#include <iostream>

using std::cout;

class SComplexListEntry;
class JoinCriterea;
class SimTimer;



class SComplexList {
public:

	SComplexList(EnergyModel *energyModel);
	~SComplexList(void);

	SComplexListEntry *addComplex(StrandComplex *newComplex);
	void initializeList(void);
	void regenerateMoves(void);
	double getTotalFlux(void);
	double getJoinFlux(void);

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
	JoinCriteria cycleForJoinChoice(double choice);
	JoinCriteria cycleForJoinChoiceArr(double choice);
	JoinCriteria findJoinNucleotides(BaseType, int, BaseCount&, SComplexListEntry*, HalfContext* = NULL);
	double doJoinChoice(double choice);
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

	double joinRate = 0.0;

}
;

class SComplexListEntry {
public:
	SComplexListEntry(StrandComplex *newComplex, int newid);
	~SComplexListEntry(void);
	void initializeComplex(void);
	void regenerateMoves(void);
	void fillData(EnergyModel *em);
	string toString(EnergyModel *em);
//	void dumpComplexEntryToPython(int *our_id, char **names, char **sequence, char **structure, double *our_energy, double *our_enthalpy);
	void dumpComplexEntryToPython(ExportData& data);
//	void dumpComplexEntryToPython(int *our_id, char **names, string* sequence, string* structure, double *our_energy, double *our_enthalpy);

	int id;
	StrandComplex* thisComplex;
	energyS ee_energy;
	double energy;
	double rate;

	SComplexListEntry *next;
};

#endif

