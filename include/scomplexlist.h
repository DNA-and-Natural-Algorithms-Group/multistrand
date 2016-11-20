/*
 Copyright (c) 2007-2008 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
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

class SComplexList {
public:

	SComplexList(EnergyModel *energyModel);
	~SComplexList(void);

	SComplexListEntry *addComplex(StrandComplex *newComplex);
	void initializeList(void);
	void regenerateMoves(void);
	double getTotalFlux(void);
	double getJoinFlux(void);
	double getJoinFluxArr(void);
	double computeArrBiRate(SComplexListEntry*);
	double cycleCrossRateArr(StrandOrdering*, StrandOrdering*);
	void computeCrossRateArr(OpenLoop*, OpenLoop*);
	int getCount(void);
	double *getEnergy(int volume_flag);
	void printComplexList();
	SComplexListEntry *getFirst(void);
	int doBasicChoice(double choice, double newtime);
	void findJoinNucleotides(BaseType, int, BaseCount&, SComplexListEntry*, JoinCriterea&);
	void doJoinChoice(double choice);
	bool checkStopComplexList(class complexItem *stoplist);
	string toString(void);
	void updateOpenInfo(void);

private:
	bool checkStopComplexList_Bound(class complexItem *stoplist);
	bool checkStopComplexList_Structure_Disassoc(class complexItem *stoplist);
	bool checkLooseStructure(char *our_struc, char *stop_struc, int count);
	bool checkCountStructure(char *our_struc, char *stop_struc, int count);

	int numentries = 0;
	int idcounter = 0;

	SComplexListEntry* first = NULL;
	EnergyModel* eModel = NULL;

	double joinRate = 0.0;

};

class SComplexListEntry {
public:
	SComplexListEntry(StrandComplex *newComplex, int newid);
	~SComplexListEntry(void);
	void initializeComplex(void);
	void regenerateMoves(void);
	void fillData(EnergyModel *em);
	string toString(EnergyModel *em);
	void dumpComplexEntryToPython(int *our_id, char **names, char **sequence, char **structure, double *our_energy);
	int id;
	StrandComplex* thisComplex;
	double energy;
	energyS ee_energy;
	double rate;
	SComplexListEntry *next;
};

#endif

