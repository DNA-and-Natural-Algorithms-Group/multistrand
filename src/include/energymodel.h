/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */

/* Energy Model class header */
#ifndef __ENERGYMODEL_H__
#define __ENERGYMODEL_H__

#include <stdio.h>
#include <Python.h>
#include <string>
#include <array>
#include <moveutil.h>
#include <sequtil.h>

#include "rapidjson/document.h"

using std::string;
using std::array;

class SimOptions;
class Loop;
class EnergyOptions;

const int VIENNA = 0;
const int MFOLD = 1;

const int pairs_vienna[5] = { 0, 4, 3, 2, 1 };
const int pairs_mfold[5] = { 0, 4, 3, 2, 1 };
extern  int pairs[5];

const int pairtypes_vienna[5][5] = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 5 }, { 0, 0, 0, 1, 0 }, { 0, 0, 2, 0, 3 }, { 0, 6, 0, 4, 0 } };
const int pairtypes_mfold[5][5] = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 1 }, { 0, 0, 0, 2, 0 }, { 0, 0, 3, 0, 5 }, { 0, 4, 0, 6, 0 } };

extern int pairtypes[5][5];

const int basepair_sw_vienna[8] = { 0, 2, 1, 4, 3, 6, 5, 7 };
const int basepair_sw_mfold[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
const int basepair_sw_mfold_actual[8] = { 0, 4, 3, 2, 1, 6, 5, 7 }; // Why do this? Vienna's parameter file stores pairings in the opposite ordering. So for one of them, we need to swap basepairs to get the correct ordering, in the other one, we don't.
extern  int basepair_sw[8]; // = {0,0,0,0,0,0,0,0};

int baseLookup(char base);

const double nupackInfinte = 100000.0;
const double gasConstant = 0.0019872041;


enum LoopType {
	openLoop, interiorLoop, bulgeLoop, stackLoop, hairpinLoop, multiLoop, LOOPTYPE_SIZE
};

class energyS {
public:
	double dH; // enthalpy
	double nTdS; // entropy - actually -TdS, such that dH + nTds = dG
};

class EnergyModel {

public:
	EnergyModel(void);
	EnergyModel(PyObject *options);

	// Implemented methods
	bool useArrhenius(void);
	double singleStrandedStacking(char* sequence, int length);
	double initializationPenalty(int, int, int);
	double arrheniusLoopEnergy(char* seq, int size);
	double saltCorrection(void);
	void setArrheniusRate(double ratesArray[], EnergyOptions* options, double temperature, int left, int right);
	void computeArrheniusRates(double temperature);
	double applyPrefactors(double tempRate, MoveType left, MoveType right);
	MoveType getPrefactorsMulti(int, int, int[]);
	MoveType prefactorOpen(int, int, int[]);
	MoveType prefactorInternal(int, int);

	void writeConstantsToFile(void);
	double fastestUniRate(void);

	// Virtual methods

	virtual ~EnergyModel(void);

	virtual double returnRate(double start_energy, double end_energy, int enth_entr_toggle) = 0;
	virtual double getJoinRate_NoVolumeTerm(void) = 0;
	virtual double getJoinRate(void) = 0;
	virtual double getVolumeEnergy(void) =0;
	virtual double getAssocEnergy(void) =0;

	virtual double StackEnergy(int i, int j, int p, int q) = 0;
	//     This is: 5' i p 3'
	//              3' j q 5'

	virtual double BulgeEnergy(int i, int j, int p, int q, int bulgesize) = 0;
	//    This is: 5' i x p 3'
	//             3' j x q 5'
	//  Where one of the x's is a sequence of length bulgesize and the other is size 0.

	virtual double InteriorEnergy(char *seq1, char *seq2, int size1, int size2) = 0;
	// This is: 5' seq1 3'
	//          3' seq2 5'

	virtual double HairpinEnergy(char *seq, int size) = 0;
	// just passing in the whole sequence

	virtual double MultiloopEnergy(int size, int *sidelen, char **sequences) = 0;
	virtual double OpenloopEnergy(int size, int *sidelen, char **sequences) = 0;

	// FD January 2018: Adding the mimicking enthalpy functions
	virtual double StackEnthalpy(int i, int j, int p, int q) = 0;
	virtual double BulgeEnthalpy(int i, int j, int p, int q, int bulgesize) = 0;
	virtual double InteriorEnthalpy(char *seq1, char *seq2, int size1, int size2) = 0;
	virtual double HairpinEnthalpy(char *seq, int size) = 0;
	virtual double MultiloopEnthalpy(int size, int *sidelen, char **sequences) = 0;
	virtual double OpenloopEnthalpy(int size, int *sidelen, char **sequences) = 0;

	SimOptions* simOptions;

	// for co-transcriptional, save the curent amount of activated nucleotides
	int numActiveNT = NULL;
	// if inspection is used, all transitions occur at rate 1.0 /s
	bool inspection = false;

	// keep track of the used paths to parameter files
	string paramFiles[2] = {"unset",  "unset"};



protected:
	long dangles;
	double arrheniusRates[MOVETYPE_SIZE * MOVETYPE_SIZE];

};

struct hairpin_energies{

	array<double, 31> basic;
	array<array<array<double, BASES>, BASES>, PAIRS_NUPACK> mismatch;
	array<double, 1024> triloop;
	array<double, 4096> tetraloop;

	hairpin_energies(){
		basic[0] = 0.0;
	}

};


struct internal_energies{

	// Internal Loop Info
	double basic[31];
	double mismatch[BASES][BASES][PAIRS_NUPACK];

	double maximum_NINIO;
	double ninio_correction[5];

	/* special internal loop lookup tables */
	double internal_1_1[PAIRS_NUPACK][PAIRS_NUPACK][BASES][BASES];
	double internal_2_1[PAIRS_NUPACK][BASES][PAIRS_NUPACK][BASES][BASES];
	double internal_2_2[PAIRS_NUPACK][PAIRS_NUPACK][BASES][BASES][BASES][BASES];

};


struct multiloop_energies{

	double base;
	double closing;
	double internal;

	// Dangle Info for multiloops and open loops
	double dangle_3[PAIRS_NUPACK][BASES];
	double dangle_5[PAIRS_NUPACK][BASES];

};


class NupackEnergyModel: public EnergyModel {

public:
	NupackEnergyModel(void);
	NupackEnergyModel(PyObject* options);
	NupackEnergyModel(SimOptions* options);
//	~NupackEnergyModel(void);

	double returnRate(double start_energy, double end_energy, int enth_entr_toggle);
	double returnRate(energyS &start_energy, energyS &end_energy);

	double getJoinRate(void);
	double getJoinRate_NoVolumeTerm(void);

	double getVolumeEnergy(void);
	double getAssocEnergy(void);

	double StackEnergy(int i, int j, int p, int q);
	double BulgeEnergy(int i, int j, int p, int q, int bulgesize);
	double InteriorEnergy(char *seq1, char *seq2, int size1, int size2);
	double HairpinEnergy(char *seq, int size);
	double MultiloopEnergy(int size, int *sidelen, char **sequences);
	double OpenloopEnergy(int size, int *sidelen, char **sequences);

	// FD jan 2018: adding the corresponding enthalpy functions
	double StackEnthalpy(int i, int j, int p, int q);
	double BulgeEnthalpy(int i, int j, int p, int q, int bulgesize);
	double InteriorEnthalpy(char *seq1, char *seq2, int size1, int size2);
	double HairpinEnthalpy(char *seq, int size);
	double MultiloopEnthalpy(int size, int *sidelen, char **sequences);
	double OpenloopEnthalpy(int size, int *sidelen, char **sequences);

private:

	void processOptions();
	FILE* openFiles(char*, string&, string&, int);

	// FD jan 2018: helper functions, now seperated out
	double HairpinEnergy(char *seq, int size, hairpin_energies&);
	double InteriorEnergy(char *seq1, char *seq2, int size1, int size2, internal_energies& internal);
	double BulgeEnergy(int i, int j, int p, int q, int bulgesize, array<double,31>, array<array<double, PAIRS_NUPACK>, PAIRS_NUPACK> );
	double MultiloopEnergy(int size, int *sidelen, char **sequences, multiloop_energies& multiloop);
	double OpenloopEnergy(int size, int *sidelen, char **sequences, multiloop_energies& multiloop);

	// JS: All energy units are integers, in units of .01 kcal/mol, as used by ViennaRNA
	// FD: In 2.0, the units changed to 0.01 kcal/mol for dH and kcal/mol for dG.
	// FD: as of jan 2018, dH and dG are now both kcal/mol.

	// Stacking Info
	array<array<double, PAIRS_NUPACK>, PAIRS_NUPACK> stack_37_dG; // Delta G's for stacks, matrix form, at 37 degrees C.
	array<array<double, PAIRS_NUPACK>, PAIRS_NUPACK> stack_37_dH; // Delta H's for stacks, matrix form, for use comparing to dG at 37 deg C.

	// Hairpin Info
	hairpin_energies hairpin_dG, hairpin_dH;

	// Bulge Info
	array<double,31> bulge_37_dG, bulge_37_dH;

	// Internal Loop Info
	internal_energies internal_dG, internal_dH;

	// Multiloop Info
	multiloop_energies multiloop_dG, multiloop_dH;

	// Terminal AU penalty for multiloops and open loops (included in mismatch penalties elsewhere
	// Appears to be pure dH.
	double terminal_AU;
	double terminal_AU_dH;

	// Logarithmic loop penalty. Doesn't seem to change for DNA/RNA?
	double log_loop_penalty_37;
	double log_loop_penalty;

	// biomolecular penalty
	double bimolecular_penalty;
	double bimolecular_penalty_dH;

	// Kinetic rate toggle. 0 = kawasaki, 1 = metropolis, 2 = entropy/enthalpy, defaults to 2.
	long kinetic_rate_method;
	double kBoltzmann;
	double current_temp;
	double _RT;
	double joinrate;
	double dG_volume;
	double dG_assoc;

	double biscale;
	double uniscale;

	bool gtenable;
	int internal;
	long logml;

	// data loading functions:
	void setupRates();

	void internal_set_stack_energies(FILE *fp, char *buffer);
	void internal_set_stack_enthalpies(FILE *fp, char *buffer);
	void internal_set_hairpin_energies(FILE *fp, char *buffer);
	void internal_set_hairpin_enthalpies(FILE *fp, char *buffer);
	void internal_set_hairpin_mismatch_energies(FILE *fp, char *buffer);
	void internal_set_hairpin_tetraloop_parameters(FILE *fp, char *buffer, hairpin_energies& );
	void internal_set_hairpin_triloop_parameters(FILE *fp, char *buffer, hairpin_energies& );
	void internal_set_hairpin_mismatch_enthalpies(FILE *fp, char *buffer);
	void internal_set_bulge_energies(FILE *fp, char *buffer);
	void internal_set_bulge_enthalpies(FILE *fp, char *buffer);
	void internal_set_interior_loop_energies(FILE *fp, char *buffer);
	void internal_set_interior_loop_enthalpies(FILE *fp, char *buffer);
	void internal_set_interior_loop_mismatch_energies(FILE *fp, char *buffer, internal_energies& );
	void internal_set_interior_1_1_energies(FILE *fp, char *buffer);
	void internal_set_interior_1_1_enthalpies(FILE *fp, char *buffer);
	void internal_set_interior_2_1_energies(FILE *fp, char *buffer);
	void internal_set_interior_2_1_enthalpies(FILE *fp, char *buffer);
	void internal_set_interior_2_2_energies(FILE *fp, char *buffer);
	void internal_set_interior_2_2_enthalpies(FILE *fp, char *buffer);
	void internal_set_multiloop_parameters(FILE *fp, char *buffer);
	void internal_set_multiloop_parameters_enthalpies(FILE *fp, char *buffer);
	void internal_set_at_penalty(FILE *fp, char *buffer);
	void internal_set_at_penalty_enthalpy(FILE *fp, char *buffer);
	void internal_set_ninio_parameters(FILE *fp, char *buffer);
	void internal_set_ninio_parameters_enthalpy(FILE *fp, char *buffer);
	void internal_set_bimolecular_penalty(FILE *fp, char *buffer);
	void internal_set_bimolecular_penalty_dH(FILE *fp, char *buffer);
	char *internal_read_array_data(FILE *fp, char *buffer, char* start_loc, double *read_loc, int size);
	void internal_set_dangle_5_energies(FILE *fp, char *buffer);
	void internal_set_dangle_3_energies(FILE *fp, char *buffer);
	void internal_set_dangle_5_enthalpies(FILE *fp, char *buffer);
	void internal_set_dangle_3_enthalpies(FILE *fp, char *buffer);

	void internal_set_stack(rapidjson::Document &d);
	void internal_set_hairpin(rapidjson::Document &d);
	void internal_set_bulge(rapidjson::Document &d);
	void internal_set_interior_loop(rapidjson::Document &d);
	void internal_set_interior_1_1(rapidjson::Document &d);
	void internal_set_interior_2_1(rapidjson::Document &d);
	void internal_set_interior_2_2(rapidjson::Document &d);
	void internal_set_dangle_5(rapidjson::Document &d);
	void internal_set_dangle_3(rapidjson::Document &d);
	void internal_set_multiloop_parameters(rapidjson::Document &d);
	void internal_set_at_penalty(rapidjson::Document &d);
	void internal_set_bimolecular_penalty(rapidjson::Document &d);
	void internal_set_ninio_parameters(rapidjson::Document &d);
	void internal_set_hairpin_tetraloop_parameters(rapidjson::Document &d);
	void internal_set_hairpin_triloop_parameters(rapidjson::Document &d);
	void internal_set_hairpin_mismatch(rapidjson::Document &d);
	void internal_set_interior_loop_mismatch(rapidjson::Document &d);








	double setWaterDensity(double temp);
};

#endif /* __ENERGYMODEL_H__ */
