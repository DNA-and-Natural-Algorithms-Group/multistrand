/*
 Copyright (c) 2007-2008 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 */

/* Energy Model class header */
#ifndef __ENERGYMODEL_H__
#define __ENERGYMODEL_H__

#include <stdio.h>
#include <python2.7/Python.h>
#include <string>
#include <moveutil.h>
#include <sequtil.h>

using std::string;

class SimOptions;
class Loop;
class EnergyOptions;

const int VIENNA = 0;
const int MFOLD = 1;

const int pairs_vienna[5] = { 0, 4, 3, 2, 1 };
const int pairs_mfold[5] = { 0, 4, 3, 2, 1 };
extern int pairs[5];

const int pairtypes_vienna[5][5] = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 5 }, { 0, 0, 0, 1, 0 }, { 0, 0, 2, 0, 3 }, { 0, 6, 0, 4, 0 } };
const int pairtypes_mfold[5][5] = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 1 }, { 0, 0, 0, 2, 0 }, { 0, 0, 3, 0, 5 }, { 0, 4, 0, 6, 0 } };

extern int pairtypes[5][5];

const int basepair_sw_vienna[8] = { 0, 2, 1, 4, 3, 6, 5, 7 };
const int basepair_sw_mfold[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
const int basepair_sw_mfold_actual[8] = { 0, 4, 3, 2, 1, 6, 5, 7 }; // Why do this? Vienna's parameter file stores pairings in the opposite ordering. So for one of them, we need to swap basepairs to get the correct ordering, in the other one, we don't.
extern int basepair_sw[8]; // = {0,0,0,0,0,0,0,0};

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
	double arrheniusLoopEnergy(char* seq, int size);
	double saltCorrection(int size);
	void setArrheniusRate(double ratesArray[], EnergyOptions* options, double temperature, int left, int right);
	void computeArrheniusRates(double temperature);
	double applyPrefactors(double tempRate, MoveType left, MoveType right);
	MoveType getPrefactorsMulti(int, int, int[]);
	MoveType prefactorOpen(int, int, int[]);
	MoveType prefactorInternal(int, int);

	void printPrecomputedArrRates(void);
	void printkBikUni(void);

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
//	virtual void eStackEnergy(int type1, int type2, energyS *energy) = 0;

	SimOptions* simOptions;

protected:
	long dangles;
	double arrheniusRates[MOVETYPE_SIZE * MOVETYPE_SIZE];

};

class NupackEnergyModel: public EnergyModel {

public:
	NupackEnergyModel(void);
	NupackEnergyModel(PyObject* options);
	NupackEnergyModel(SimOptions* options);
	~NupackEnergyModel(void);
	void processOptions();

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

//	void eStackEnergy(int type1, int type2, energyS *energy);

private:

	// All energy units are integers, in units of .01 kcal/mol, as used by ViennaRNA

	// Stacking Info
	double stack_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASEPAIRS_NUPACK]; // Delta G's for stacks, matrix form, at 37 degrees C.
	int stack_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASEPAIRS_NUPACK]; // Delta H's for stacks, matrix form, for use comparing to dG at 37 deg C.

	// Hairpin Info
	double hairpin_37_dG[31];
	int hairpin_37_dH[31];
	double hairpin_mismatch_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];
	int hairpin_mismatch_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];

	double hairpin_triloop_37_dG[1024];
	int hairpin_triloop_37_dH[1024];
	// This needs about 1024 doubles to store the entire matrix.
	// lookups on this are then 4 shifts, 5adds, 5 dec 1s, but better than a strstr.

	double hairpin_tetraloop_37_dG[4096];
	int hairpin_tetraloop_37_dH[4096];

	// Bulge Info
	double bulge_37_dG[31];
	int bulge_37_dH[31];

	// Internal Loop Info
	double internal_37_dG[31];
	double internal_mismatch_37_dG[NUM_BASES][NUM_BASES][NUM_BASEPAIRS_NUPACK];
	int internal_37_dH[31];
	int internal_mismatch_37_dH[NUM_BASES][NUM_BASES][NUM_BASEPAIRS_NUPACK];

	double maximum_NINIO;
	int maximum_NINIO_dH;
	double ninio_correction_37[5];
	int ninio_correction_37_dH[5];

	/* special internal loop lookup tables */
	double internal_1_1_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];
	int internal_1_1_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];

	double internal_2_1_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];
	int internal_2_1_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];

	double internal_2_2_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES][NUM_BASES][NUM_BASES];
	int internal_2_2_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES][NUM_BASES][NUM_BASES];

	// Multiloop Info
	double multiloop_base;
	double multiloop_closing;
	double multiloop_internal;

	int multiloop_base_dH;
	int multiloop_closing_dH;
	int multiloop_internal_dH;

	// should this really not be tied to one of the others?
	// it actually appears to be the scaling term for internal, hairpin and (the outdated) multiloop mismatches.

	// Dangle Info for multiloops and open loops
	double dangle_3_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASES];
	double dangle_5_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASES];
	int dangle_3_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASES];
	int dangle_5_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASES];

	// Terminal AU penalty for multiloops and open loops (included in mismatch penalties elsewhere
	// Appears to be pure dH.
	double terminal_AU;
	int terminal_AU_dH;

	// Logarithmic loop penalty. Doesn't seem to change for DNA/RNA?
	double log_loop_penalty_37;
	double log_loop_penalty;

	// biomolecular penalty
	double bimolecular_penalty;
	int bimolecular_penalty_dH;

	// Kinetic rate toggle. 0 = kawasaki, 1 = metropolis, 2 = entropy/enthalpy, defaults to 2.
	long kinetic_rate_method;
	double kBoltzmann;
	double current_temp;
	double _RT;
	double joinrate;
	double dG_volume;
	double dG_assoc;

	double waterdensity;
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
	void internal_set_hairpin_tetraloop_parameters(FILE *fp, char *buffer);
	void internal_set_hairpin_triloop_parameters(FILE *fp, char *buffer);
	void internal_set_hairpin_mismatch_enthalpies(FILE *fp, char *buffer);
	void internal_set_hairpin_tetraloop_parameters_enthalpy(FILE *fp, char *buffer);
	void internal_set_hairpin_triloop_parameters_enthalpy(FILE *fp, char *buffer);
	void internal_set_bulge_energies(FILE *fp, char *buffer);
	void internal_set_bulge_enthalpies(FILE *fp, char *buffer);
	void internal_set_interior_loop_energies(FILE *fp, char *buffer);
	void internal_set_interior_loop_enthalpies(FILE *fp, char *buffer);
	void internal_set_interior_loop_mismatch_energies(FILE *fp, char *buffer);
	void internal_set_interior_loop_mismatch_enthalpies(FILE *fp, char *buffer);
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
	char *internal_read_array_data(FILE *fp, char *buffer, char* start_loc, int *read_loc, int size);
	char *internal_read_array_data(FILE *fp, char *buffer, char* start_loc, double *read_loc, int size);
	void internal_set_dangle_5_energies(FILE *fp, char *buffer);
	void internal_set_dangle_3_energies(FILE *fp, char *buffer);
	void internal_set_dangle_5_enthalpies(FILE *fp, char *buffer);
	void internal_set_dangle_3_enthalpies(FILE *fp, char *buffer);

	double setWaterDensity(double temp);
};

#endif /* __ENERGYMODEL_H__ */
