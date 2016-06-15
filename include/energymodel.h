/*
 Copyright (c) 2007-2008 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 */

/* Energy Model class header */
#ifndef __ENERGYMODEL_H__
#define __ENERGYMODEL_H__

#include <stdio.h>
#include <python2.7/Python.h>

class SimOptions;

#define NUM_BASEPAIRS_VIENNA 8
// Vienna: 0 is invalid, then CG, GC, GU, UG, AU, UA, and Special are 1-7
#define NUM_BASEPAIRS_NUPACK 6
// MFold/Nupack:  0 is AT, then CG, GC, TA, GT, TG
#define NUM_BASES 5
// 0 is invalid, then A, C, G, U

#define BASE_A 1
#define BASE_C 2
#define BASE_G 3
#define BASE_T 4
#define BASE_U 4

#ifndef VIENNA
#define VIENNA 0
//#define ENERGYMODEL_VIENNA 0
#endif
#ifndef MFOLD
#define MFOLD 1
//#define ENERGYMODEL_NUPACK  1
#endif

const int pairs_vienna[5] = { 0, 4, 3, 2, 1 };
const int pairs_mfold[5] = { 0, 4, 3, 2, 1 };
extern int pairs[5]; // = {0,0,0,0,0} ; // move to inside class?

const int pairtypes_vienna[5][5] = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 5 }, { 0,
		0, 0, 1, 0 }, { 0, 0, 2, 0, 3 }, { 0, 6, 0, 4, 0 } };

const int pairtypes_mfold[5][5] = { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 1 }, { 0,
		0, 0, 2, 0 }, { 0, 0, 3, 0, 5 }, { 0, 4, 0, 6, 0 } };

extern int pairtypes[5][5];/* = {
 {0,0,0,0,0},
 {0,0,0,0,0},
 {0,0,0,0,0},
 {0,0,0,0,0},
 {0,0,0,0,0}
 };*/

const int basepair_sw_vienna[8] = { 0, 2, 1, 4, 3, 6, 5, 7 };
const int basepair_sw_mfold[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
const int basepair_sw_mfold_actual[8] = { 0, 4, 3, 2, 1, 6, 5, 7 }; // Why do this? Yeah, because vienna's parameter file sucks and stores pairings in the opposite ordering. So for one of them, we need to swap basepairs to get the correct ordering, in the other one, we don't.
extern int basepair_sw[8]; // = {0,0,0,0,0,0,0,0};

int baseLookup(char base);

#define INF 100000

class energyS {
public:
	double dH; // enthalpy
	double nTdS; // entropy - actually -TdS, such that dH + nTds = dG
};

class EnergyModel {
public:
	EnergyModel(void);
	EnergyModel(PyObject *options);
	virtual ~EnergyModel(void);

	virtual double returnRate(double start_energy, double end_energy,
			int enth_entr_toggle) = 0;
	//  virtual double returnRate( energyS &start_energy, energyS &end_energy );
	virtual double getJoinRate(void) = 0;
	virtual double getJoinRate_NoVolumeTerm(void) = 0;
	//  virtual double getJoinEnergy( void ) = 0;
	virtual double getVolumeEnergy(void) =0;
	virtual double getAssocEnergy(void) =0;
	virtual double StackEnergy(int i, int j, int p, int q) = 0;
	//     This is: 5' i p 3'
	//              3' j q 5'

	virtual double BulgeEnergy(int i, int j, int p, int q, int bulgesize) = 0;
	//    This is: 5' i x p 3'
	//             3' j x q 5'
	//  Where one of the x's is a sequence of length bulgesize and the other is size 0.

	virtual double InteriorEnergy(char *seq1, char *seq2, int size1,
			int size2) = 0;
	// This is: 5' seq1 3'
	//          3' seq2 5'

	virtual double HairpinEnergy(char *seq, int size) = 0;
	// just passing in the whole sequence

	virtual double MultiloopEnergy(int size, int *sidelen,
			char **sequences) = 0;
	virtual double OpenloopEnergy(int size, int *sidelen, char **sequences) = 0;

	virtual void eStackEnergy(int type1, int type2, energyS *energy) = 0;
	//  virtual void eBulgeEnergy( int type1, int type2, int bulgesize, energyS *energy );
	// virtual void eInteriorEnergy( int type1, int type2, int size1, int size2, char mismatch[4], energyS *energy );
	//virtual void eHairpinEnergy( int type1, int size, char *special , energyS *energy);
	//virtual void eMultiloopEnergy( int size, int *pairtypes, int *sidelen, char **sequences, energyS *energy   );
	//virtual void eOpenloopEnergy( int size, int *pairtypes, int *sidelen, char **sequences, energyS *energy   );

protected:
	long dangles;

};

class ViennaEnergyModel: public EnergyModel {
public:
	ViennaEnergyModel(void);
	ViennaEnergyModel(PyObject *options);
	~ViennaEnergyModel(void);

	double returnRate(double start_energy, double end_energy,
			int enth_entr_toggle);
	double returnRate(energyS &start_energy, energyS &end_energy);
	double getJoinRate(void);
	double getJoinRate_NoVolumeTerm(void);

	double getVolumeEnergy(void);
	double getAssocEnergy(void);

	//double getJoinEnergy( void );

	double StackEnergy(int i, int j, int p, int q);
	double BulgeEnergy(int i, int j, int p, int q, int bulgesize);
	double InteriorEnergy(char *seq1, char *seq2, int size1, int size2);
	double HairpinEnergy(char *seq, int size);
	double MultiloopEnergy(int size, int *sidelen, char **sequences);
	double OpenloopEnergy(int size, int *sidelen, char **sequences);

	void eStackEnergy(int type1, int type2, energyS *energy);
	void eBulgeEnergy(int type1, int type2, int bulgesize, energyS *energy);
	void eInteriorEnergy(int type1, int type2, int size1, int size2,
			char mismatch[4], energyS *energy);
	void eHairpinEnergy(int type1, int size, char *special, energyS *energy);
	void eMultiloopEnergy(int size, int *pairtypes, int *sidelen,
			char **sequences, energyS *energy);
	void eOpenloopEnergy(int size, int *pairtypes, int *sidelen,
			char **sequences, energyS *energy);

private:

	// All energy units are integers, in units of .01 kcal/mol, as used by ViennaRNA

	// Stacking Info
	int stack_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASEPAIRS_VIENNA]; // Delta G's for stacks, matrix form, at 37 degrees C.
	int stack_37_dH[NUM_BASEPAIRS_VIENNA][NUM_BASEPAIRS_VIENNA]; // Delta H's for stacks, matrix form, for use comparing to dG at 37 deg C.

	// Hairpin Info, note no enthalpies
	int hairpin_37_dG[31];
	int hairpin_mismatch_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES];

	int hairpin_triloop_37_dG[1024];
	int hairpin_triloop_37_dH;
	// This needs about 1024 ints to store the entire matrix.
	// lookups on this are then 4 shifts, 5adds, 5 dec 1s, but better than a strstr.

	int hairpin_tetraloop_37_dG[4096];
	int hairpin_tetraloop_37_dH;

	// NOTES: the DNA dataset (dna.par) has 110 tetraloops, but their code seems to only load the first 80 of them. We load all of them, so there may be slight energy differences in exactly those 30.
	// Also, their dataset only uses a single dH, rather than one for each tri/tetraloop.

	// Bulge Info
	int bulge_37_dG[31];

	// Internal Loop Info
	int internal_37_dG[31];
	int internal_mismatch_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES];

	int maximum_NINIO;
	int ninio_correction_37[5]; // Why? Must look into this.
								// It appears only [2] is used. Perhaps change
								// this to a single parameter.

	/* special internal loop lookup tables */
	int internal_1_1_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES];
	int internal_1_1_37_dH[NUM_BASEPAIRS_VIENNA][NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES];

	int internal_2_1_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES][NUM_BASES];
	int internal_2_1_37_dH[NUM_BASEPAIRS_VIENNA][NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES][NUM_BASES];

	int internal_2_2_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES][NUM_BASES][NUM_BASES];
	int internal_2_2_37_dH[NUM_BASEPAIRS_VIENNA][NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES][NUM_BASES][NUM_BASES];

	// Multiloop Info
	int multiloop_mismatch_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES]; // this data is in the files, but Vienna doesn't use it. Loading it for posterity.

	int multiloop_base;
	int multiloop_closing;
	int multiloop_internal;

	// should this really not be tied to one of the others?
	// it actually appears to be the scaling term for internal, hairpin and (the outdated) multiloop mismatches.

	int mismatch_37_dH[NUM_BASEPAIRS_VIENNA][NUM_BASES][NUM_BASES];

	// Dangle Info for multiloops and open loops
	int dangle_3_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASES];
	int dangle_5_37_dG[NUM_BASEPAIRS_VIENNA][NUM_BASES];
	int dangle_3_37_dH[NUM_BASEPAIRS_VIENNA][NUM_BASES];
	int dangle_5_37_dH[NUM_BASEPAIRS_VIENNA][NUM_BASES];

	// Terminal AU penalty for multiloops and open loops (included in mismatch penalties elsewhere
	// Appears to be pure dH.
	int terminal_AU;

	// Logarithmic loop penalty. Doesn't seem to change for DNA/RNA?
	double log_loop_penalty_37;

	// parameter type, used for the internal data loading.
	int ptype;

	// biomolecular penalty
	int bimolecular_penalty;

	// Kinetic rate toggle. 0 = kawasaki, 1 = metropolis, 2 = entropy/enthalpy, defaults to 2.
	long kinetic_rate_method;
	double _RT;
	double joinrate;

	// internal data loading functions:
	void internal_set_stack_energies(FILE *fp, char *buffer);
	void internal_set_stack_enthalpies(FILE *fp, char *buffer);
	void internal_set_hairpin_energies(FILE *fp, char *buffer);
	void internal_set_hairpin_mismatch_energies(FILE *fp, char *buffer);
	void internal_set_hairpin_tetraloop_parameters(FILE *fp, char *buffer);
	void internal_set_hairpin_triloop_parameters(FILE *fp, char *buffer);
	void internal_set_bulge_energies(FILE *fp, char *buffer);
	void internal_set_interior_loop_energies(FILE *fp, char *buffer);
	void internal_set_interior_loop_mismatch_energies(FILE *fp, char *buffer);
	void internal_set_interior_1_1_energies(FILE *fp, char *buffer);
	void internal_set_interior_1_1_enthalpies(FILE *fp, char *buffer);
	void internal_set_interior_2_1_energies(FILE *fp, char *buffer);
	void internal_set_interior_2_1_enthalpies(FILE *fp, char *buffer);
	void internal_set_interior_2_2_energies(FILE *fp, char *buffer);
	void internal_set_interior_2_2_enthalpies(FILE *fp, char *buffer);
	void internal_set_multiloop_parameters(FILE *fp, char *buffer);
	void internal_set_at_penalty(FILE *fp, char *buffer);
	void internal_set_ninio_parameters(FILE *fp, char *buffer);
	void internal_set_multiloop_mismatch_energies(FILE *fp, char *buffer);
	void internal_set_bimolecular_penalty(FILE *fp, char *buffer);
	void internal_set_mismatch_enthalpies(FILE *fp, char *buffer);
	char *internal_read_array_data(FILE *fp, char *buffer, char* start_loc,
			int *read_loc, int size);
	void internal_set_dangle_5_energies(FILE *fp, char *buffer);
	void internal_set_dangle_3_energies(FILE *fp, char *buffer);
	void internal_set_dangle_5_enthalpies(FILE *fp, char *buffer);
	void internal_set_dangle_3_enthalpies(FILE *fp, char *buffer);
};

class NupackEnergyModel: public EnergyModel {
public:
	NupackEnergyModel(void);
	NupackEnergyModel(PyObject* options);
	NupackEnergyModel(SimOptions* options);
	~NupackEnergyModel(void);
	void processOptions();

	double returnRate(double start_energy, double end_energy,
			int enth_entr_toggle);
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

	void eStackEnergy(int type1, int type2, energyS *energy);
	//void eBulgeEnergy( int type1, int type2, int bulgesize, energyS *energy );
	//void eInteriorEnergy( int type1, int type2, int size1, int size2, char mismatch[4], energyS *energy );
	//void eHairpinEnergy( int type1, int size, char *special , energyS *energy);
	//void eMultiloopEnergy( int size, int *pairtypes, int *sidelen, char **sequences, energyS *energy   );
	//  void eOpenloopEnergy( int size, int *pairtypes, int *sidelen, char **sequences, energyS *energy   );

private:
	SimOptions* simOptions;

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
	//[NUM_BASES-1][NUM_BASES-1][NUM_BASES-1][NUM_BASES-1][NUM_BASES-1];
	int hairpin_triloop_37_dH[1024];
	// This needs about 1024 doubles to store the entire matrix.
	// lookups on this are then 4 shifts, 5adds, 5 dec 1s, but better than a strstr.

	double hairpin_tetraloop_37_dG[4096];
	int hairpin_tetraloop_37_dH[4096];

	/* char hairpin_tetraloops[(7*120) + 1]; // 120 tetraloops */
	/* double hairpin_tetraloop_37_dG[120]; // why 80 here? Check Vienna Code/Datasets. */
	/* // Ok, changed them to 120 each, as the DNA dataset has at least 110 tetraloops listed... very odd that it was tossing out the last 30.  */
	/* //  int hairpin_tetraloop_37_dH;  */
	/* int hairpin_tetraloop_37_dH[120];  */

	/* char hairpin_triloops[(6*40) + 1]; */
	/* double hairpin_triloop_37_dG[40]; */
	/* int hairpin_triloop_37_dH[40]; */

	// Bulge Info
	double bulge_37_dG[31];
	int bulge_37_dH[31];

	// Internal Loop Info
	double internal_37_dG[31];
	//  int internal_mismatch_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];
	double internal_mismatch_37_dG[NUM_BASES][NUM_BASES][NUM_BASEPAIRS_NUPACK];
	int internal_37_dH[31];
	//  int internal_mismatch_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];
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
	//  int multiloop_mismatch_37_dG[NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES]; // this data is in the files, but Vienna doesn't use it. Loading it for posterity.

	double multiloop_base;
	double multiloop_closing;
	double multiloop_internal;

	int multiloop_base_dH;
	int multiloop_closing_dH;
	int multiloop_internal_dH;

	// should this really not be tied to one of the others?
	// it actually appears to be the scaling term for internal, hairpin and (the outdated) multiloop mismatches.

	//int mismatch_37_dH[NUM_BASEPAIRS_NUPACK][NUM_BASES][NUM_BASES];

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

	// parameter type, used for the internal data loading.
	//  int ptype;
	// JS: no longer needed as of rebuilt options object.

	// biomolecular penalty
	double bimolecular_penalty;
	int bimolecular_penalty_dH;

	// Kinetic rate toggle. 0 = kawasaki, 1 = metropolis, 2 = entropy/enthalpy, defaults to 2.
	long kinetic_rate_method;
	double kBoltzmann;
	double current_temp;
	double _RT;
	double joinrate;
	//  double joinenergy;
	double dG_volume;
	double dG_assoc;

	double waterdensity;
	double biscale;
	double uniscale;

	bool gtenable;
	int internal;
	long logml;

	// data loading functions:
	//void setupRates(PyObject *opt);
	void setupRates();

	void internal_set_stack_energies(FILE *fp, char *buffer);
	void internal_set_stack_enthalpies(FILE *fp, char *buffer);
	void internal_set_hairpin_energies(FILE *fp, char *buffer);
	void internal_set_hairpin_enthalpies(FILE *fp, char *buffer);
	void internal_set_hairpin_mismatch_energies(FILE *fp, char *buffer);
	void internal_set_hairpin_tetraloop_parameters(FILE *fp, char *buffer);
	void internal_set_hairpin_triloop_parameters(FILE *fp, char *buffer);
	void internal_set_hairpin_mismatch_enthalpies(FILE *fp, char *buffer);
	void internal_set_hairpin_tetraloop_parameters_enthalpy(FILE *fp,
			char *buffer);
	void internal_set_hairpin_triloop_parameters_enthalpy(FILE *fp,
			char *buffer);
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
	char *internal_read_array_data(FILE *fp, char *buffer, char* start_loc,
			int *read_loc, int size);
	char *internal_read_array_data(FILE *fp, char *buffer, char* start_loc,
			double *read_loc, int size);
	void internal_set_dangle_5_energies(FILE *fp, char *buffer);
	void internal_set_dangle_3_energies(FILE *fp, char *buffer);
	void internal_set_dangle_5_enthalpies(FILE *fp, char *buffer);
	void internal_set_dangle_3_enthalpies(FILE *fp, char *buffer);

	double setWaterDensity(double temp);
};

#endif /* __ENERGYMODEL_H__ */
