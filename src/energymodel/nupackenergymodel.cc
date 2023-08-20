/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "simoptions.h"
#include "options.h"
#include "basetype.h"


#include <time.h>

// For json
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/filereadstream.h"
#include <iostream>

#undef DEBUG
//#define DEBUG

using namespace rapidjson;


// FD: Jan 2018, moved dH over to doubles
static double T_scale(double dG, double dH, double T) {

	return (dH + T * (dG - dH) / 310.15);

}

const double CELSIUS37_IN_KELVIN = 310.15;
const double TEMPERATURE_ZERO_CELSIUS_IN_KELVIN = 273.15;

extern  int pairs[5];
extern  int pairtypes[5][5];
extern  int basepair_sw[8];

// helper function to convert to numerical base format.
extern BaseType baseLookup(char base);

double NupackEnergyModel::returnRate(double start_energy, double end_energy,
									 int enth_entr_toggle) {
	if (inspection) {
		return 1.0;
	}

	double dE = end_energy - start_energy;

	if (enth_entr_toggle == 3) {
		return biscale * exp(-(dE - dG_assoc) / _RT);
	}

	// dG_assoc, if it were included in (start_energy, end_energy), would need
	// to be deleted here. However, it never gets added into any energies except
	// for display purposes. So it gets used in the join move rate, but not
	// here.
	//
	// OLD: dG_assoc is typically a negative number, and included as part of the
	// complex before disassociation. Thus it must be subtracted from the dE
	// (leading to a typically slower disassociation rate.).
	if (kinetic_rate_method == RATE_METHOD_KAWASAKI) {

		return uniscale * exp(-0.5 * dE / _RT);

	} else if (kinetic_rate_method == RATE_METHOD_METROPOLIS or
			   kinetic_rate_method == RATE_METHOD_ARRHENIUS) {
		// Metropolis
		if (dE < 0) {
			return uniscale * 1.0;
		} else {
			return uniscale * exp(-dE / _RT);
		}
	}

	fprintf(stderr, "ERROR: Invalid rate method in NupackEnergyModel::returnRate()\n");
	exit(1);

}

double NupackEnergyModel::getJoinRate(void) {

	if (inspection) {
		return (double) 1.0;
	}

	return joinrate; // replace with the passed in rate
					 // joinrate includes biscale (via setupRates();)
}

double NupackEnergyModel::getJoinRate_NoVolumeTerm(void) {

	if (inspection) {
		return 1.0;
	}

	return biscale;
}

double NupackEnergyModel::getVolumeEnergy(void) {
	return dG_volume;
}

double NupackEnergyModel::getAssocEnergy(void) {
	return dG_assoc;
}

// non entropy/enthalpy energy functions
double NupackEnergyModel::StackEnergy(int i, int j, int p, int q) {

	return stack_37_dG[pairtypes[i][j] - 1][pairtypes[p][q] - 1];

}

// non entropy/enthalpy energy functions
double NupackEnergyModel::StackEnthalpy(int i, int j, int p, int q) {

	return stack_37_dH[pairtypes[i][j] - 1][pairtypes[p][q] - 1];

}

double NupackEnergyModel::BulgeEnergy(int i, int j, int p, int q, int bulgesize) {

	return this->BulgeEnergy(i, j, p, q, bulgesize, bulge_37_dG, stack_37_dG);

}

double NupackEnergyModel::BulgeEnthalpy(int i, int j, int p, int q, int bulgesize) {

	return this->BulgeEnergy(i, j, p, q, bulgesize, bulge_37_dH, stack_37_dH);
}

double NupackEnergyModel::BulgeEnergy(int i, int j, int p, int q, int bulgesize, array<double, 31> bulge,
		array<array<double, PAIRS_NUPACK>, PAIRS_NUPACK> stack) {

	double energy = 0.0;

	if (bulgesize <= 30) {

		energy = bulge[bulgesize];

	} else {

		energy = bulge[30] + (log((double) bulgesize / 30.0) * log_loop_penalty / 100.0);

	}

	if (bulgesize == 1) { // add stacking term for single-base bulges.

		energy += stack[pairtypes[i][j] - 1][pairtypes[p][q] - 1];

	} else { // AU penalty doesn't apply if they stack.

		if (pairtypes[i][j] == 1 || pairtypes[i][j] > 3) // AT penalty applies
			energy += terminal_AU;
		if (pairtypes[q][p] == 1 || pairtypes[q][p] > 3) // AT penalty applies
			energy += terminal_AU;

	}

	return energy;
}

double NupackEnergyModel::InteriorEnergy(BaseType *seq1, BaseType *seq2, int size1, int size2) {

	return this->InteriorEnergy(seq1, seq2, size1, size2, internal_dG);

}

double NupackEnergyModel::InteriorEnthalpy(BaseType *seq1, BaseType *seq2, int size1, int size2) {

	return this->InteriorEnergy(seq1, seq2, size1, size2, internal_dH);

}

double NupackEnergyModel::InteriorEnergy(BaseType *seq1, BaseType *seq2, int size1, int size2, internal_energies& internal) {

	double energy, ninio;
	int type1 = pairtypes[seq1[0]][seq2[size2 + 1]] - 1;
	int type2 = pairtypes[seq1[size1 + 1]][seq2[0]] - 1;

    if(type1 == -1 || type2 == -1){
        fprintf(stderr, "ERROR: Initial structure contains invalid bindings.\n");
		exit(1);
	}
	// special case time. 1x1, 2x1 and 2x2's all get special cases.
	if (size1 == 1 && size2 == 1){
		return internal.internal_1_1[type1][type2][seq1[1]][seq2[size2]]; //seq1[1][0] seq2[1][0]
	}
	if (size1 <= 2 && size2 <= 2)
		if (size1 == 1 || size2 == 1) {
			if (size1 == 1)
				return internal.internal_2_1[type1][seq1[1]][type2][seq2[1]][seq2[size2]];
			else
				return internal.internal_2_1[basepair_sw_mfold_actual[type2 + 1] - 1][seq2[1]][basepair_sw_mfold_actual[type1 + 1] - 1][seq1[1]][seq1[size1]];
		}
	if (size1 == 2 && size2 == 2)
		return internal.internal_2_2[type1][type2][seq1[1]][seq1[size1]][seq2[1]][seq2[size2]];

	// Generic case.

	if (size1 + size2 <= 30) {
		energy = internal.basic[size1 + size2];
	} else {
		energy = internal.basic[30] + (log((double) (size1 + size2) / 30.0) * log_loop_penalty / 100.0);
	}

	// NINIO term...
	int asym;
	if (size1 < size2) {
		asym = size1;
	} else {
		asym = size2;
	}
	if (asym > 4) {
		asym = 4;
	}

	ninio = abs(size2 - size1) * internal.ninio_correction[asym - 1];

	if (internal.maximum_NINIO < ninio) {
		energy += internal.maximum_NINIO;
	} else {
		energy += ninio;
	}

	// try gail params?
	if (size1 == 1 || size2 == 1) {
		energy += internal.mismatch[1][1][type1] + internal.mismatch[1][1][basepair_sw_mfold_actual[type2 + 1] - 1];
	} else {
		energy += internal.mismatch[seq1[1]][seq2[size2]][type1] + internal.mismatch[seq2[1]][seq1[size1]][basepair_sw_mfold_actual[type2 + 1] - 1];
	}

	// FD: adding the singlestranded stacking term.
	energy += singleStrandedStacking(seq1, size1); // TODO: decide if these are dH or dS contributions.
	energy += singleStrandedStacking(seq2, size2); // they are considered dH right now.

	return energy;
}

double NupackEnergyModel::HairpinEnergy(BaseType *seq, int size) {

	return HairpinEnergy(seq, size, hairpin_dG);

}

double NupackEnergyModel::HairpinEnthalpy(BaseType *seq, int size) {

	return HairpinEnergy(seq, size, hairpin_dH);

}

double NupackEnergyModel::HairpinEnergy(BaseType *seq, int size, hairpin_energies& hairpin) {

	double energy = 0.0;
	int lookup_index = 0;

	if (size <= 30) {
		energy = hairpin.basic[size];
	} else {
		energy = hairpin.basic[30];
		energy += (log((double) size / 30.0) * log_loop_penalty / 100.0);

	}

	if (size == 3) { // triloop bonuses

		// We now always do the lookup, if no entry then it is 0.0.
		lookup_index = ((seq[0] - 1) << 8) + ((seq[1] - 1) << 6) + ((seq[2] - 1) << 4) + ((seq[3] - 1) << 2) + (seq[4] - 1);

		energy += hairpin.triloop[lookup_index];

		if ((seq[0] == baseT) || (seq[size + 1] == baseT)) {
			energy += terminal_AU;
		}
	}

	if (size == 4) {
		lookup_index = ((seq[0] - 1) << 10) + ((seq[1] - 1) << 8) + ((seq[2] - 1) << 6) + ((seq[3] - 1) << 4) + ((seq[4] - 1) << 2) + (seq[5] - 1);

		energy += hairpin.tetraloop[lookup_index];
	}

	if (size >= 4)
		energy += hairpin.mismatch[(pairtypes[seq[0]][seq[size + 1]] - 1)][seq[1]][seq[size]];

	// FD: single stranded stacks.
	energy += singleStrandedStacking(seq, size);

	return energy;
}

double NupackEnergyModel::MultiloopEnergy(int size, int *sidelen, BaseType **sequences) {

	return MultiloopEnergy(size, sidelen, sequences, multiloop_dG);

}

double NupackEnergyModel::MultiloopEnthalpy(int size, int *sidelen, BaseType **sequences) {

	return MultiloopEnergy(size, sidelen, sequences, multiloop_dH);

}

double NupackEnergyModel::MultiloopEnergy(int size, int *sidelen, BaseType **sequences, multiloop_energies& multiloop) {

	// no dangle terms yet, this is equiv to dangles = 0;
	int totallength = 0;
	double energy = 0.0, dangle3, dangle5;
	int pt, rt_pt;
	int loopminus1 = size - 1;

	for (int loop = 0; loop < size; loop++) {

		totallength += sidelen[loop];
		pt = pairtypes[sequences[loopminus1][sidelen[loopminus1] + 1]][sequences[loop][0]] - 1;

		if ((pt == 0) || (pt > 2)) { // AT penalty applies
			energy += terminal_AU;
		}
		if (!gtenable && (pt > 3)) { // GT penalty applies
			energy += 100000.00;
		}

		loopminus1++;

		if (loopminus1 == size) {
			loopminus1 = 0;
		}

		// FD: single stranded stacks.
		energy += singleStrandedStacking(sequences[loop], sidelen[loop]);
		// FD: initialization of branch migration penalty.
		energy += initializationPenalty(sidelen[loop], loop, size);

	}

	if (simOptions->debug)
		cout << " Mid MultiLoop -- Energy is now " << energy << endl;

	energy += size * multiloop.internal;
	energy += multiloop.closing;

	if (!logml) {
		energy += multiloop.base * totallength;
	} else if (totallength <= 6) {
		energy += multiloop.base * totallength;
	} else {
		energy += multiloop.base * 6 + (log((double) totallength / 6.0) * log_loop_penalty / 100.0);
	}

	if (dangles == DANGLES_NONE) {

		return energy;

	} else {

		loopminus1 = size - 1;
		pt = pairtypes[sequences[loopminus1][0]][sequences[loopminus1 - 1][sidelen[loopminus1 - 1] + 1]] - 1;

		for (int loop = 0; loop < size; loop++) {

			rt_pt = pairtypes[sequences[loop][0]][sequences[loopminus1][sidelen[loopminus1] + 1]] - 1;

			if (!(dangles == DANGLES_SOME && sidelen[loopminus1] == 0)) {

				dangle5 = multiloop.dangle_5[pt][sequences[loopminus1][1]];
				dangle3 = multiloop.dangle_3[rt_pt][sequences[loopminus1][sidelen[loopminus1]]];

				if (dangles == DANGLES_SOME && sidelen[loopminus1] == 1) {
					energy += ((dangle3 < dangle5) ? dangle3 : dangle5); // minimum of two terms.

				} else {
					energy += dangle3 + dangle5;
				}
			}

			loopminus1++;
			if (loopminus1 == size) {
				loopminus1 = 0;
			}
			pt = rt_pt;
		}

		return energy;
	}

}

double NupackEnergyModel::OpenloopEnergy(int size, int *sidelen, BaseType **sequences) {
	return OpenloopEnergy(size, sidelen, sequences, multiloop_dG);

}

double NupackEnergyModel::OpenloopEnthalpy(int size, int *sidelen, BaseType **sequences) {

	return OpenloopEnergy(size, sidelen, sequences, multiloop_dH);

}

double NupackEnergyModel::OpenloopEnergy(int size, int *sidelen, BaseType **sequences, multiloop_energies& multiloop) {

	if (simOptions->debug)
		cout << "Computing OpenLoopEnergy, size = " << size << endl;

	// no dangle terms yet, this is equiv to dangles = 0;
	double energy = 0.0;
	int pt, loop, rt_pt;

	for (loop = 0; loop < size; loop++) {
		pt = pairtypes[sequences[loop][sidelen[loop] + 1]][sequences[loop + 1][0]] - 1;
		// TODO: slight efficiency gain if we wrap this into the dangles version separately, rather than doing a double pass in the dangle case.

		if ((pt == 0) || (pt > 2)) { // AT penalty applies
			energy += terminal_AU;
		}
		if (!gtenable && (pt > 3)) { // GT penalty applies
			energy += 100000.0;
		}
		// FD: adding singlestranded stacking.
		energy += singleStrandedStacking(sequences[loop], sidelen[loop]);
		// FD: initialization of branch migration penalty.
		energy += initializationPenalty(sidelen[loop], loop, size);

	}

	if (simOptions->debug)
		cout << " Mid OpenLoop -- Energy is now " << energy << endl;

	if (dangles == DANGLES_NONE || size == 0) {
		return energy;
	} else {

		double dangle3 = 0.0, dangle5 = 0.0;

		// 5' most sequence's dangle3 component.
		pt = pairtypes[sequences[1][0]][sequences[0][sidelen[0] + 1]] - 1;

		if (sidelen[0] == 0)
			dangle3 = 0.0;
		else
			dangle3 = multiloop.dangle_3[pt][sequences[0][sidelen[0]]];

		energy += dangle3; // added for either dangle version.

		if (simOptions->debug)
			cout << " Mid2 OpenLoop -- Energy is now " << energy << endl;

		for (loop = 0; loop < size - 1; loop++) {

			rt_pt = pairtypes[sequences[loop + 2][0]][sequences[loop + 1][sidelen[loop + 1] + 1]] - 1;
			dangle5 = multiloop.dangle_5[pt][sequences[loop + 1][1]];
			dangle3 = multiloop.dangle_3[rt_pt][sequences[loop + 1][sidelen[loop + 1]]];
			if (dangles == DANGLES_SOME && sidelen[loop + 1] == 1) {
				energy += (dangle3 < dangle5 ? dangle3 : dangle5); // minimum of the two terms.
			} else if (dangles == DANGLES_SOME && sidelen[loop + 1] == 0) {
				energy += 0.0; // dangles=DANGLES_SOME has no stacking when 0 bases between.
							   // dangles=DANGLES_ALL, however, does. Weird, eh?
			} else {
				energy += dangle3 + dangle5;
			}

			pt = rt_pt;
		}

		if (simOptions->debug)
			cout << " Mid3 OpenLoop -- Energy is now " << energy << endl;

		// 3' most sequence's dangle5 component.
		if (sidelen[size] == 0) {
			dangle5 = 0.0;
		} else {
			dangle5 = multiloop.dangle_5[pt][sequences[size][1]];
		}

		energy += dangle5; // added for either dangle version.
	}

	if (simOptions->debug)
		cout << " End OpenLoop -- Energy is now " << energy << endl;

	return energy;
}

// constructors, internal functions

NupackEnergyModel::NupackEnergyModel(void) :
	// Check references for this loop penalty term.
	log_loop_penalty_37(107.856), kinetic_rate_method(RATE_METHOD_KAWASAKI),
	bimolecular_penalty(1.96), kBoltzmann(.00198717), current_temp(310.15)
{
}

NupackEnergyModel::NupackEnergyModel(PyObject* energy_options) :
	NupackEnergyModel()
{
	initOptions(new PSimOptions(energy_options));
}

NupackEnergyModel::NupackEnergyModel(SimOptions* options) :
	NupackEnergyModel()
{
	initOptions(options);
}

void NupackEnergyModel::initOptions(SimOptions* options)
{
	simOptions = options;
	processOptions();
	if (simOptions->energyOptions->usingArrhenius()) {
		computeArrheniusRates(current_temp);
//		printPrecomputedArrRates();
	}
}


// returns a FILE pointer or prints an error message.
FILE* NupackEnergyModel::openFiles(char* nupackhome, string& paramPath, string& fileName, int select) {

	string& fullpath = paramFiles[select];
	FILE* fp = NULL;

	if (nupackhome == NULL) {
		fullpath = fileName;
	} else {
		fullpath = nupackhome;
		fullpath += paramPath;
		fullpath += fileName;
	}

	fp = fopen(fullpath.c_str(), "rt");

	if (fp == NULL) {
		fp = fopen(fileName.c_str(), "rt");
		if (fp == NULL) {

			string errorm = string("ERROR: nupack parameter file not found:") + string(nupackhome);
			errorm += paramPath + fileName + string(" or ") + fileName;
			errorm += string(" in current directory.\n");
			fprintf(stderr, "%s", errorm.c_str()); //fd: without "%s" potentially insecure, clang
			exit(0);
		}

	}

	return fp;

}

void NupackEnergyModel::processOptions() {

	int loop, loop2, loop3, loop4, loop5, loop6;
	double temperature;
	FILE *fp = NULL; // fp is json file.

	EnergyOptions* myEnergyOptions = simOptions->getEnergyOptions();

	// why use global variables when you can call an object?
	current_temp = myEnergyOptions->getTemperature();

	temperature = myEnergyOptions->getTemperature();
	dangles = myEnergyOptions->getDangles();
	logml = myEnergyOptions->getLogml();
	gtenable = myEnergyOptions->getGtenable();
	kinetic_rate_method = myEnergyOptions->getKineticRateMethod();

	for (loop = 0; loop < BASES; loop++)
		pairs[loop] = pairs_mfold[loop];
	for (loop = 0; loop < PAIRS_VIENNA; loop++)
		basepair_sw[loop] = basepair_sw_mfold[loop];
	for (loop = 0; loop < BASES; loop++)
		for (loop2 = 0; loop2 < BASES; loop2++)
			pairtypes[loop][loop2] = pairtypes_mfold[loop][loop2];

	if (myEnergyOptions->compareSubstrateType(SUBSTRATE_INVALID)) {

		PyObject *tmpStr = NULL;
		char* tmp = NULL;

		myEnergyOptions->getParameterFile(tmp, tmpStr);

		if (tmp != NULL) {
			fp = fopen(tmp, "rt");
			if (fp == NULL) {
				fprintf(stderr, "ERROR: Bad Parameter Filename: %s not found in path.\n", tmp);
				exit(1);
			}
			Py_DECREF(tmpStr);
			tmp = NULL;
		} else {
			fprintf(stderr, "ERROR: Invalid substrate chosen, and no parameter file given. Try the #Energymodel option!\n");
			exit(1);
		}

	} else if (myEnergyOptions->compareSubstrateType(SUBSTRATE_DNA) || myEnergyOptions->compareSubstrateType(SUBSTRATE_RNA)) {

		char *nupackhome;
		std::string file;

		nupackhome = getenv("NUPACKHOME");

		if(nupackhome == NULL){
		    fprintf(stderr, "ERROR: NUPACKHOME environment variable not set properly. Review setup instructions\n");
			exit(1);
		}

		if (myEnergyOptions->compareSubstrateType(SUBSTRATE_DNA)) {
			file = "dna04-nupack3.json";
		} else {
			file = "rna06-nupack3.json";
		}
		std::string paramPath = "/source/parameters/";
		fp = openFiles(nupackhome, paramPath, file, 0);
	}

    char buffer[65536];
    FileReadStream is(fp, buffer, sizeof(buffer));
    Document d;
    d.ParseStream(is);

    if (!d.HasMember("dH")) {
		if (temperature < CELSIUS37_IN_KELVIN - .0001 || temperature > CELSIUS37_IN_KELVIN + .0001) {
			fprintf(stderr,
					"ERROR: Temperature was set to %0.2lf C, but only dG type data files could be found. Please ensure that the requested parameter set has both dG and dH parameters!\n",
					temperature);
			exit(1);
		}
		return;
	}



    clock_t t;
    t = clock();

    internal_set_stack(d);
    internal_set_hairpin(d);
    internal_set_bulge(d);
    internal_set_interior_loop(d);
    internal_set_interior_1_1(d);
    internal_set_interior_2_1(d);
    internal_set_interior_2_2(d);
    internal_set_dangle_5(d);
    internal_set_dangle_3(d);
    internal_set_multiloop_parameters(d);
    internal_set_at_penalty(d);
    internal_set_bimolecular_penalty(d);
    internal_set_ninio_parameters(d);
    internal_set_hairpin_tetraloop_parameters(d);
    internal_set_hairpin_triloop_parameters(d);
    internal_set_hairpin_mismatch(d);
    internal_set_interior_loop_mismatch(d);
    fclose(fp);

// Temperature change section.

	_RT = kBoltzmann * current_temp;

	log_loop_penalty = 100.0 * 1.75 * kBoltzmann * current_temp;

	double saltCorrection = 0.368 * log(myEnergyOptions->sodium + 3.3 * sqrt(myEnergyOptions->magnesium));

	for (loop = 0; loop < PAIRS_NUPACK; loop++) {
		for (loop2 = 0; loop2 < PAIRS_NUPACK; loop2++) {

			stack_37_dG[loop][loop2] = T_scale(stack_37_dG[loop][loop2], stack_37_dH[loop][loop2], temperature);
			// now adjusting for a single salt correction term.
			stack_37_dG[loop][loop2] += (saltCorrection * -temperature) / 1000.0;
		}
	}

	for (loop = 0; loop < 31; loop++)
		hairpin_dG.basic[loop] = T_scale(hairpin_dG.basic[loop], hairpin_dH.basic[loop], temperature);

	for (loop = 0; loop < PAIRS_NUPACK; loop++)
		for (loop2 = 0; loop2 < BASES; loop2++)
			for (loop3 = 0; loop3 < BASES; loop3++)
				hairpin_dG.mismatch[loop][loop2][loop3] = T_scale(hairpin_dG.mismatch[loop][loop2][loop3], hairpin_dH.mismatch[loop][loop2][loop3],
						temperature);

	for (loop = 0; loop < 4096; loop++)
		hairpin_dG.tetraloop[loop] = T_scale(hairpin_dG.tetraloop[loop], hairpin_dH.tetraloop[loop], temperature);

	for (loop = 0; loop < 1024; loop++)
		hairpin_dG.triloop[loop] = T_scale(hairpin_dG.triloop[loop], hairpin_dH.triloop[loop], temperature);

	for (loop = 0; loop < 31; loop++)
		bulge_37_dG[loop] = T_scale(bulge_37_dG[loop], bulge_37_dH[loop], temperature);

	for (loop = 0; loop < 31; loop++)
		internal_dG.basic[loop] = T_scale(internal_dG.basic[loop], internal_dH.basic[loop], temperature);

	for (loop = 0; loop < BASES; loop++)
		for (loop2 = 0; loop2 < BASES; loop2++)
			for (loop3 = 0; loop3 < PAIRS_NUPACK; loop3++)
				internal_dG.mismatch[loop][loop2][loop3] = T_scale(internal_dG.mismatch[loop][loop2][loop3], internal_dH.mismatch[loop][loop2][loop3],
						temperature);

	internal_dG.maximum_NINIO = T_scale(internal_dG.maximum_NINIO, internal_dH.maximum_NINIO, temperature);

	for (loop = 0; loop < 5; loop++)
		internal_dG.ninio_correction[loop] = T_scale(internal_dG.ninio_correction[loop], internal_dH.ninio_correction[loop], temperature);

	for (loop = 0; loop < PAIRS_NUPACK; loop++)
		for (loop2 = 0; loop2 < PAIRS_NUPACK; loop2++)
			for (loop3 = 0; loop3 < BASES; loop3++)
				for (loop4 = 0; loop4 < BASES; loop4++)
					internal_dG.internal_1_1[loop][loop2][loop3][loop4] = T_scale(internal_dG.internal_1_1[loop][loop2][loop3][loop4],
							internal_dH.internal_1_1[loop][loop2][loop3][loop4], temperature);

	for (loop = 0; loop < PAIRS_NUPACK; loop++)
		for (loop5 = 0; loop5 < BASES; loop5++)
			for (loop2 = 0; loop2 < PAIRS_NUPACK; loop2++)
				for (loop3 = 0; loop3 < BASES; loop3++)
					for (loop4 = 0; loop4 < BASES; loop4++)
						internal_dG.internal_2_1[loop][loop5][loop2][loop3][loop4] = T_scale(internal_dG.internal_2_1[loop][loop5][loop2][loop3][loop4],
								internal_dH.internal_2_1[loop][loop5][loop2][loop3][loop4], temperature);

	for (loop = 0; loop < PAIRS_NUPACK; loop++)
		for (loop2 = 0; loop2 < PAIRS_NUPACK; loop2++)
			for (loop3 = 0; loop3 < BASES; loop3++)
				for (loop4 = 0; loop4 < BASES; loop4++)
					for (loop5 = 0; loop5 < BASES; loop5++)
						for (loop6 = 0; loop6 < BASES; loop6++)
							internal_dG.internal_2_2[loop][loop2][loop3][loop4][loop5][loop6] = T_scale(
									internal_dG.internal_2_2[loop][loop2][loop3][loop4][loop5][loop6],
									internal_dH.internal_2_2[loop][loop2][loop3][loop4][loop5][loop6], temperature);

	multiloop_dG.base = T_scale(multiloop_dG.base, multiloop_dH.base, temperature);
	multiloop_dG.closing = T_scale(multiloop_dG.closing, multiloop_dH.closing, temperature);
	multiloop_dG.internal = T_scale(multiloop_dG.internal, multiloop_dH.internal, temperature);

	for (loop = 0; loop < PAIRS_NUPACK; loop++)
		for (loop2 = 0; loop2 < BASES; loop2++) {
			multiloop_dG.dangle_3[loop][loop2] = T_scale(multiloop_dG.dangle_3[loop][loop2], multiloop_dH.dangle_3[loop][loop2], temperature);
			multiloop_dG.dangle_5[loop][loop2] = T_scale(multiloop_dG.dangle_5[loop][loop2], multiloop_dH.dangle_5[loop][loop2], temperature);
		}

	terminal_AU = T_scale(terminal_AU, terminal_AU_dH, temperature);

	bimolecular_penalty = T_scale(bimolecular_penalty, bimolecular_penalty_dH, temperature);
// need additional conversion as well

	_RT = kBoltzmann * temperature;

	current_temp = temperature;

	//FD: adding cotranscriptional initialziation
	numActiveNT = simOptions->initialActiveNT;

	setupRates();
}

/* ------------------------------------------------------------------------


 Private functions for use in loading and setting data by the main constructor


 ------------------------------------------------------------------------ */
void NupackEnergyModel::internal_set_stack(Document &d){
    int loop, loop2;
    Value &dG = d["dG"]["stack"],
          &dH = d["dH"]["stack"];
    for(loop = 0; loop < PAIRS_NUPACK; loop++){
        for(loop2 = 0; loop2 < PAIRS_NUPACK; loop2++){
            char name[] = {
                basepairString[loop + 1][0], basepairString[loop2 + 1][0],
                basepairString[loop2 + 1][2], basepairString[loop + 1][2], '\0'};
            stack_37_dG[loop][loop2] = dG[name].GetDouble();
            stack_37_dH[loop][loop2] = dH[name].GetDouble();
        }
    }
}

void NupackEnergyModel::internal_set_hairpin(Document &d) {
    Value &dG = d["dG"]["hairpin_size"].GetArray(),
          &dH = d["dH"]["hairpin_size"].GetArray();
    for(int i = 0; i < 30; i++){
        hairpin_dG.basic[i + 1] = dG[i].GetDouble();
        hairpin_dH.basic[i + 1] = dH[i].GetDouble();
    }
}

void NupackEnergyModel::internal_set_bulge(Document &d) {
	bulge_37_dG[0] = 0;
	bulge_37_dH[0] = 0;
    Value &dG = d["dG"]["bulge_size"].GetArray(),
          &dH = d["dH"]["bulge_size"].GetArray();
    for(int i = 0; i < 30; i++){
        bulge_37_dG[i + 1] = dG[i].GetDouble();
        bulge_37_dH[i + 1] = dH[i].GetDouble();
    }
}


void NupackEnergyModel::internal_set_interior_loop(Document &d) {
    internal_dG.basic[0] = 0;
    internal_dH.basic[0] = 0;
    Value &dG = d["dG"]["interior_size"].GetArray(),
          &dH = d["dH"]["interior_size"].GetArray();
    for(int i = 0; i < 30; i++){
        internal_dG.basic[i + 1] = dG[i].GetDouble();
        internal_dH.basic[i + 1] = dH[i].GetDouble();
    }
}

void NupackEnergyModel::internal_set_interior_1_1(Document &d) {
    // CG..AU CXAUYG order
    int loop, loop2, loop3, loop4;
    Value &dG = d["dG"]["interior_1_1"],
          &dH = d["dH"]["interior_1_1"];
    for (loop = 0; loop < PAIRS_NUPACK; loop++){
        for (loop2 = 0; loop2 < PAIRS_NUPACK; loop2++) {
            for (loop3 = 0; loop3 < BASES; loop3++) {
                internal_dG.internal_1_1[loop][loop2][loop3][0] = 0;
                internal_dG.internal_1_1[loop][loop2][0][loop3] = 0;

                internal_dH.internal_1_1[loop][loop2][loop3][0] = 0;
                internal_dH.internal_1_1[loop][loop2][0][loop3] = 0;
            }

            for (loop3 = 1; loop3 < BASES; loop3++) {
                for(loop4 = 1; loop4 < BASES; loop4++){
                    char name[] = {
                        basepairString[loop + 1][0], baseTypeString[loop3][0],
                    basepairString[loop2 + 1][0], basepairString[loop2 + 1][2],
                    baseTypeString[loop4][0], basepairString[loop + 1][2], '\0'};
                    internal_dG.internal_1_1[loop][loop2][loop3][loop4] = dG[name].GetDouble();
                    internal_dH.internal_1_1[loop][loop2][loop3][loop4] = dH[name].GetDouble();
                }
            }
        }
    }
}

void NupackEnergyModel::internal_set_interior_2_1(Document &d) {
    int loop, loop2, loop3, loop4, loop5;
    Value &dG = d["dG"]["interior_1_2"],
          &dH = d["dH"]["interior_1_2"];
    for (loop = 0; loop < PAIRS_NUPACK; loop++) {
        for (loop2 = 0; loop2 < PAIRS_NUPACK; loop2++) {
            for (loop3 = 1; loop3 < BASES; loop3++) {
                for (loop4 = 1; loop4 < BASES; loop4++) {
                    for(loop5 = 1; loop5 < BASES; loop5++){
                        char name[] = {
                            basepairString[loop + 1][0], baseTypeString[loop3][0],
                            basepairString[loop2 + 1][0], basepairString[loop2 + 1][2],
                            baseTypeString[loop4][0], baseTypeString[loop5][0],
                            basepairString[loop + 1][2], '\0'};
                        internal_dG.internal_2_1[loop][loop3][loop2][loop4][loop5] = dG[name].GetDouble();
                        internal_dH.internal_2_1[loop][loop3][loop2][loop4][loop5] = dH[name].GetDouble();
                    }
                }
            }
        }
    }
}

void NupackEnergyModel::internal_set_interior_2_2(Document &d) {
    Value &dG = d["dG"]["interior_2_2"],
          &dH = d["dH"]["interior_2_2"];
    int loop, loop2, loop3, loop4, loop5, loop6;
    for (loop = 0; loop < PAIRS_NUPACK; loop++) {
        for (loop2 = 0; loop2 < PAIRS_NUPACK; loop2++) {
            for (loop3 = 1; loop3 < BASES; loop3++) {
                for (loop4 = 1; loop4 < BASES; loop4++) {
                    for (loop5 = 1; loop5 < BASES; loop5++) {
                        for(loop6 = 1; loop6 < BASES; loop6++){ // eek 22,500 reads
                            char name[] = {
                                basepairString[loop + 1][0], baseTypeString[loop3][0],
                                baseTypeString[loop4][0], basepairString[loop2 + 1][0],
                                basepairString[loop2 + 1][2], baseTypeString[loop5][0],
                                baseTypeString[loop6][0], basepairString[loop + 1][2], '\0'};
                            internal_dG.internal_2_2[loop][loop2][loop3][loop4][loop5][loop6] = dG[name].GetDouble();
                            internal_dH.internal_2_2[loop][loop2][loop3][loop4][loop5][loop6] = dH[name].GetDouble();
                        }
                    }
                }
            }
        }
    }
}

void NupackEnergyModel::internal_set_dangle_5(Document &d) {
    // X2 X1 Y
    int loop, loop2;
    for (loop = 0; loop < PAIRS_NUPACK; loop++){
        for (loop2 = 0; loop2 < BASES; loop2++){
            multiloop_dG.dangle_5[loop][loop2] = 0.0;
            multiloop_dH.dangle_5[loop][loop2] = 0.0;
        }
    }
    Value &dG = d["dG"]["dangle_5"],
          &dH = d["dH"]["dangle_5"];
    for (loop = 0; loop < PAIRS_NUPACK; loop++){
         for (loop2 = 1; loop2 < BASES; loop2++){
            char name[] = {
                basepairString[loop + 1][2], basepairString[loop + 1][0],
                baseTypeString[loop2][0], '\0'};
            multiloop_dG.dangle_5[loop][loop2] = dG[name].GetDouble();
            multiloop_dH.dangle_5[loop][loop2] = dH[name].GetDouble();
         }
    }

}

void NupackEnergyModel::internal_set_dangle_3(Document &d) {
    // This reads as Y X2 X1
    int loop, loop2;
    for (loop = 0; loop < PAIRS_NUPACK; loop++){
        for (loop2 = 0; loop2 < BASES; loop2++){
            multiloop_dG.dangle_3[loop][loop2] = 0.0;
            multiloop_dH.dangle_3[loop][loop2] = 0.0;
        }
    }

    Value &dG = d["dG"]["dangle_3"],
          &dH = d["dH"]["dangle_3"];
    for (loop = 0; loop < PAIRS_NUPACK; loop++){
         for (loop2 = 1; loop2 < BASES; loop2++){
            char name[] = {
                baseTypeString[loop2][0], basepairString[loop + 1][2],
                basepairString[loop + 1][0], '\0'};
            if(!d["dG"]["dangle_3"].HasMember(name)){
                multiloop_dG.dangle_3[loop][loop2] = 0.0;
                multiloop_dH.dangle_3[loop][loop2] = 0.0;
            } else {
                multiloop_dG.dangle_3[loop][loop2] = dG[name].GetDouble();
                multiloop_dH.dangle_3[loop][loop2] = dH[name].GetDouble();
            }
         }
    }
}

void NupackEnergyModel::internal_set_multiloop_parameters(Document &d) {
    multiloop_dG.base = d["dG"]["multiloop_base"].GetDouble();
    multiloop_dG.closing = d["dG"]["multiloop_init"].GetDouble();
    multiloop_dG.internal = d["dG"]["multiloop_pair"].GetDouble();

    multiloop_dH.base = d["dH"]["multiloop_base"].GetDouble();
    multiloop_dH.closing = d["dH"]["multiloop_init"].GetDouble();
    multiloop_dH.internal = d["dH"]["multiloop_pair"].GetDouble();

}

void NupackEnergyModel::internal_set_at_penalty(Document &d) {
    // Accounts for GT wobble pairs
    terminal_AU = d["dG"]["terminal_penalty"]["AT"].GetDouble();
    terminal_AU_dH = d["dH"]["terminal_penalty"]["AT"].GetDouble();
}

void NupackEnergyModel::internal_set_bimolecular_penalty(Document &d) {
    bimolecular_penalty = d["dG"]["join_penalty"].GetDouble();
    bimolecular_penalty_dH = d["dH"]["join_penalty"].GetDouble();
}

void NupackEnergyModel::internal_set_ninio_parameters(Document &d) {
    Value &dG = d["dG"]["asymmetry_ninio"].GetArray(),
          &dH = d["dH"]["asymmetry_ninio"].GetArray();
    for(int i = 0; i < 4; i++){
        internal_dG.ninio_correction[i] = dG[i].GetDouble();
        internal_dH.ninio_correction[i] = dH[i].GetDouble();
    }
    internal_dG.maximum_NINIO = dG[4].GetDouble();
    internal_dH.maximum_NINIO = dH[4].GetDouble();
}

void NupackEnergyModel::internal_set_hairpin_tetraloop_parameters(Document &d) {
    for (int loop = 0; loop < 4096; loop++) {
        hairpin_dG.tetraloop[loop] = 0.0;
        hairpin_dH.tetraloop[loop] = 0.0;
    }
    int lookup_index;
    for (Value::ConstMemberIterator itr = d["dG"]["hairpin_tetraloop"].MemberBegin();
         itr != d["dG"]["hairpin_tetraloop"].MemberEnd();
         ++itr)
    {
        string name = itr->name.GetString();
        lookup_index = ((baseLookup(name[0]) - 1) << 10) + ((baseLookup(name[1]) - 1) << 8)
            + ((baseLookup(name[2]) - 1) << 6) + ((baseLookup(name[3]) - 1) << 4)
            + ((baseLookup(name[4]) - 1) << 2) + (baseLookup(name[5]) - 1);
        hairpin_dG.tetraloop[lookup_index] = itr->value.GetDouble();
    }
    for (Value::ConstMemberIterator itr = d["dH"]["hairpin_tetraloop"].MemberBegin();
         itr != d["dH"]["hairpin_tetraloop"].MemberEnd();
         ++itr)
    {
        string name = itr->name.GetString();
        lookup_index = ((baseLookup(name[0]) - 1) << 10) + ((baseLookup(name[1]) - 1) << 8)
            + ((baseLookup(name[2]) - 1) << 6) + ((baseLookup(name[3]) - 1) << 4)
            + ((baseLookup(name[4]) - 1) << 2) + (baseLookup(name[5]) - 1);
        hairpin_dH.tetraloop[lookup_index] = itr->value.GetDouble();
    }
}

void NupackEnergyModel::internal_set_hairpin_triloop_parameters(Document &d) {
    // triloop
    for (int loop = 0; loop < 1024; loop++) {
        hairpin_dG.triloop[loop] = 0.0;
        hairpin_dH.triloop[loop] = 0.0;
    }
    int lookup_index = 0;
    for (Value::ConstMemberIterator itr = d["dG"]["hairpin_triloop"].MemberBegin();
         itr != d["dG"]["hairpin_triloop"].MemberEnd();
         ++itr)
    {
        string name = itr->name.GetString();
        lookup_index = ((baseLookup(name[0]) - 1) << 8)
            + ((baseLookup(name[1]) - 1) << 6) + ((baseLookup(name[2]) - 1) << 4)
            + ((baseLookup(name[3]) - 1) << 2) + (baseLookup(name[4]) - 1);
        hairpin_dG.triloop[lookup_index] = itr->value.GetDouble();
    }
    for (Value::ConstMemberIterator itr = d["dH"]["hairpin_triloop"].MemberBegin();
         itr != d["dH"]["hairpin_triloop"].MemberEnd();
         ++itr)
    {
        string name = itr->name.GetString();
        lookup_index = ((baseLookup(name[0]) - 1) << 8)
            + ((baseLookup(name[1]) - 1) << 6) + ((baseLookup(name[2]) - 1) << 4)
            + ((baseLookup(name[3]) - 1) << 2) + (baseLookup(name[4]) - 1);
        hairpin_dH.triloop[lookup_index] = itr->value.GetDouble();
    }
}

void NupackEnergyModel::internal_set_hairpin_mismatch(Document &d) {
    // hairpin mismatch possibly rework this
    //int j = ((baseLookup(name[2]) - 1) * 2) + (baseLookup(name[1]) - 1) - 3; // col
    // bottom right, upper right, upper left, bottom left
    int loop, loop2, loop3;
    for (loop = 0; loop < PAIRS_NUPACK; loop++){
        for (loop2 = 0; loop2 < BASES; loop2++){
            for (loop3 = 0; loop3 < BASES; loop3++){
                hairpin_dG.mismatch[loop][loop2][loop3] = 0;
                hairpin_dH.mismatch[loop][loop2][loop3] = 0;
            }
        }
    }
    for (Value::ConstMemberIterator itr = d["dG"]["hairpin_mismatch"].MemberBegin();
         itr != d["dG"]["hairpin_mismatch"].MemberEnd();
         ++itr)
    {
        string name = itr->name.GetString();
        int j = ((baseLookup(name[2]) - 1) * 2) + (baseLookup(name[1]) - 1) - 3; // col
        hairpin_dG.mismatch[j][baseLookup(name[3])][baseLookup(name[0])] = itr->value.GetDouble();
    }
    for (Value::ConstMemberIterator itr = d["dH"]["hairpin_mismatch"].MemberBegin();
         itr != d["dH"]["hairpin_mismatch"].MemberEnd();
         ++itr)
    {
        string name = itr->name.GetString();
        int j = ((baseLookup(name[2]) - 1) * 2) + (baseLookup(name[1]) - 1) - 3; // col
        hairpin_dH.mismatch[j][baseLookup(name[3])][baseLookup(name[0])] = itr->value.GetDouble();
    }
}

void NupackEnergyModel::internal_set_interior_loop_mismatch(Document &d) {
    // interior mismatch again maybe rework this
    int loop, loop2, loop3;
    for (loop = 0; loop < BASES; loop++){
        for (loop2 = 0; loop2 < BASES; loop2++){
            for (loop3 = 0; loop3 < PAIRS_NUPACK; loop3++){
                internal_dG.mismatch[loop][loop2][loop3] = 0;
                internal_dH.mismatch[loop][loop2][loop3] = 0;
            }
        }
    }
    for (Value::ConstMemberIterator itr = d["dG"]["interior_mismatch"].MemberBegin();
         itr != d["dG"]["interior_mismatch"].MemberEnd();
         ++itr)
    {
        string name = itr->name.GetString();
        int j = ((baseLookup(name[2]) - 1) * 2) + (baseLookup(name[1]) - 1) - 3; // col
        internal_dG.mismatch[baseLookup(name[3])][baseLookup(name[0])][j] = itr->value.GetDouble();
    }
    for (Value::ConstMemberIterator itr = d["dH"]["interior_mismatch"].MemberBegin();
         itr != d["dH"]["interior_mismatch"].MemberEnd();
         ++itr)
    {
        string name = itr->name.GetString();
        int j = ((baseLookup(name[2]) - 1) * 2) + (baseLookup(name[1]) - 1) - 3; // col
        internal_dH.mismatch[baseLookup(name[3])][baseLookup(name[0])][j] = itr->value.GetDouble();
    }
}

void NupackEnergyModel::setupRates() {
//void NupackEnergyModel::setupRates(PyObject *energy_options) {
// input concentration is in molar (M) units
// future versions of the parser will convert to these units.
//
// To get the energy difference (and thus the join rate)
// we must find the number of solvent molecules M_s in the volume
// which would contain a single molecule at this concentration.
// the energy term is then kT ln M_s + join term.

// Roughly, this conversion is as follows:
// Convert concentration to moles per cubic meter.
// Convert concentration to molecules per cubic meter
// invert to get cubic meters per molecule
// multiply by water density to get grams water per molecule
// divide by water molecular weight (g/mol) to get mols water per (strand) molecule
// multiply by avogadro's number to get molecules water per (strand) molecule
// This result is our M_s, the number of solvent molecules in our volume.

// Shortened form, to try and avoid floating point error:
// joinenergy = kT log( W / C)
// W = 55.6 mol/L (molarity of water)

// joinrate = exp( -dG / kT )
//          = exp( - kT log( W / C) / kT)
//          = exp( -log( W / C))
//          = C / W
	EnergyOptions* eOptions = simOptions->getEnergyOptions();

	biscale = eOptions->getBiScale();
	uniscale = eOptions->getUniScale();

// Two components to the join rate, the dG_assoc from mass action, and the dG_volume term which is related to the concentration. We compute each individually, following my derivation for the volume term, and the method used in Nupack for the dG_assoc term.

// the concentration is now in M units, rather than the previous mM. The input parser converts to these.
	dG_volume = _RT * log(1.0 / eOptions->getJoinConcentration());

	dG_assoc = bimolecular_penalty; // already computed and scaled for water density.

	joinrate = biscale * eOptions->getJoinConcentration();

}
