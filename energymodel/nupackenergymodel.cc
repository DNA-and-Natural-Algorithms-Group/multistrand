/*
  Copyright (c) 2007-2010 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "../include/energymodel.h"
#include "../include/python_options.h"

#undef DEBUG
//#define DEBUG
#define T_scale( dG, dH, T ) ((double) ((((dG) - ((dH)/100.0)) * (T) / 310.15) + ((dH)/100.0)))

#define CELSIUS37_IN_KELVIN 310.15
#define TEMPERATURE_ZERO_CELSIUS_IN_KELVIN 273.15

extern int pairs[5];// = {0,0,0,0,0};
extern int pairtypes[5][5];// = {
//  {0,0,0,0,0},
//  {0,0,0,0,0},
//  {0,0,0,0,0},
//  {0,0,0,0,0},
//  {0,0,0,0,0}
//};
extern int basepair_sw[8];// = {0,0,0,0,0,0,0,0};


extern int lookuphelper[26];// = {1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,4,4,0,0,0,0,0};
//                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z

// helper function to convert to numerical base format.
extern int baseLookup( char base ); //
/*{
  char temp = toupper(base);
  if( temp < 'A' || temp > 'Z' )
  return base;
  return lookuphelper[temp-'A'];
  }*/


NupackEnergyModel::~NupackEnergyModel( void )
{
  // TODO: is anything allocated now? Don't think so, all arrays are static still.
  // nothing is allocated within an energy model.
}


double NupackEnergyModel::returnRate( double start_energy, double end_energy, int enth_entr_toggle)
{
  double dE = end_energy - start_energy;
  //  double dE = ((double) ( end_energy - start_energy )) / 100.;
  if( enth_entr_toggle == 3 ) // O/O delete, needs bi scale and full dE.
    return biscale * exp( -(dE) / _RT ); 
  // dG_assoc, if it were included in (start_energy, end_energy), would need to be deleted here. However, it never gets added into any energies except for display purposes. So it gets used in the join move rate, but not here.
  // OLD: dG_assoc is typically a negative number, and included as part of the complex before disassociation. Thus it must be subtracted from the dE (leading to a typically slower disassociation rate.).
  if( kinetic_rate_method == RATE_METHOD_KAWASAKI )  // Kawasaki
    return uniscale * exp( -0.5 * dE / _RT );
  else if( kinetic_rate_method == RATE_METHOD_METROPOLIS ) // Metropolis
    if( dE < 0 )
      return uniscale * 1;
    else
      return uniscale * exp( -dE / _RT );
}

double NupackEnergyModel::getJoinRate( void )
{
  //  printf("%lf\n",joinrate);
  return joinrate; // replace with the passed in rate
                   // joinrate includes biscale (via setupRates();)
}

double NupackEnergyModel::getJoinRate_NoVolumeTerm( void )
{
  return biscale;
}

double NupackEnergyModel::getVolumeEnergy( void )
{
  return dG_volume;
}

double NupackEnergyModel::getAssocEnergy( void )
{
  return dG_assoc;
}

// entropy/enthalpy energy parameters

void NupackEnergyModel::eStackEnergy( int type1, int type2, energyS *energy )
{
  energy->dH = stack_37_dH[type1][basepair_sw[type2]];
  energy->nTdS = stack_37_dG[type1][basepair_sw[type2]] - energy->dH;
  //  return stack_37_dG[type1][basepair_sw[type2]];
}

/* change these to ent/enth versions.
   int NupackEnergyModel::BulgeEnergy( int type1, int type2, int bulgesize )
   {
   int energy = 0;
   if( bulgesize <= 30 )
   energy = bulge_37_dG[bulgesize];
   else
   energy = bulge_37_dG[30] + (int) (log((double)bulgesize / 30.0) * log_loop_penalty);

   if( bulgesize == 1) // add stacking term for single-base bulges.
   energy += stack_37_dG[type1][basepair_sw[type2]];
   else // AU penalty doesn't apply if they stack.
   {
   if( ptype == VIENNA )
   {
   if( type1 > 2 ) // AU penalty applies
   energy += terminal_AU;
   if( type2 > 2 ) // AU penalty applies
   energy += terminal_AU;
   }
   }

   return energy;
   }

   int NupackEnergyModel::InteriorEnergy( int type1, int type2, int size1, int size2, char mismatch[4] )
   {
   double energy,ninio;
   //  for( int loop = 0; loop < 4 ; loop++ )
   //  mismatch[loop] = baseLookup(mismatch[loop]);

   // special case time. 1x1, 2x1 and 2x2's all get special cases.
   if( size1 == 1 && size2 == 1)
   return internal_1_1_37_dG[type1][basepair_sw[type2]][mismatch[0]][mismatch[1]];
   if( size1 <= 2 && size2 <= 2)
   if( size1 == 1 || size2 == 1)
   {
   if( size1 == 1)
   return internal_2_1_37_dG[type1][basepair_sw[type2]][mismatch[0]][mismatch[3]][mismatch[1]];
   else
   {
   if( ptype == MFOLD)
   return internal_2_1_37_dG[basepair_sw_mfold_actual[type2]][basepair_sw_mfold_actual[type1]][mismatch[3]][mismatch[0]][mismatch[2]];
   else
   return internal_2_1_37_dG[basepair_sw[type2]][type1][mismatch[3]][mismatch[0]][mismatch[2]];
   }
   }
   if( size1 == 2 && size2 == 2)
   return internal_2_2_37_dG[type1][basepair_sw[type2]][mismatch[0]][mismatch[2]][mismatch[3]][mismatch[1]];


   // Generic case.

   if( size1 + size2 <= 30 )
   energy = internal_37_dG[ size1 + size2 ];
   else
   energy = internal_37_dG[30] + (int) (log((double)(size1+size2) / 30.0) * log_loop_penalty);

   // NINIO term... no idea why this is used.
   ninio = abs(size2 - size1) * ninio_correction_37[2]; // don't ask me why, this is according to the Vienna energy model and this term is in it. 
  
   if( maximum_NINIO < ninio )
   energy += maximum_NINIO;
   else
   energy += ninio;

   energy += internal_mismatch_37_dG[type1][mismatch[0]][mismatch[1]] + 
   internal_mismatch_37_dG[basepair_sw[type2]][mismatch[3]][mismatch[2]];

   return energy;

   }

   int NupackEnergyModel::HairpinEnergy( int type1, int size, char *special )
   {
   double energy;
   char cmpstring[7] = {0,0,0,0,0,0,0};
   char *temp;
   if( size <= 30 )
   energy = hairpin_37_dG[size];
   else
   energy = hairpin_37_dG[30] + (log((double)size / 30.0) * log_loop_penalty / 100.0);
  
   if( size == 3 ) // triloop bonuses
   {
   strncpy( cmpstring, special, size+2 );
   temp = strstr( hairpin_triloops  ,cmpstring);
   if( temp != NULL )
   energy += hairpin_triloop_37_dG[(temp - hairpin_triloops ) / 6];
   if( ptype == VIENNA && type1 > 2) // pairs which qualify for the AT penalty
   energy += terminal_AU;
   else if( ptype == MFOLD && ((type1 > 3) || (type1 == 1)))
   energy += terminal_AU;

   }
   if( size == 4)
   {
   strncpy( cmpstring, special, size+2 );
   temp = strstr( hairpin_tetraloops  ,cmpstring);
   if( temp != NULL )
   energy += hairpin_tetraloop_37_dG[(temp - hairpin_tetraloops ) / 7];
   }
   if( size >= 4 )
   energy += hairpin_mismatch_37_dG[type1][special[1]][special[size]];

   return energy;
   }

   int NupackEnergyModel::MultiloopEnergy( int size, int *pairtypes, int *sidelen, char **sequences)
   {
   // no dangle terms yet, this is equiv to dangles = 0;
   int totallength=0,energy=0;
   for( int loop = 0; loop < size; loop++)
   {
   totallength += sidelen[loop];
   if( pairtypes[loop] > 2)
   energy += terminal_AU;
   }
   energy += size * multiloop_internal;
   energy += multiloop_closing;
   if( totallength <= 6 )
   energy += multiloop_base * totallength;
   else
   energy += multiloop_base * 6 + (int) (log((double)totallength / 6.0) * log_loop_penalty);

   return energy;
   }


   int NupackEnergyModel::OpenloopEnergy( int size, int *pairtypes, int *sidelen, char **sequences)
   {
   // no dangle terms yet, this is equiv to dangles = 0;
   int energy=0;
   for( int loop = 0; loop < size; loop++)
   {
   if( pairtypes[loop] > 2)
   energy += terminal_AU;
   }
   return energy;
   }
*/

// non entropy/enthalpy energy functions
double NupackEnergyModel::StackEnergy( int i, int j, int p, int q)
{
  return stack_37_dG[pairtypes[i][j]-1][pairtypes[p][q]-1];
}


double NupackEnergyModel::BulgeEnergy( int i, int j, int p, int q, int bulgesize )
{
  double energy = 0;
  if( bulgesize <= 30 )
    energy = bulge_37_dG[bulgesize];
  else
    energy = bulge_37_dG[30] + (log((double)bulgesize / 30.0) * log_loop_penalty / 100.0);

  if( bulgesize == 1) // add stacking term for single-base bulges.
    energy += stack_37_dG[pairtypes[i][j]-1][pairtypes[p][q]-1];
  else // AU penalty doesn't apply if they stack.
    {
      if( pairtypes[i][j] == 1 || pairtypes[i][j] > 3) // AT penalty applies
		energy += terminal_AU;
      if( pairtypes[q][p] == 1 || pairtypes[q][p] > 3) // AT penalty applies
		energy += terminal_AU;
    }

  return energy;
}

double NupackEnergyModel::InteriorEnergy( char *seq1, char *seq2, int size1, int size2 )
{
  double energy,ninio;
  //  for( int loop = 0; loop < 4 ; loop++ )
  //  mismatch[loop] = baseLookup(mismatch[loop]);
  int type1 = pairtypes[seq1[0]][seq2[size2+1]] -1;
  int type2 = pairtypes[seq1[size1+1]][seq2[0]] -1;

  // mismatch0 = seq1[1]
  // mismatch1 = seq2[size2]
  // mismatch2 = seq1[size1]
  // mismatch3 = seq2[1]

  // special case time. 1x1, 2x1 and 2x2's all get special cases.
  if( size1 == 1 && size2 == 1)
    return internal_1_1_37_dG[type1][type2][seq1[1]][seq2[size2]] ;
  if( size1 <= 2 && size2 <= 2)
	if( size1 == 1 || size2 == 1)
	  {
		if( size1 == 1)
		  return internal_2_1_37_dG[type1][seq1[1]][type2][seq2[1]][seq2[size2]];
		else
		  return internal_2_1_37_dG[basepair_sw_mfold_actual[type2+1]-1][seq2[1]][basepair_sw_mfold_actual[type1+1]-1][seq1[1]][seq1[size1]];
	  }
  if( size1 == 2 && size2 == 2)
    return internal_2_2_37_dG[type1][type2][seq1[1]][seq1[size1]][seq2[1]][seq2[size2]];


  // Generic case.

  if( size1 + size2 <= 30 )
    energy = internal_37_dG[ size1 + size2 ];
  else
    energy = internal_37_dG[30] + (log((double)(size1+size2) / 30.0) * log_loop_penalty/100.0);

  // NINIO term... 
  int asym;
  if ( size1 < size2)
    asym = size1;
  else
    asym = size2;
  if( asym > 4 )
    asym = 4;
  ninio = abs(size2 - size1) * ninio_correction_37[asym-1]; 

  if( maximum_NINIO < ninio )
    energy += maximum_NINIO;
  else
    energy += ninio;

  // TODO: derotate this in the input rather than the output. see ene.c from nupack and the parameter file for details.
  //    energy += internal_mismatch_37_dG[type1][seq2[size2]][seq1[1]] + 
  //  internal_mismatch_37_dG[basepair_sw_mfold_actual[type2+1]-1][seq1[size1]][seq2[1]];
  // OLD VERSION:
  //  energy += internal_mismatch_37_dG[type1][seq1[1]][seq2[size2]] + 
  //  internal_mismatch_37_dG[basepair_sw_mfold_actual[type2+1]-1][seq2[1]][seq1[size1]];

  // try gail params?
  if(size1 == 1 || size2 == 1)
    energy += internal_mismatch_37_dG[1][1][type1] + internal_mismatch_37_dG[1][1][basepair_sw_mfold_actual[type2+1]-1];
  else
    energy += internal_mismatch_37_dG[seq1[1]][seq2[size2]][type1] + internal_mismatch_37_dG[seq2[1]][seq1[size1]][basepair_sw_mfold_actual[type2+1]-1];

  // Yet another try to match the wacky parameters.
  //energy += internal_mismatch_37_dG[seq2[size2]][seq1[1]][type1] +
  //  internal_mismatch_37_dG[seq2[1]][seq1[size1]][type2];
  //  int type1 = pairtypes[seq1[0]][seq2[size2+1]] -1;
  //  int type2 = pairtypes[seq1[size1+1]][seq2[0]] -1;

  return energy;
}

double NupackEnergyModel::HairpinEnergy( char *seq, int size )
{
  double energy = 0.0;
  char cmpstring[7] = {0,0,0,0,0,0,0};
  char *temp;
  int type1 = pairtypes[seq[0]][seq[size+1]]-1;

  if( size <= 30 )
    energy = hairpin_37_dG[size];
  else
    energy = hairpin_37_dG[30] + (log((double)size / 30.0) * log_loop_penalty/100.0);
  
  if( size == 3 ) // triloop bonuses
    {
      strncpy( cmpstring, seq, size+2 );
      temp = strstr( hairpin_triloops  ,cmpstring);
      if( temp != NULL )
		energy += hairpin_triloop_37_dG[(temp - hairpin_triloops ) / 6];
      if((type1 > 2) || (type1 == 0))
		energy += terminal_AU;
    }
  if( size == 4)
    {
      strncpy( cmpstring, seq, size+2 );
      temp = strstr( hairpin_tetraloops  ,cmpstring);
      if( temp != NULL )
		energy += hairpin_tetraloop_37_dG[(temp - hairpin_tetraloops ) / 7];
    }
  if( size >= 4 )
    energy += hairpin_mismatch_37_dG[type1][seq[1]][seq[size]];

  return energy;
}

double NupackEnergyModel::MultiloopEnergy( int size, int *sidelen, char **sequences)
{
  // no dangle terms yet, this is equiv to dangles = 0;
  int totallength=0;
  double energy=0.0, dangle3, dangle5;
  int pt, rt_pt;
  int loopminus1 = size-1;

  // pairtypes[loop] = pairtypes[ ][ ]
  // pairtypes[loop] is the last of sequence[loop-1] and first of sequence[loop]
  // pairtypes[loop] = pairtypes[sequences[loop-1][sidelen[loop-1]+1]][sequences[loop][0]];
  // pairtypes[sequences[loop-1][sidelen[loop-1]+1]][sequences[loop][0]];

  for( int loop = 0; loop < size; loop++)
    {
      totallength += sidelen[loop];
      pt = pairtypes[sequences[loopminus1][sidelen[loopminus1]+1]][sequences[loop][0]] -1;

      if((pt == 0) || (pt > 2)) // AT penalty applies
		energy += terminal_AU;
      if( !gtenable && (pt > 3)) // GT penalty applies
		energy += 100000.00;
      loopminus1++;
      if( loopminus1 == size ) loopminus1 = 0;
    }
  energy += size * multiloop_internal;
  energy += multiloop_closing;
  if( !logml )
    energy += multiloop_base * totallength;
  else if( totallength <= 6 )
    energy += multiloop_base * totallength;
  else
    energy += multiloop_base * 6 + (log((double)totallength / 6.0) * log_loop_penalty / 100.0);
  if( dangles == DANGLES_NONE )
    return energy;
  else
    {
      loopminus1 = size-1;
      //      pt = pairtypes[ sequences[ ][sidelen[ ] ][ sequences[][] ] - 1;
      pt = pairtypes[
					 sequences[loopminus1][0]
					 ][
					   sequences[loopminus1-1][sidelen[loopminus1-1]+1]
					   ] -1;
      for( int loop = 0; loop < size; loop++)
		{
		  rt_pt = pairtypes[
							sequences[loop][0]			  
							][
							  sequences[loopminus1][sidelen[loopminus1]+1]
							  ] -1;
	  
		  if( !(dangles == DANGLES_SOME && sidelen[loopminus1] == 0 ))
			{
			  dangle5 = dangle_5_37_dG[ pt ][ sequences[loopminus1][1] ];
			  dangle3 = dangle_3_37_dG[ rt_pt ][ sequences[loopminus1][sidelen[loopminus1]] ];
			  if (dangles == DANGLES_SOME && sidelen[loopminus1] == 1 )
				energy += ((dangle3<dangle5)?dangle3:dangle5); // minimum of two terms.
			  else
				energy += dangle3 + dangle5;
			}
	    
		  loopminus1++;
		  if( loopminus1 == size ) loopminus1 = 0;
		  pt = rt_pt;
		}
      
      return energy;
    }
}


double NupackEnergyModel::OpenloopEnergy( int size, int *sidelen, char **sequences)
{
  // no dangle terms yet, this is equiv to dangles = 0;
  double energy=0.0;
  int pt, loop, rt_pt;

  for( loop = 0; loop < size; loop++)
    {
      pt = pairtypes[sequences[loop][sidelen[loop]+1]][sequences[loop+1][0]] -1;
      // TODO: slight efficiency gain if we wrap this into the dangles version seperately, rather than doing a double pass in the dangle case.

      if((pt == 0) || (pt > 2)) // AT penalty applies
		energy += terminal_AU;
      if( !gtenable && (pt > 3)) // GT penalty applies
		energy += 100000.0;
    }
  if( dangles == DANGLES_NONE || size == 0)
    return energy;
  else
    {
      double dangle3=0.0,dangle5=0.0;

      // 5' most sequence's dangle3 component. 
      pt = pairtypes[ sequences[1][0] ][ sequences[0][sidelen[0]+1] ] -1;
      
      if( sidelen[0] == 0 )
		dangle3 = 0.0;
      else	  
		dangle3 = dangle_3_37_dG[ pt ][ sequences[0][sidelen[0]] ];

      energy += dangle3; // added for either dangle version.

      for( loop = 0; loop < size-1; loop++ )
		{
		  rt_pt = pairtypes[ sequences[loop+2][0] ][ sequences[loop+1][sidelen[loop+1]+1] ] - 1;
		  dangle5 = dangle_5_37_dG[ pt ][ sequences[loop+1][1] ];
		  dangle3 = dangle_3_37_dG[ rt_pt ][ sequences[loop+1][sidelen[loop+1]] ];
		  if( dangles == DANGLES_SOME && sidelen[loop+1] == 1)
			energy += (dangle3<dangle5?dangle3:dangle5); // minimum of the two terms.
		  else if( dangles == DANGLES_SOME && sidelen[loop+1] == 0 )
			energy += 0; // dangles=DANGLES_SOME has no stacking when 0 bases between.
                         // dangles=DANGLES_ALL, however, does. Weird, eh?
		  else
			energy += dangle3 + dangle5;

		  pt = rt_pt;
		}

      // 3' most sequence's dangle5 component.
      if( sidelen[size] == 0 )
		dangle5 = 0.0;
      else
		dangle5 = dangle_5_37_dG[ pt ][ sequences[size][1] ];

      energy += dangle5; // added for either dangle version.
    }
  return energy;
}


// constructors, internal functions


NupackEnergyModel::NupackEnergyModel( PyObject *energy_options ) : log_loop_penalty_37(107.856) , kinetic_rate_method( RATE_METHOD_KAWASAKI) , bimolecular_penalty(196), kBoltzmann(.00198717),current_temp(310.15) // Check references for this loop penalty term.
{
  // This is the tough part, performing all read/input duties.
  char in_buffer[2048];
  int loop, loop2,loop3, loop4, loop5, loop6 ;
  double temperature;
  FILE *fp = NULL, *fp2 = NULL; // fp is dG energy file, fp2 is dH.

  getDoubleAttr(energy_options, temperature,&temperature);
  current_temp =  temperature;

  getLongAttr(energy_options, dangles,&dangles);
  //  getLongAttr(energy_options, parameter_type,&ptype);
  getLongAttr(energy_options, log_ml,&logml);
  getBoolAttr(energy_options, gt_enable ,&gtenable);
  getLongAttr(energy_options, rate_method,&kinetic_rate_method);

  //  if( logml == -1 ) logml = 0;
  // JS: Should now be correct via new python options object.
  
  waterdensity = setWaterDensity( temperature - TEMPERATURE_ZERO_CELSIUS_IN_KELVIN);


  //  printf("Loop Pen: %lf %lf\n",log_loop_penalty_37, log_loop_penalty );

  for( loop = 0; loop < NUM_BASES; loop++ )
    pairs[loop] = pairs_mfold[loop];
  for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
    basepair_sw[loop] = basepair_sw_mfold[loop];
  for( loop = 0; loop < NUM_BASES; loop++)
    for( loop2 = 0; loop2 < NUM_BASES; loop2++)
      pairtypes[loop][loop2] = pairtypes_mfold[loop][loop2];

  if( testLongAttr(energy_options, substrate_type ,=, SUBSTRATE_INVALID ) )
    {
      PyObject *tmpStr=NULL;  
	  
      char *tmp = (char *) getStringAttr( energy_options, parameter_file,tmpStr);
	  if( tmp != NULL )
		{
		  fp = fopen( tmp, "rt");
		  if( fp == NULL )
			{
			  fprintf(stderr,"ERROR: Bad Parameter Filename: %s not found in path.\n", tmp );
			  exit(1);
			}
          Py_DECREF( tmpStr );
          tmp = NULL;
		}
	  else
		{
		  fprintf(stderr,"ERROR: Invalid substrate chosen, and no parameter file given. Try the #Energymodel option!\n");
		  exit(1);
		}
    }
  else if( testLongAttr(energy_options, substrate_type ,=, SUBSTRATE_DNA ) )
    {
      char *nupackhome;
      char fullpath[512];
      char fullpath_dH[512];
      
      nupackhome = getenv("NUPACKHOME");
     
      if( nupackhome == NULL )
		{
		  strcpy(fullpath,"DNA_mfold2.3.dG");
		  strcpy(fullpath_dH,"DNA_mfold2.3.dH");
		}
      else
		{
		  strcpy(fullpath,nupackhome);
		  strcpy(fullpath_dH,nupackhome);
		  strcat(fullpath,"/parameters/DNA_mfold2.3.dG");
		  strcat(fullpath_dH,"/parameters/DNA_mfold2.3.dH");
		}
      //printf("fp: %s\n fpdh: %s\n",fullpath, fullpath_dH);

      fp = fopen( fullpath, "rt");
      fp2 = fopen(fullpath_dH,"rt");
      
      if( fp == NULL )
		{
		  fp = fopen("DNA_mfold2.3.dG","rt");
		  if( fp == NULL )
			{
			  fprintf(stderr,"ERROR: nupack/mfold parameter file not found: $NUPACKHOME/parameters/DNA_mfold2.3.dG or DNA_mfold2.3.dG in current directory.\n");
			  exit(0);
			}

		}

      if( fp2 == NULL )
		{
		  fp2 = fopen("DNA_mfold2.3.dH","rt");
		  if( fp2 == NULL )
			{
			  fprintf(stderr,"ERROR: nupack/mfold parameter file not found: $NUPACKHOME/parameters/DNA_mfold2.3.dH or DNA_mfold2.3.dH in current directory.\n");
			  exit(0);
			}

		}
    }
  else if( testLongAttr(energy_options, substrate_type ,=, SUBSTRATE_RNA ) )
    {
      char *nupackhome;
      char fullpath[512];
      char fullpath_dH[512];
      
      nupackhome = getenv("NUPACKHOME");
      //      printf("nupackhome = %s\n",nupackhome);
     
      if( nupackhome == NULL )
		{
		  strcpy(fullpath,"RNA_mfold2.3.dG");
		  strcpy(fullpath_dH,"RNA_mfold2.3.dH");
		}
      else
		{
		  strcpy(fullpath,nupackhome);
		  strcpy(fullpath_dH,nupackhome);
		  strcat(fullpath,"/parameters/RNA_mfold2.3.dG");
		  strcat(fullpath_dH,"/parameters/RNA_mfold2.3.dH");
		}
      //      printf("fp: %s\n fpdh: %s\n",fullpath, fullpath_dH);

      fp = fopen( fullpath, "rt");
      fp2 = fopen(fullpath_dH,"rt");
      
      if( fp == NULL )
		{
		  fp = fopen("RNA_mfold2.3.dG","rt");
		  if( fp == NULL )
			{
			  fprintf(stderr,"ERROR: nupack/mfold parameter file not found: $NUPACKHOME/parameters/RNA_mfold2.3.dG or RNA_mfold2.3.dG in current directory.\n");
			  exit(0);
			}

		}

      if( fp2 == NULL )
		{
		  fp2 = fopen("RNA_mfold2.3.dH","rt");
		  if( fp2 == NULL )
			{
			  fprintf(stderr,"ERROR: nupack/mfold parameter file not found: $NUPACKHOME/parameters/RNA_mfold2.3.dH or RNA_mfold2.3.dH in current directory.\n");
			  exit(0);
			}

		}
    }
  
  
  fgets( in_buffer, 2048, fp );
  while( !feof(fp) )
    {
	  if( in_buffer[0] == '>' ) // data area or comment (mfold)
	    {
	      if( strncmp( in_buffer, ">Stacking 5' X1 Y1 3'", 21) == 0)
			{
#ifdef DEBUG
			  printf("Loading Stack Energies (MFOLD).\n");
#endif
			  internal_set_stack_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			  if( feof(fp) )
				continue;
			}
	      
	      else if( strncmp( in_buffer, ">Hairpin Loop Energies:", 23) == 0)
			{
#ifdef DEBUG
			  printf("Loading Hairpin Energies (MFOLD).\n");
#endif
			  internal_set_hairpin_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			  if( feof(fp) )
				continue;
			}
	      
	      else if( strncmp( in_buffer, ">Bulge loop Energies:", 21) == 0)
			{
#ifdef DEBUG
			  printf("Loading Bulge Energies (MFOLD).\n");
#endif
			  internal_set_bulge_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			  if( feof(fp) )
				continue;	    
			}
	      else if( strncmp( in_buffer, ">Interior Loop Energies:", 24) == 0)
			{
#ifdef DEBUG
			  printf("Loading Internal Loop Energies (MFOLD).\n");
#endif
			  internal_set_interior_loop_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			  if( feof(fp) )
				continue;	    
			}
	      
	      else if( strncmp( in_buffer, ">NINIO asymmetry", 16) == 0)
			{
#ifdef DEBUG
			  printf("Loading NINIO parameters (MFOLD).\n");
#endif
			  internal_set_ninio_parameters( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			  if( feof(fp) )
				continue;	    
			}
	      
	      else if( strncmp( in_buffer, ">Triloops ", 10) == 0)
			{
#ifdef DEBUG
			  printf("Loading Hairpin Triloop parameters (MFOLD).\n");
#endif
			  internal_set_hairpin_triloop_parameters( fp, in_buffer );
			}
	      else if( strncmp( in_buffer, ">Tetraloops ", 12) == 0)
			{
#ifdef DEBUG
			  printf("Loading Hairpin Tetraloop parameters (MFOLD).\n");
#endif
			  internal_set_hairpin_tetraloop_parameters( fp, in_buffer );
			}
	      else if( strncmp( in_buffer, ">Mismatch HP", 12) == 0)
			{
#ifdef DEBUG
			  printf("Loading Hairpin Mismatch Energies (MFOLD).\n");
#endif
			  internal_set_hairpin_mismatch_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			}	    
	      else if( strncmp( in_buffer, ">Mismatch Interior", 18) == 0)
			{
#ifdef DEBUG
			  printf("Loading Interior Loop Mismatch Energies (MFOLD).\n");
#endif
			  internal_set_interior_loop_mismatch_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			}
	      else if( strncmp( in_buffer, ">Dangle Energies: 5' X1 . 3'", 28) == 0)
			{
#ifdef DEBUG
			  printf("Loading Dangle 3' Energies (MFOLD).\n");
#endif
			  internal_set_dangle_3_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			}
	      else if( strncmp( in_buffer, ">Dangle Energies: 5' X1 Y 3'", 28) == 0)
			{
#ifdef DEBUG
			  printf("Loading Dangle 5' Energies (MFOLD).\n");
#endif
			  internal_set_dangle_5_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp );
			}
	      else if( strncmp( in_buffer, ">Multiloop terms:", 17) == 0)
			{
#ifdef DEBUG
			  printf("Loading Multiloop parameters (MFOLD).\n");
#endif
			  internal_set_multiloop_parameters( fp, in_buffer );
			  fgets( in_buffer, 2048, fp);
			}
	      else if( strncmp( in_buffer, ">AT_PENALTY:", 12) == 0)
			{
#ifdef DEBUG
			  printf("Loading AT (AU) Penalty parameter (MFOLD).\n");
#endif
			  internal_set_at_penalty( fp, in_buffer );
			  fgets( in_buffer, 2048, fp);
			}
	      else if( strncmp( in_buffer, ">Interior Loops 1x1", 19) == 0)
			{
#ifdef DEBUG
			  printf("Loading Internal 1-1 mismatch Energies (MFOLD).\n");
#endif
			  internal_set_interior_1_1_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp);
			}
	      else if( strncmp( in_buffer, ">Interior Loops 2x2", 19) == 0)
			{
#ifdef DEBUG
			  printf("Loading Internal 2-2 mismatch Energies (MFOLD).\n");
#endif
			  internal_set_interior_2_2_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp);
			}
	      else if( strncmp( in_buffer, ">Interior Loops 1x2", 19) == 0)
			{
#ifdef DEBUG
			  printf("Loading Internal 2-1 mismatch Energies (MFOLD).\n");
#endif
			  internal_set_interior_2_1_energies( fp, in_buffer );
			  fgets( in_buffer, 2048, fp);
			}
	      else if( strncmp( in_buffer, ">BIMOLECULAR", 12) == 0)
			{
#ifdef DEBUG
			  printf("Loading Bimolecular Association Penalty (MFOLD).\n");
#endif
			  internal_set_bimolecular_penalty( fp, in_buffer );
			  fgets( in_buffer, 2048, fp);
			}
	      else
			{
			  fgets( in_buffer, 2048, fp );
			}
	      
	    }
	  else
	    fgets( in_buffer, 2048, fp );
	}
  fclose(fp);
  /* Enthalpy loading section for all those pesky dH terms. */
  
  if(fp2 == NULL )
    {
      if( temperature < CELSIUS37_IN_KELVIN - .0001 || temperature > CELSIUS37_IN_KELVIN + .0001)
		{
		  fprintf(stderr,"ERROR: Temperature was set to %0.2lf C, but only dG type data files could be found. Please ensure that the requested parameter set has both .dG and .dH files!\n",temperature);
		  exit(0);
		}
      return;
    }
  
  fgets( in_buffer, 2048, fp2 );
  while( !feof(fp2) )
    {
	  if( in_buffer[0] == '>' ) // data area or comment (mfold)
		{
		  if( strncmp( in_buffer, ">Stacking 5' X1 Y1 3'", 21) == 0)
			{
#ifdef DEBUG
			  printf("Loading Stack Enthalpies (MFOLD).\n");
#endif
			  internal_set_stack_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2 );
			  if( feof(fp2) )
				continue;
			}
		  
		  else if( strncmp( in_buffer, ">Hairpin Loop Energies:", 23) == 0)
			{
#ifdef DEBUG
			  printf("Loading Hairpin Enthalpies (MFOLD).\n");
#endif
			  internal_set_hairpin_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2 );
			  if( feof(fp2) )
				continue;
			}
		  
		  else if( strncmp( in_buffer, ">Bulge loop Energies:", 21) == 0)
			{
#ifdef DEBUG
			  printf("Loading Bulge Enthalpies (MFOLD).\n");
#endif
			  internal_set_bulge_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2);
			  if( feof(fp2) )
				continue;	    
			}
		  else if( strncmp( in_buffer, ">Interior Loop Energies:", 24) == 0)
			{
#ifdef DEBUG
			  printf("Loading Internal Loop Enthalpies (MFOLD).\n");
#endif
			  internal_set_interior_loop_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2 );
			  if( feof(fp2) )
				continue;	    
			}
		  
		  else if( strncmp( in_buffer, ">NINIO asymmetry", 16) == 0)
			{
#ifdef DEBUG
			  printf("Loading NINIO parameters - enthalpy (MFOLD).\n");
#endif
			  internal_set_ninio_parameters_enthalpy( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2 );
			  if( feof(fp2) )
				continue;	    
			}
		  
		  else if( strncmp( in_buffer, ">Triloops ", 10) == 0)
			{
#ifdef DEBUG
			  printf("Loading Hairpin Triloop parameters - enthalpy (MFOLD).\n");
#endif
			  internal_set_hairpin_triloop_parameters_enthalpy( fp2, in_buffer );
			}
		  else if( strncmp( in_buffer, ">Tetraloops ", 12) == 0)
			{
#ifdef DEBUG
			  printf("Loading Hairpin Tetraloop parameters - enthalpy (MFOLD).\n");
#endif
			  internal_set_hairpin_tetraloop_parameters_enthalpy( fp2, in_buffer );
			}
		  else if( strncmp( in_buffer, ">Mismatch HP", 12) == 0)
			{
#ifdef DEBUG
			  printf("Loading Hairpin Mismatch Enthalpies (MFOLD).\n");
#endif
			  internal_set_hairpin_mismatch_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2 );
			}	    
		  else if( strncmp( in_buffer, ">Mismatch Interior", 18) == 0)
			{
#ifdef DEBUG
			  printf("Loading Interior Loop Mismatch Enthalpies (MFOLD).\n");
#endif
			  internal_set_interior_loop_mismatch_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2 );
			}
		  else if( strncmp( in_buffer, ">Dangle Energies: 5' X1 . 3'", 28) == 0)
			{
#ifdef DEBUG
			  printf("Loading Dangle 3' Enthalpies (MFOLD).\n");
#endif
			  internal_set_dangle_3_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2);
			}
		  else if( strncmp( in_buffer, ">Dangle Energies: 5' X1 Y 3'", 28) == 0)
			{
#ifdef DEBUG
			  printf("Loading Dangle 5' Enthalpies (MFOLD).\n");
#endif
			  internal_set_dangle_5_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2 );
			}
		  else if( strncmp( in_buffer, ">Multiloop terms:", 17) == 0)
			{
#ifdef DEBUG
			  printf("Loading Multiloop parameters (dH) (MFOLD).\n");
#endif
			  internal_set_multiloop_parameters_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2);
			}
		  else if( strncmp( in_buffer, ">AT_PENALTY:", 12) == 0)
			{
#ifdef DEBUG
			  printf("Loading AT (AU) Penalty parameter (dH) (MFOLD).\n");
#endif
			  internal_set_at_penalty_enthalpy( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2);
			}
		  else if( strncmp( in_buffer, ">Interior Loops 1x1", 19) == 0)
			{
#ifdef DEBUG
			  printf("Loading Internal 1-1 mismatch enthalpies (MFOLD).\n");
#endif
			  internal_set_interior_1_1_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2);
			}
		  else if( strncmp( in_buffer, ">Interior Loops 2x2", 19) == 0)
			{
#ifdef DEBUG
			  printf("Loading Internal 2-2 mismatch enthalpies (MFOLD).\n");
#endif
			  internal_set_interior_2_2_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2);
			}
		  else if( strncmp( in_buffer, ">Interior Loops 1x2", 19) == 0)
			{
#ifdef DEBUG
			  printf("Loading Internal 2-1 mismatch enthalpies (MFOLD).\n");
#endif
			  internal_set_interior_2_1_enthalpies( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2);
			}
		  else if( strncmp( in_buffer, ">BIMOLECULAR", 12) == 0)
			{
#ifdef DEBUG
			  printf("Loading Bimolecular Association Penalty (dH) (MFOLD).\n");
#endif
			  internal_set_bimolecular_penalty_dH( fp2, in_buffer );
			  fgets( in_buffer, 2048, fp2);
			}
		  else
			{
			  fgets( in_buffer, 2048, fp2 );
			}
		  
		}
	  else
		fgets( in_buffer, 2048, fp2 );
	}
  
  fclose(fp2);

  // Temperature change section.
  //  double getDoubleAttr(energy_options, temperature,&temperature);
  
  // Note: #define T_scale( dG, dH, T ) ((((dG) - (dH)) * (T) / 310.15) + dH)

  _RT = kBoltzmann * current_temp;
  log_loop_penalty = 100.0 * 1.75 * kBoltzmann * current_temp;

  if(!((temperature < CELSIUS37_IN_KELVIN - .00001) || (temperature > CELSIUS37_IN_KELVIN + .00001) ))
    {
      current_temp = CELSIUS37_IN_KELVIN;
      bimolecular_penalty = bimolecular_penalty - kBoltzmann * current_temp * log( waterdensity);
      setupRates( energy_options );  
      return;
    }

  // NO LONGER NEEDED:
  //    temperature = temperature + 273.15; // converted to kelvin.
  // Temerature is always stored in Kelvin and is converted to Celsius in the one spot where it needs that.
  // temperature conversions, in the order in which they show up in the class.

  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
      stack_37_dG[loop][loop2] = T_scale( stack_37_dG[loop][loop2], stack_37_dH[loop][loop2], temperature);

  for( loop = 0; loop < 31; loop++ )
    hairpin_37_dG[loop] = T_scale( hairpin_37_dG[loop], hairpin_37_dH[loop], temperature );

  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
      for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		hairpin_mismatch_37_dG[loop][loop2][loop3] = T_scale( hairpin_mismatch_37_dG[loop][loop2][loop3], hairpin_mismatch_37_dH[loop][loop2][loop3], temperature );

  for( loop = 0; loop < 120; loop++ )
    hairpin_tetraloop_37_dG[loop] = T_scale( hairpin_tetraloop_37_dG[loop], hairpin_tetraloop_37_dH[loop], temperature );


  for( loop = 0; loop < 40; loop++ )
    hairpin_triloop_37_dG[loop] = T_scale( hairpin_triloop_37_dG[loop], hairpin_triloop_37_dH[loop], temperature );

  for( loop = 0; loop < 31; loop++ )
    bulge_37_dG[loop] = T_scale( bulge_37_dG[loop], bulge_37_dH[loop], temperature );

  for( loop = 0; loop < 31; loop++ )
    internal_37_dG[loop] = T_scale( internal_37_dG[loop], internal_37_dH[loop], temperature );
  
  for( loop = 0; loop < NUM_BASES; loop++ )
    //  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
      for( loop3 = 0; loop3 < NUM_BASEPAIRS_NUPACK; loop3++ )
		//      for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		internal_mismatch_37_dG[loop][loop2][loop3] = T_scale( internal_mismatch_37_dG[loop][loop2][loop3], internal_mismatch_37_dH[loop][loop2][loop3], temperature );

  maximum_NINIO = T_scale( maximum_NINIO, maximum_NINIO_dH, temperature );

  for( loop = 0; loop < 5; loop++ )
    ninio_correction_37[loop] = T_scale( ninio_correction_37[loop], ninio_correction_37_dH[loop], temperature );


  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
      for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		for( loop4 = 0; loop4 < NUM_BASES; loop4++ )
		  internal_1_1_37_dG[loop][loop2][loop3][loop4] = T_scale( internal_1_1_37_dG[loop][loop2][loop3][loop4], internal_1_1_37_dH[loop][loop2][loop3][loop4], temperature );

  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop5 = 0; loop5 < NUM_BASES; loop5++ )
      for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
		for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		  for( loop4 = 0; loop4 < NUM_BASES; loop4++ )
			internal_2_1_37_dG[loop][loop5][loop2][loop3][loop4] = T_scale( internal_2_1_37_dG[loop][loop5][loop2][loop3][loop4], internal_2_1_37_dH[loop][loop5][loop2][loop3][loop4], temperature );

  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
	for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
	  for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		for( loop4 = 0; loop4 < NUM_BASES; loop4++ )
		  for( loop5 = 0; loop5 < NUM_BASES; loop5++ )
			for( loop6 = 0; loop6 < NUM_BASES; loop6++ )
			  internal_2_2_37_dG[loop][loop2][loop3][loop4][loop5][loop6] = T_scale( internal_2_2_37_dG[loop][loop2][loop3][loop4][loop5][loop6], internal_2_2_37_dH[loop][loop2][loop3][loop4][loop5][loop6], temperature );

  multiloop_base = T_scale( multiloop_base , multiloop_base_dH, temperature );
  multiloop_closing = T_scale( multiloop_closing , multiloop_closing_dH, temperature );
  multiloop_internal = T_scale( multiloop_internal , multiloop_internal_dH, temperature );

  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
      dangle_3_37_dG[loop][loop2] = T_scale( dangle_3_37_dG[loop][loop2], dangle_3_37_dH[loop][loop2], temperature );

  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
      dangle_5_37_dG[loop][loop2] = T_scale( dangle_5_37_dG[loop][loop2], dangle_5_37_dH[loop][loop2], temperature );

  terminal_AU = T_scale( terminal_AU, terminal_AU_dH, temperature );

  bimolecular_penalty = T_scale( bimolecular_penalty, bimolecular_penalty_dH, temperature ) - kBoltzmann * temperature * log( waterdensity );
  // need additional conversion as well

  _RT = kBoltzmann * temperature;
  current_temp = temperature;
  setupRates( energy_options );  
}



/* ------------------------------------------------------------------------
   

   Private functions for use in loading and setting data by the main constructor


   ------------------------------------------------------------------------ */


void NupackEnergyModel::internal_set_stack_energies( FILE *fp, char *buffer )
{
  int loop;
  char *cur_bufspot;
  
  while( buffer[0] == '>' )
    fgets( buffer, 2048, fp );
  //      for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
  //	{
  //	  stack_37_dG[loop][NUM_BASEPAIRS_NUPACK-1] = 0;
  //	  stack_37_dG[NUM_BASEPAIRS_NUPACK-1][loop] = 0;
  //	}
  //    }
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    {
      cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &stack_37_dG[loop][0], NUM_BASEPAIRS_NUPACK);
    }
}


void NupackEnergyModel::internal_set_stack_enthalpies( FILE *fp, char *buffer )
{
  int loop;
  char *cur_bufspot;
  
  while( buffer[0] == '>' )
    fgets( buffer, 2048, fp );

  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    {
      cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &stack_37_dH[loop][0], NUM_BASEPAIRS_NUPACK);
    }
}


void NupackEnergyModel::internal_set_hairpin_energies( FILE *fp, char *buffer )
{
  hairpin_37_dG[0] = 0;
  internal_read_array_data( fp, buffer, buffer, &hairpin_37_dG[1], 30);
}

void NupackEnergyModel::internal_set_hairpin_enthalpies( FILE *fp, char *buffer )
{
  hairpin_37_dH[0] = 0;
  internal_read_array_data( fp, buffer, buffer, &hairpin_37_dH[1], 30);
}

void NupackEnergyModel::internal_set_bulge_energies( FILE *fp, char *buffer )
{
  bulge_37_dG[0] = 0;
  internal_read_array_data( fp, buffer, buffer, &bulge_37_dG[1], 30 );
}

void NupackEnergyModel::internal_set_bulge_enthalpies( FILE *fp, char *buffer )
{
  bulge_37_dH[0] = 0;
  internal_read_array_data( fp, buffer, buffer, &bulge_37_dH[1], 30 );
}

void NupackEnergyModel::internal_set_interior_loop_energies( FILE *fp, char *buffer )
{
  internal_37_dG[0] = 0;
  internal_read_array_data( fp, buffer, buffer, &internal_37_dG[1], 30 );
}

void NupackEnergyModel::internal_set_interior_loop_enthalpies( FILE *fp, char *buffer )
{
  internal_37_dH[0] = 0;
  internal_read_array_data( fp, buffer, buffer, &internal_37_dH[1], 30 );
}

void NupackEnergyModel::internal_set_interior_1_1_energies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3;
  char *cur_bufspot;

  while( buffer[0] == '>' ) fgets( buffer, 2048, fp );
      
  cur_bufspot = buffer;
  
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
      {
		if( buffer[0] == 'A' || buffer[0] == 'T' || buffer[0] == 'G' || buffer[0] == 'C') 
		  fgets( buffer, 2048, fp); // eat the lines with AT..AT, etc, just in case.
	
		for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		  {
			internal_1_1_37_dG[loop][loop2][loop3][0] = 0;
			internal_1_1_37_dG[loop][loop2][0][loop3] = 0;
		  }
		for( loop3 = 1; loop3 < NUM_BASES; loop3++ )
		  cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_1_1_37_dG[loop][loop2][loop3][1], NUM_BASES-1);
	
      }
}

void NupackEnergyModel::internal_set_interior_1_1_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3;
  char *cur_bufspot;

  while( buffer[0] == '>' ) fgets( buffer, 2048, fp );
      
  cur_bufspot = buffer;
  
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
      {
		if( buffer[0] == 'A' || buffer[0] == 'T' || buffer[0] == 'G' || buffer[0] == 'C') 
		  fgets( buffer, 2048, fp); // eat the lines with AT..AT, etc, just in case.
	
		for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		  {
			internal_1_1_37_dG[loop][loop2][loop3][0] = 0;
			internal_1_1_37_dG[loop][loop2][0][loop3] = 0;
		  }
		for( loop3 = 1; loop3 < NUM_BASES; loop3++ )
		  cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_1_1_37_dH[loop][loop2][loop3][1], NUM_BASES-1);
	
      }
}

void NupackEnergyModel::internal_set_interior_2_1_energies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3,loop4;
  char *cur_bufspot;

  while( buffer[0] == '>' ) fgets( buffer, 2048, fp );  

  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    {
      for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
		{	
		  if( buffer[0] == 'A' || buffer[0] == 'T' || buffer[0] == 'G' || buffer[0] == 'C')
			fgets( buffer, 2048, fp); // eat the lines with AT..AT, etc, just in case. 
	  
		  for( loop3 = 1; loop3 < NUM_BASES; loop3++ )
			{
			  for( loop4 = 1; loop4 < NUM_BASES; loop4++ )
				cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_2_1_37_dG[loop][loop3][loop2][loop4][1], (NUM_BASES-1));
			}
		}
    }
}

void NupackEnergyModel::internal_set_interior_2_1_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3,loop4;
  char *cur_bufspot;

  while( buffer[0] == '>' ) fgets( buffer, 2048, fp );  

  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    {
      for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
		{	
		  if( buffer[0] == 'A' || buffer[0] == 'T' || buffer[0] == 'G' || buffer[0] == 'C')
			fgets( buffer, 2048, fp); // eat the lines with AT..AT, etc, just in case. 
	  
		  for( loop3 = 1; loop3 < NUM_BASES; loop3++ )
			{
			  for( loop4 = 1; loop4 < NUM_BASES; loop4++ )
				cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_2_1_37_dH[loop][loop3][loop2][loop4][1], (NUM_BASES-1));
			}
		}
    }
}

void NupackEnergyModel::internal_set_interior_2_2_energies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3,loop4,loop5;
  char *cur_bufspot;
  while( buffer[0] == '>')
	fgets( buffer, 2048, fp );


  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
	{
	  for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
	    {
	      if( buffer[0] == 'A' || buffer[0] == 'T' || buffer[0] == 'G' || buffer[0] == 'C') 
			fgets( buffer, 2048, fp ); // eat the description lines.
	      
	      //cur_bufspot=buffer;
	      for( loop3 = 1; loop3 < NUM_BASES; loop3++ )
			{
			  for( loop4 = 1; loop4 < NUM_BASES; loop4++ )
				{
				  for( loop5 = 1; loop5 < NUM_BASES; loop5++ )
					{
					  cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_2_2_37_dG[loop][loop2][loop3][loop4][loop5][1], NUM_BASES-1);
					}
				}
			}
	    }
	} 
}

void NupackEnergyModel::internal_set_interior_2_2_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3,loop4,loop5;
  char *cur_bufspot;
  while( buffer[0] == '>')
	fgets( buffer, 2048, fp );


  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
	{
	  for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
	    {
	      if( buffer[0] == 'A' || buffer[0] == 'T' || buffer[0] == 'G' || buffer[0] == 'C') 
			fgets( buffer, 2048, fp ); // eat the description lines.
	      
	      //cur_bufspot=buffer;
	      for( loop3 = 1; loop3 < NUM_BASES; loop3++ )
			{
			  for( loop4 = 1; loop4 < NUM_BASES; loop4++ )
				{
				  for( loop5 = 1; loop5 < NUM_BASES; loop5++ )
					{
					  cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_2_2_37_dH[loop][loop2][loop3][loop4][loop5][1], NUM_BASES-1);
					}
				}
			}
	    }
	} 
}

void NupackEnergyModel::internal_set_dangle_5_energies( FILE *fp, char *buffer )
{
  char *cur_bufspot;
  while( buffer[0] == '>' ) fgets(buffer, 2048, fp);

  int loop,loop2;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++)
      dangle_5_37_dG[loop][loop2] = 0.0;

  cur_bufspot=buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++)
    cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &dangle_5_37_dG[loop][1], NUM_BASES -1);      
}

void NupackEnergyModel::internal_set_dangle_5_enthalpies( FILE *fp, char *buffer )
{
  char *cur_bufspot;
  while( buffer[0] == '>' ) fgets(buffer, 2048, fp);

  int loop,loop2;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++)
      dangle_5_37_dH[loop][loop2] = 0;
  
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++)
    cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &dangle_5_37_dH[loop][1], NUM_BASES -1);      
}


void NupackEnergyModel::internal_set_dangle_3_energies( FILE *fp, char *buffer )
{
  char *cur_bufspot;
  while( buffer[0] == '>' ) fgets(buffer, 2048, fp);
  
  int loop,loop2;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++)
      dangle_3_37_dG[loop][loop2] = 0.0;
  
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++)
    cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &dangle_3_37_dG[loop][1], NUM_BASES -1);      
}

void NupackEnergyModel::internal_set_dangle_3_enthalpies( FILE *fp, char *buffer )
{
  char *cur_bufspot;
  while( buffer[0] == '>' ) fgets(buffer, 2048, fp);
  
  int loop,loop2;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++)
      dangle_3_37_dH[loop][loop2] = 0;
  
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++)
    cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &dangle_3_37_dH[loop][1], NUM_BASES -1);      
}

void NupackEnergyModel::internal_set_multiloop_parameters( FILE *fp, char *buffer )
{
  double temp[4];
  while( buffer[0] == '>') fgets(buffer, 2048, fp);

  internal_read_array_data( fp, buffer, buffer, temp, 3 );
  multiloop_base = temp[2];
  multiloop_closing = temp[0];
  multiloop_internal = temp[1];
}

void NupackEnergyModel::internal_set_multiloop_parameters_enthalpies( FILE *fp, char *buffer )
{
  int temp[4];
 
  while( buffer[0] == '>') fgets(buffer, 2048, fp);
  
  internal_read_array_data( fp, buffer, buffer, temp, 3 );
  multiloop_base_dH = temp[2];
  multiloop_closing_dH = temp[0];
  multiloop_internal_dH = temp[1];
}

void NupackEnergyModel::internal_set_at_penalty( FILE *fp, char *buffer )
{
  while( buffer[0] == '>' ) fgets( buffer, 2048, fp);
      
  internal_read_array_data( fp, buffer, buffer, &terminal_AU, 1 );
}

void NupackEnergyModel::internal_set_at_penalty_enthalpy( FILE *fp, char *buffer )
{
  while( buffer[0] == '>' ) fgets( buffer, 2048, fp);
      
  internal_read_array_data( fp, buffer, buffer, &terminal_AU_dH, 1 );
}


void NupackEnergyModel::internal_set_bimolecular_penalty( FILE *fp, char *buffer )
{
  // Vienna parameter set doesn't have this term.
  while( buffer[0] == '>' ) fgets( buffer, 2048, fp);

  internal_read_array_data( fp, buffer, buffer, &bimolecular_penalty, 1 );
}

void NupackEnergyModel::internal_set_bimolecular_penalty_dH( FILE *fp, char *buffer )
{
  while( buffer[0] == '>' ) fgets( buffer, 2048, fp);

  internal_read_array_data( fp, buffer, buffer, &bimolecular_penalty_dH, 1 );
}

void NupackEnergyModel::internal_set_ninio_parameters( FILE *fp, char *buffer )
{
  double temp[5];
  while( buffer[0] == '>')
    fgets(buffer, 2048, fp);
  internal_read_array_data( fp, buffer, buffer, temp, 5 );
  maximum_NINIO = temp[4];
  ninio_correction_37[0] = temp[0];
  ninio_correction_37[1] = temp[1];
  ninio_correction_37[2] = temp[2];
  ninio_correction_37[3] = temp[3];
}

void NupackEnergyModel::internal_set_ninio_parameters_enthalpy( FILE *fp, char *buffer )
{
  int temp[5];
  while( buffer[0] == '>')
    fgets(buffer, 2048, fp);
  internal_read_array_data( fp, buffer, buffer, temp, 5 );
  maximum_NINIO_dH = temp[4];
  ninio_correction_37_dH[0] = temp[0];
  ninio_correction_37_dH[1] = temp[1];
  ninio_correction_37_dH[2] = temp[2];
  ninio_correction_37_dH[3] = temp[3];
}


void NupackEnergyModel::internal_set_hairpin_tetraloop_parameters( FILE *fp, char *buffer )
{
  int tetra_index = 0;
  int buf_index = 0;
  fgets( buffer, 2048, fp );
  while( strlen(buffer) > 7  && buffer[0] != '>')
    {
      buf_index = 0;
      while( isspace(buffer[buf_index])) buf_index++;
      
      for( int loop = 0; loop < 6; loop ++ )
		{
		  if(buffer[buf_index+loop] == 'U')
			hairpin_tetraloops[ (tetra_index * 7) + loop] = baseLookup('T');
		  else
			hairpin_tetraloops[ (tetra_index * 7) + loop] = baseLookup(buffer[buf_index + loop]);
		}
      hairpin_tetraloops[ (tetra_index * 7) + 6] = ' ';
      hairpin_tetraloops[ (tetra_index * 7) + 7] = 0;
      hairpin_tetraloop_37_dG[tetra_index] = atof( &buffer[buf_index + 6]) / 100.0;

      tetra_index++;
      fgets( buffer, 2048, fp );
    }
#ifdef DEBUG
  fprintf(stderr,"Tetraloop Paramaters (MFOLD): %d read.\n",tetra_index);
#endif
}

void NupackEnergyModel::internal_set_hairpin_tetraloop_parameters_enthalpy( FILE *fp, char *buffer )
{
  char *index;
  int i_index;
  int value,loop;
  for( loop = 0; loop < 120; loop++ )
    hairpin_tetraloop_37_dH[loop]=0;
  fgets( buffer, 2048, fp );
  while( strlen(buffer) > 7 && buffer[0] != '>' )
    {
      value = atoi(buffer+6);
      buffer[6] = '\0';
      for( loop = 0; loop < 6; loop ++ )
		{
		  if(buffer[loop] == 'U')
			buffer[loop] = baseLookup('T');
		  else
			buffer[loop] = baseLookup(buffer[loop]);
		}     
      index = strstr(hairpin_tetraloops,buffer);
      if( index != NULL )
		{
		  i_index = (int) (index - hairpin_tetraloops) / 7;

		  hairpin_tetraloop_37_dH[i_index] = value;
		}
      fgets(buffer, 2048, fp );
    }
}

void NupackEnergyModel::internal_set_hairpin_triloop_parameters( FILE *fp, char *buffer )
{
  int tri_index = 0;
  int buf_index = 0;
  fgets( buffer, 2048, fp );
  while( strlen(buffer) > 6 && buffer[0] != '>')
    {
      buf_index = 0;
      while( isspace(buffer[buf_index])) buf_index++;
      
      for( int loop = 0; loop < 5; loop ++ )
		{
		  if(buffer[buf_index+loop] == 'U')
			hairpin_triloops[ (tri_index * 6) + loop] = baseLookup('T');
		  else
			hairpin_triloops[ (tri_index * 6) + loop] = baseLookup(buffer[buf_index + loop]);
		}

      hairpin_triloops[ (tri_index * 6) + 5] = ' ';
      hairpin_triloops[ (tri_index * 6) + 6] = 0;
      hairpin_triloop_37_dG[tri_index] = atof( &buffer[buf_index + 5]) / 100.0;

      tri_index++;
      fgets( buffer, 2048, fp );
    }
#ifdef DEBUG
  fprintf(stderr,"Triloop Paramaters (MFOLD): %d read.\n",tri_index);
#endif
}

void NupackEnergyModel::internal_set_hairpin_triloop_parameters_enthalpy( FILE *fp, char *buffer )
{
  char *index;
  int i_index;
  int value,loop;  
  for( loop = 0; loop < 40; loop++ )
    hairpin_triloop_37_dH[loop]=0;

  fgets( buffer, 2048, fp );
  while( strlen(buffer) > 6 && buffer[0] != '>' )
    {
      value = atoi(buffer+5);
      buffer[5] = '\0';
      for( loop = 0; loop < 5; loop ++ )
		{
		  if(buffer[loop] == 'U')
			buffer[loop] = baseLookup('T');
		  else
			buffer[loop] = baseLookup(buffer[loop]);
		}
      
      index = strstr(hairpin_triloops,buffer);
      if( index != NULL )
		{
		  i_index = (int) (index - hairpin_triloops) / 6;

		  hairpin_triloop_37_dH[i_index] = value;
		}
      fgets(buffer, 2048, fp );
    }

}

void NupackEnergyModel::internal_set_hairpin_mismatch_energies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3;
  char *cur_bufspot;
  
  double temp[NUM_BASEPAIRS_NUPACK];
  cur_bufspot = buffer;
  while( buffer[0] == '>' ) fgets(buffer, 2048, fp);
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
      for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		hairpin_mismatch_37_dG[loop][loop2][loop3] = 0;
  
  for( loop = 0; loop < (NUM_BASES-1) * (NUM_BASES-1); loop++)
    {
      cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &temp[0], NUM_BASEPAIRS_NUPACK);
      loop3 = (loop-(loop % (NUM_BASES-1))) / (NUM_BASES-1);
      for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
		hairpin_mismatch_37_dG[loop2][loop3+1][(loop % (NUM_BASES-1))+1] = temp[loop2];
      
    }
}


void NupackEnergyModel::internal_set_hairpin_mismatch_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3;
  char *cur_bufspot;
  
  
  int temp[NUM_BASEPAIRS_NUPACK];
  cur_bufspot = buffer;
  while( buffer[0] == '>' ) fgets(buffer, 2048, fp);
  for( loop = 0; loop < NUM_BASEPAIRS_NUPACK; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
      for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
		hairpin_mismatch_37_dH[loop][loop2][loop3] = 0;
  
  for( loop = 0; loop < (NUM_BASES-1) * (NUM_BASES-1); loop++)
    {
      cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &temp[0], NUM_BASEPAIRS_NUPACK);
      loop3 = (loop-(loop % (NUM_BASES-1))) / (NUM_BASES-1);
      for( loop2 = 0; loop2 < NUM_BASEPAIRS_NUPACK; loop2++ )
		hairpin_mismatch_37_dH[loop2][loop3+1][(loop % (NUM_BASES-1))+1] = temp[loop2];
      
    }
}

void NupackEnergyModel::internal_set_interior_loop_mismatch_energies( FILE *fp, char *buffer )
{
  int loop, loop2, loop3;
  char *cur_bufspot;

  while( buffer[0] == '>' ) fgets(buffer, 2048, fp);
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASES; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
      for( loop3 = 0; loop3 < NUM_BASEPAIRS_NUPACK; loop3++ )
		internal_mismatch_37_dG[loop][loop2][loop3] = 0;

  for( loop = 1; loop < NUM_BASES; loop++)
    for( loop2 = 1; loop2 < NUM_BASES; loop2++ )
	  {
		cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_mismatch_37_dG[loop][loop2][0], NUM_BASEPAIRS_NUPACK);
	  }

}

void NupackEnergyModel::internal_set_interior_loop_mismatch_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2, loop3;
  char *cur_bufspot;

  while( buffer[0] == '>' ) fgets(buffer, 2048, fp);
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASES; loop++ )
    for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
      for( loop3 = 0; loop3 < NUM_BASEPAIRS_NUPACK; loop3++ )
		internal_mismatch_37_dH[loop][loop2][loop3] = 0;

  for( loop = 1; loop < NUM_BASES; loop++)
    for( loop2 = 1; loop2 < NUM_BASES; loop2++ )
	  {
		cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_mismatch_37_dH[loop][loop2][0], NUM_BASEPAIRS_NUPACK);
	  }
}

char *NupackEnergyModel::internal_read_array_data( FILE *fp, char *buffer, char *start_loc, int *read_loc, int size )
{
  int loop,loop2;
  char *cur_bufspot, *temp_char = NULL;
  int temp_int;

  cur_bufspot = start_loc;
  if( cur_bufspot[0] == '>' )
    {
      fgets( buffer, 2048, fp );
      cur_bufspot = buffer;
    }
  for( loop = 0; loop < size; loop++ )
    {
      temp_int = strtol(cur_bufspot, &temp_char, 10 );
      if( cur_bufspot == temp_char ) // we didn't read any value
		{
		  if( strstr( cur_bufspot, "/*") != NULL )
			{
			  cur_bufspot = strstr( cur_bufspot, "*/")+2;
			  loop--;
			}
		  // check to see if the value is INF
		  else if( (temp_char = strstr( cur_bufspot, "INF")) != NULL )
			{
			  // and if it is, set it correctly
			  read_loc[loop] = INF;
			  cur_bufspot = temp_char + 3;
			}
		  else if( (temp_char = strstr( cur_bufspot, "x")) != NULL )
			{
			  if( loop == 0 ) // we have an error. ERROR
				;
			  else
				read_loc[loop] = read_loc[loop-1] + (int) rint( log_loop_penalty * log(((double) loop) / ((double) (loop - 1)))) ;
			  cur_bufspot = temp_char +1;
			}
		  else // we need to check for leading characters
			{
			  for( loop2 = 0; loop2 < strlen(cur_bufspot); loop2++)
				if( isdigit( cur_bufspot[loop2] ) )
				  {
					break;
				  }
	      
			  if( loop2 != strlen(cur_bufspot))
				cur_bufspot = cur_bufspot + loop2 - 1;
			  else
				{
				  // otherwise get more data from the stream and reset counters so we reread to fill the same value.
				  fgets( buffer, 2048, fp );
				  cur_bufspot = buffer;
				}
			  loop--;
			}
		}
      else
		{
		  // we got a value to fill our current spot, so increment
		  cur_bufspot = temp_char;
		  // and fill the value into our array
		  read_loc[loop] = temp_int;
		}
    }

  return cur_bufspot;
}


char *NupackEnergyModel::internal_read_array_data( FILE *fp, char *buffer, char *start_loc, double *read_loc, int size )
{
  int loop,loop2;
  char *cur_bufspot, *temp_char;
  double temp_double;

  cur_bufspot = start_loc;
  if( cur_bufspot[0] == '>' )
    {
      fgets( buffer, 2048, fp );
      cur_bufspot = buffer;
    }
  for( loop = 0; loop < size; loop++ )
    {
      temp_double = strtod(cur_bufspot, &temp_char) / 100.0;
      if( cur_bufspot == temp_char ) // we didn't read any value
		{
		  if( strstr( cur_bufspot, "/*") != NULL )
			{
			  cur_bufspot = strstr( cur_bufspot, "*/")+2;
			  loop--;
			}
		  // check to see if the value is INF
		  else if( (temp_char = strstr( cur_bufspot, "INF")) != NULL )
			{
			  // and if it is, set it correctly
			  read_loc[loop] = INF;
			  cur_bufspot = temp_char + 3;
			}
		  else if( (temp_char = strstr( cur_bufspot, "x")) != NULL )
			{
			  if( loop == 0 ) // we have an error. ERROR
				;
			  else
				read_loc[loop] = read_loc[loop-1] + (double) rint( log_loop_penalty * log(((double) loop) / ((double) (loop - 1)))) ;
			  cur_bufspot = temp_char +1;
			}
		  else // we need to check for leading characters
			{
			  for( loop2 = 0; loop2 < strlen(cur_bufspot); loop2++)
				if( isdigit( cur_bufspot[loop2] ) )
				  {
					break;
				  }
	      
			  if( loop2 != strlen(cur_bufspot))
				cur_bufspot = cur_bufspot + loop2 - 1;
			  else
				{
				  // otherwise get more data from the stream and reset counters so we reread to fill the same value.
				  fgets( buffer, 2048, fp );
				  cur_bufspot = buffer;
				}
			  loop--;
			}
		}
      else
		{
		  // we got a value to fill our current spot, so increment
		  cur_bufspot = temp_char;
		  // and fill the value into our array
		  read_loc[loop] = temp_double;
		}
    }

  return cur_bufspot;
}

void NupackEnergyModel::setupRates( PyObject *energy_options )
{
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


  double joinconc,joinrate_volume;
  getDoubleAttr(energy_options, join_concentration,&joinconc);
  // the concentration is now in M units, rather than the previous mM. The input parser converts to these.

  getDoubleAttr(energy_options, bimolecular_scaling,&biscale);
  getDoubleAttr(energy_options, unimolecular_scaling,&uniscale);
  // Two components to the join rate, the dG_assoc from mass action, and the dG_volume term which is related to the concentration. We compute each individually, following my derivation for the volume term, and the method used in Nupack for the dG_assoc term.

  // Please note, waterdensity is actually the Mols H2O / L, not the actual density of water (g/mol).
  dG_volume = _RT * log( waterdensity / joinconc);
  joinrate_volume = joinconc / waterdensity;
  //              = exp( -dG_volume / _RT )

  // dG_volume = _RT * log( waterdensity / (.001 * joinconc));
  // joinrate_volume = .001 * joinconc / waterdensity;
  // old, for units of mM


  
  dG_assoc = bimolecular_penalty; // already computed and scaled for water density. 

  // joinenergy = dG_volume + dG_assoc;
  joinrate = biscale * joinrate_volume; // * exp( - dG_assoc / _RT );

}

double NupackEnergyModel::setWaterDensity( double temp )
{
  // returns the density of water at a given temperature (currently passed in degrees C)
  // Need to use data via appropriate table - nupack references:
  // "Recommended table for the density of water between 0C and 40C based on recent experimental reports", Tanaka M., Girard, G. et al, Metrologia, 2001 38, 301-309. 
  // Look this up and enter the full formula from there. The paper is not available for free, probably need to use a caltech login or order it.

  // uses formulation via Nupack:
  double values[6] = { -3.983035, 301.797, 522528.9, 69.34881, 999.97495, 18.0152 };
  
  values[0] = temp + values[0];
  values[1] = temp + values[1];
  values[3] = temp + values[3];
  return values[4] * (1 - values[0] * values[0] * values[1] / values[2] / values[3]) / values[5];

  //  return 55.6; // default for now, in units of mols/liter.
}





