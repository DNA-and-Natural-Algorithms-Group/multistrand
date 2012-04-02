/*
  Copyright (c) 2007-2008 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#include "options.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <ctype.h>
#include "energymodel.h"

#undef DEBUG
//#define DEBUG

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


ViennaEnergyModel::~ViennaEnergyModel( void )
{
  // TODO: is anything allocated now? Don't think so, all arrays are static still.
  // nothing is allocated within an energy model.
}


double ViennaEnergyModel::returnRate( double start_energy, double end_energy, int enth_entr_toggle)
{
  double dE = ( end_energy - start_energy );
  //  double dE = ((double) ( end_energy - start_energy )) / 100.;
  if( kinetic_rate_method == 2 )
    return exp( -0.5 * dE / _RT );
}

double ViennaEnergyModel::getVolumeEnergy( void )
{
  return 0.0;
}

double ViennaEnergyModel::getAssocEnergy( void )
{
  return 0.0;
}

double ViennaEnergyModel::getJoinRate( void )
{
  //  printf("%lf\n",joinrate);
  return joinrate; // replace with the passed in rate
}

double ViennaEnergyModel::getJoinRate_NoVolumeTerm( void )
{
  // TODO: update Vienna model rates.
  return joinrate;
}

// entropy/enthalpy energy parameters

void ViennaEnergyModel::eStackEnergy( int type1, int type2, energyS *energy )
{
  energy->dH = stack_37_dH[type1][basepair_sw[type2]];
  energy->nTdS = stack_37_dG[type1][basepair_sw[type2]] - energy->dH;
  //  return stack_37_dG[type1][basepair_sw[type2]];
}

// regular parameters

double ViennaEnergyModel::StackEnergy( int i, int j, int p, int q)
{
  return (double) stack_37_dG[pairtypes[i][j]][pairtypes[q][p]] / 100.;
}


double ViennaEnergyModel::BulgeEnergy( int i, int j, int p, int q, int bulgesize )
{
  int energy = 0;
  if( bulgesize <= 30 )
    energy = bulge_37_dG[bulgesize];
  else
    energy = bulge_37_dG[30] + (int) (log((double)bulgesize / 30.0) * log_loop_penalty_37);

  if( bulgesize == 1) // add stacking term for single-base bulges.
    energy += stack_37_dG[pairtypes[i][j]][pairtypes[q][p]];
  else // AU penalty doesn't apply if they stack.
    {
      if( pairtypes[i][j] > 2 ) // AU penalty applies
        energy += terminal_AU;
      if( pairtypes[q][p] > 2 ) // AU penalty applies
        energy += terminal_AU;
    }

  return (double) energy / 100.;
}

double ViennaEnergyModel::InteriorEnergy( char *seq1, char *seq2, int size1, int size2 )
{
  int energy,ninio;

  int type1 = pairtypes[seq1[0]][seq2[size2+1]];
  int type2 = pairtypes[seq2[0]][seq1[size1+1]];

  // mismatch0 = seq1[1]
  // mismatch1 = seq2[size2]
  // mismatch2 = seq1[size1]
  // mismatch3 = seq2[1]

  // special case time. 1x1, 2x1 and 2x2's all get special cases.
  if( size1 == 1 && size2 == 1)
    return (double)internal_1_1_37_dG[type1][type2][seq1[1]][seq2[size2]] / 100.0;
  if( size1 <= 2 && size2 <= 2)
    if( size1 == 1 || size2 == 1)
      {
        if( size1 == 1)
          return (double) internal_2_1_37_dG[type1][type2][seq1[1]][seq2[1]][seq2[size2]] / 100.0;
        else
          return (double) internal_2_1_37_dG[type2][type1][seq2[1]][seq1[1]][seq1[size1]] / 100.0;
      }
  if( size1 == 2 && size2 == 2)
    return (double) internal_2_2_37_dG[type1][type2][seq1[1]][seq1[size1]][seq2[1]][seq2[size2]] / 100.0;


  // Generic case.

  if( size1 + size2 <= 30 )
    energy = internal_37_dG[ size1 + size2 ];
  else
    energy = internal_37_dG[30] + (int) (log((double)(size1+size2) / 30.0) * log_loop_penalty_37);

  // NINIO term... no idea why this is used.
  ninio = abs(size2 - size1) * ninio_correction_37[2]; // don't ask me why, this is according to the Vienna energy model and this term is in it. 
  
  if( maximum_NINIO < ninio )
    energy += maximum_NINIO;
  else
    energy += ninio;

  energy += internal_mismatch_37_dG[type1][seq1[1]][seq2[size2]] + 
    internal_mismatch_37_dG[type2][seq2[1]][seq1[size1]];

  return (double) energy / 100.0;

}

double ViennaEnergyModel::HairpinEnergy( char *seq, int size )
{
  int energy = 0;
  int lookup_index = 0;

  if( size <= 30 )
    energy = hairpin_37_dG[size];
  else
    energy = hairpin_37_dG[30] + (int) (log((double)size / 30.0) * log_loop_penalty_37);
  
  if( size == 3 ) // triloop bonuses
    {
      // We now always do the lookup, if no entry then it is 0.0.
      lookup_index = ((seq[0]-1) << 8) + 
        ((seq[1]-1) << 6) + 
        ((seq[2]-1) << 4) + 
        ((seq[3]-1) << 2) +
        (seq[4]-1);

	  energy += hairpin_triloop_37_dG[lookup_index];

      if( (seq[0] == BASE_T) || (seq[size+1] == BASE_T) )
        energy += terminal_AU;

    }
  if( size == 4)
    {
      lookup_index = ((seq[0]-1) << 10) + 
        ((seq[1]-1) << 8) + 
        ((seq[2]-1) << 6) + 
        ((seq[3]-1) << 4) +
        ((seq[4]-1) << 2) +
        (seq[5]-1);

      energy += hairpin_tetraloop_37_dG[lookup_index];
    }
  // in the below, note that (pairtypes[seq[0]][seq[size+1]]-1) is the base pair type
  // forming the hairpin, and the terminal_AU energy gets included via the mismatch calc.
  if( size >= 4 )
    energy += hairpin_mismatch_37_dG[pairtypes[seq[0]][seq[size+1]]][seq[1]][seq[size]];

  return (double) energy / 100.0;
}

double ViennaEnergyModel::MultiloopEnergy( int size, int *sidelen, char **sequences)
{
  // no dangle terms yet, this is equiv to dangles = 0;
  int totallength=0,energy=0;
  int pt;
  int loopminus1 = size-1;

  // pairtypes[loop] = pairtypes[ ][ ]
  // pairtypes[loop] is the last of sequence[loop-1] and first of sequence[loop]
  // pairtypes[loop] = pairtypes[sequences[loop-1][sidelen[loop-1]+1]][sequences[loop][0]];
  // pairtypes[sequences[loop-1][sidelen[loop-1]+1]][sequences[loop][0]];

  for( int loop = 0; loop < size; loop++)
    {
      totallength += sidelen[loop];
      pt = pairtypes[sequences[loopminus1][sidelen[loopminus1]+1]][sequences[loop][0]];

      if( pt > 2) // AT penalty applies
        energy += terminal_AU;
      loopminus1++;
      if( loopminus1 == size ) loopminus1 = 0;
    }
  energy += size * multiloop_internal;
  energy += multiloop_closing;
  if( totallength <= 6 )
    energy += multiloop_base * totallength;
  else
    energy += multiloop_base * 6 + (int) (log((double)totallength / 6.0) * log_loop_penalty_37);

  return (double) energy / 100.0;
}


double ViennaEnergyModel::OpenloopEnergy( int size, int *sidelen, char **sequences)
{
  // no dangle terms yet, this is equiv to dangles = 0;
  int energy=0;
  int pt;

  for( int loop = 0; loop < size; loop++)
    {
      pt = pairtypes[sequences[loop][sidelen[loop]+1]][sequences[loop+1][0]];
      if( pt > 2)  // AT penalty applies
        energy += terminal_AU;
    }
  return (double) energy / 100.0;
}


// constructors, internal functions


ViennaEnergyModel::ViennaEnergyModel( PyObject *energy_options ) : log_loop_penalty_37(107.856) , kinetic_rate_method(2) , _RT(.6) , bimolecular_penalty(196) // Check references for this loop penalty term.
{
  // This is the tough part, performing all read/input duties.
  char in_buffer[2048];
  int loop, loop2 ;
  FILE *fp = NULL;


  if( testLongAttr(energy_options, parameter_type, =, ENERGYMODEL_VIENNA) && testLongAttr(energy_options, substrate_type, =, SUBSTRATE_DNA) )
    {
      fp = fopen( "dna.par", "rt");
      if( fp == NULL )
        {
          fprintf(stderr,"ERROR: Could not find Vienna DNA parameter file \"dna.par\" in the path.\n" );
          exit(1);
        }
    }
  else if( testLongAttr(energy_options, parameter_type, =, ENERGYMODEL_VIENNA) && testLongAttr(energy_options, substrate_type, =, SUBSTRATE_RNA) )
    {
      fp = fopen( "default.par", "rt");
      if( fp == NULL )
        {
          fprintf(stderr,"ERROR: Could not find Vienna RNA parameter file \"default.par\" in the path.\n" );
          exit(1);
        }
    }
      
  getDoubleAttr(energy_options, join_concentration, &joinrate);
  joinrate = .001 * joinrate / 55.6;
  getLongAttr(energy_options, dangles, &dangles);
  
  getLongAttr(energy_options, parameter_type, &ptype);
  if( ptype == VIENNA ) // dna.par (Vienna)
    {
      for( loop = 0; loop < NUM_BASES; loop++ )
        pairs[loop] = pairs_vienna[loop];
      for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        basepair_sw[loop] = basepair_sw_vienna[loop];
      for( loop = 0; loop < NUM_BASES; loop++)
        for( loop2 = 0; loop2 < NUM_BASES; loop2++)
          pairtypes[loop][loop2] = pairtypes_vienna[loop][loop2];
    }

  fgets( in_buffer, 2048, fp );
  while( !feof(fp) )
    {
      ////////////////////////////////////////////////////////
      //                 VIENNA param inputs                //
      ////////////////////////////////////////////////////////
	  if( in_buffer[0] == '#' ) // signal for a data area (vienna)
	    {
	      
	      // Switch to each member function for handling each data type.
	      if( strncmp( in_buffer, "# stack_energies", 16) == 0)
            {
#ifdef DEBUG
              printf("Loading Stack Energies.\n");
#endif
              internal_set_stack_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# stack_enthalpies", 18) == 0)
            {
#ifdef DEBUG
              printf("Loading Stack Enthalpies.\n");
#endif
              internal_set_stack_enthalpies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# hairpin", 9) == 0)
            {
#ifdef DEBUG
              printf("Loading Hairpin Energies.\n");
#endif
              internal_set_hairpin_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# mismatch_hairpin", 18) == 0)
            {
#ifdef DEBUG
              printf("Loading Hairpin Mismatch Energies.\n");
#endif
              internal_set_hairpin_mismatch_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# mismatch_interior", 19) == 0)
            {
#ifdef DEBUG
              printf("Loading Interior Loop Mismatch Energies.\n");
#endif
              internal_set_interior_loop_mismatch_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# mismatch_multi", 16) == 0)
            {
#ifdef DEBUG
              printf("Loading Multiloop Mismatch Energies.\n");
#endif
              internal_set_multiloop_mismatch_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# mismatch_enthalpies", 18) == 0)
            {
#ifdef DEBUG
              printf("Loading Hairpin Mismatch Energies.\n");
#endif
              internal_set_mismatch_enthalpies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# bulge", 7) == 0)
            {
#ifdef DEBUG
              printf("Loading Bulge Energies.\n");
#endif
              internal_set_bulge_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# internal_loop", 15) == 0)
            {
#ifdef DEBUG
              printf("Loading Internal Loop Energies.\n");
#endif
              internal_set_interior_loop_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# int11_energies", 16) == 0)
            {
#ifdef DEBUG
              printf("Loading Internal 1-1 mismatch Energies.\n");
#endif
              internal_set_interior_1_1_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# int11_enthalpies", 15) == 0)
            {
#ifdef DEBUG
              printf("Loading Internal 1-1 mismatch Enthalpies.\n");
#endif
              internal_set_interior_1_1_enthalpies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# int21_energies", 16) == 0)
            {
#ifdef DEBUG
              printf("Loading Internal 2-1 mismatch Energies.\n");
#endif
              internal_set_interior_2_1_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# int21_enthalpies", 15) == 0)
            {
#ifdef DEBUG
              printf("Loading Internal 2-1 mismatch Enthalpies.\n");
#endif
              internal_set_interior_2_1_enthalpies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# int22_energies", 16) == 0)
            {
#ifdef DEBUG
              printf("Loading Internal 2-2 mismatch Energies.\n");
#endif
              internal_set_interior_2_2_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# int22_enthalpies", 15) == 0)
            {
#ifdef DEBUG
              printf("Loading Internal 2-2 mismatch Enthalpies.\n");
#endif
              internal_set_interior_2_2_enthalpies( fp, in_buffer );
            }
	  
	      // The following two must be BEFORE the dangle5 and dangle3 checks.
	      // Otherwise, it will incorrectly recognize them (the same first 9 chars 
          // are in both.
	      if( strncmp( in_buffer, "# dangle5_enthalpies", 19) == 0)
            {
#ifdef DEBUG
              printf("Loading Dangle 5' Enthalpies.\n");
#endif
              internal_set_dangle_5_enthalpies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# dangle3_enthalpies", 19) == 0)
            {
#ifdef DEBUG
              printf("Loading Dangle 3' Enthalpies.\n");
#endif
              internal_set_dangle_3_enthalpies( fp, in_buffer );
            }
	      
	      
	      if( strncmp( in_buffer, "# dangle5", 9) == 0)
            {
#ifdef DEBUG
              printf("Loading Dangle 5' Energies.\n");
#endif
              internal_set_dangle_5_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# dangle3", 9) == 0)
            {
#ifdef DEBUG
              printf("Loading Dangle 3' Energies.\n");
#endif
              internal_set_dangle_3_energies( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# ML_params", 11) == 0)
            {
#ifdef DEBUG
              printf("Loading Multiloop parameters.\n");
#endif
              internal_set_multiloop_parameters( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# NINIO", 7) == 0)
            {
#ifdef DEBUG
              printf("Loading NINIO parameters.\n");
#endif
              internal_set_ninio_parameters( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# Tetraloops", 12) == 0)
            {
#ifdef DEBUG
              printf("Loading Hairpin Tetraloop parameters.\n");
#endif
              internal_set_hairpin_tetraloop_parameters( fp, in_buffer );
            }
	      if( strncmp( in_buffer, "# Triloops", 10) == 0)
            {
#ifdef DEBUG
              printf("Loading Hairpin Triloop parameters.\n");
#endif
              internal_set_hairpin_triloop_parameters( fp, in_buffer );
            }
	    }
	  
	  fgets( in_buffer, 2048, fp );
	  if( feof(fp) )
	    continue;	    
	}
}



/* ------------------------------------------------------------------------
   

   Private functions for use in loading and setting data by the main constructor


   ------------------------------------------------------------------------ */


void ViennaEnergyModel::internal_set_stack_energies( FILE *fp, char *buffer )
{
  int loop, loop2;
  char *cur_bufspot;
  
  if( ptype == MFOLD )
    {
      while( buffer[0] == '>' )
        fgets( buffer, 2048, fp );
      for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        {
          stack_37_dG[loop][NUM_BASEPAIRS_VIENNA-1] = 0;
          stack_37_dG[NUM_BASEPAIRS_VIENNA-1][loop] = 0;
        }
    }
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
    {
      if( loop == 0 )
        for( loop2 = 0; loop2 < NUM_BASEPAIRS_VIENNA; loop2++ )
          stack_37_dG[loop][loop2] = INF;
      else
        {
          stack_37_dG[loop][0] = INF;
          if( ptype == VIENNA )
            cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &stack_37_dG[loop][1], NUM_BASEPAIRS_VIENNA - 1);
          else if( ptype == MFOLD )
            cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &stack_37_dG[loop][1], NUM_BASEPAIRS_VIENNA - 2);
        }
      if( ptype == MFOLD && loop == NUM_BASEPAIRS_VIENNA-2) loop = NUM_BASEPAIRS_VIENNA;
    }
}


void ViennaEnergyModel::internal_set_stack_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2;
  char *cur_bufspot;
  
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
    {
      if( loop == 0 )
        for( loop2 = 0; loop2 < NUM_BASEPAIRS_VIENNA; loop2++ )
          stack_37_dH[loop][loop2] = INF;
      else
        {
          stack_37_dH[loop][0] = INF;
          cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &stack_37_dH[loop][1], NUM_BASEPAIRS_VIENNA - 1);
        }
    }
}




void ViennaEnergyModel::internal_set_hairpin_energies( FILE *fp, char *buffer )
{
  if( ptype == VIENNA)
    internal_read_array_data( fp, buffer, buffer, hairpin_37_dG, 31 );
  if( ptype == MFOLD )
    {
      hairpin_37_dG[0] = 0;
      internal_read_array_data( fp, buffer, buffer, &hairpin_37_dG[1], 30);
    }
}

void ViennaEnergyModel::internal_set_bulge_energies( FILE *fp, char *buffer )
{
  if( ptype == VIENNA)
    internal_read_array_data( fp, buffer, buffer, bulge_37_dG, 31 );
  else if( ptype == MFOLD )
    {
      hairpin_37_dG[0] = 0;
      internal_read_array_data( fp, buffer, buffer, &bulge_37_dG[1], 30 );
    }
}

void ViennaEnergyModel::internal_set_interior_loop_energies( FILE *fp, char *buffer )
{
  if( ptype == VIENNA)
    internal_read_array_data( fp, buffer, buffer, internal_37_dG, 31 );
  else if( ptype == MFOLD )
    {
      hairpin_37_dG[0] = 0;
      internal_read_array_data( fp, buffer, buffer, &internal_37_dG[1], 30 );
    }
}

void ViennaEnergyModel::internal_set_interior_1_1_energies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3;
  char *cur_bufspot;

  if( ptype == VIENNA )
    {
      fgets( buffer, 2048, fp );
      
      cur_bufspot = buffer;
      for( loop = 1; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        {
          for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA; loop2++ )
            {	
              cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_1_1_37_dG[loop][loop2][0][0], NUM_BASES * NUM_BASES);
            }
        }
    }
  else if( ptype == MFOLD )
    {
      while( buffer[0] == '>' ) fgets( buffer, 2048, fp );
      
      cur_bufspot = buffer;

      for( loop = 1; loop < NUM_BASEPAIRS_VIENNA -1; loop++ )
        for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA -1; loop2++ )
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
}

void ViennaEnergyModel::internal_set_interior_1_1_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2;
  char *cur_bufspot;
  fgets( buffer, 2048, fp );
  cur_bufspot = buffer;
  for( loop = 1; loop < NUM_BASEPAIRS_VIENNA; loop++ )
    {
      for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA; loop2++ )
        {	
          cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_1_1_37_dH[loop][loop2][0][0], NUM_BASES * NUM_BASES);
        }
    }
}

void ViennaEnergyModel::internal_set_interior_2_1_energies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3,loop4;
  char *cur_bufspot;
  
  if( ptype == VIENNA)
    {

      fgets( buffer, 2048, fp );
      cur_bufspot = buffer;
      for( loop = 1; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        {
          for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA; loop2++ )
            {	
              for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
                {
                  cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_2_1_37_dG[loop][loop2][loop3][0][0], NUM_BASES * NUM_BASES);
                }
            }
        }
    }
  else if( ptype == MFOLD )
    {
      while( buffer[0] == '>' ) fgets( buffer, 2048, fp );
      
      cur_bufspot = buffer;

      for( loop = 1; loop < NUM_BASEPAIRS_VIENNA -1; loop++ )
        for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA -1; loop2++ )
          {
            if( buffer[0] == 'A' || buffer[0] == 'T' || buffer[0] == 'G' || buffer[0] == 'C') 
              fgets( buffer, 2048, fp); // eat the lines with AT..AT, etc, just in case.
            for( loop3 = 1; loop3 < NUM_BASES; loop3++ )
              {
		
		
                for( loop4 = 0; loop4 < NUM_BASES; loop4++ )
                  {
                    internal_2_1_37_dG[loop][loop2][loop3][loop4][0] = 0;
                    internal_2_1_37_dG[loop][loop2][loop3][0][loop4] = 0;
                  }
                for( loop4 = 1; loop4 < NUM_BASES; loop4++ )
                  cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_2_1_37_dG[loop][loop2][loop3][loop4][1], NUM_BASES-1);
		
              }
          } 
    }
}

void ViennaEnergyModel::internal_set_interior_2_1_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3;
  char *cur_bufspot;
  fgets( buffer, 2048, fp );
  cur_bufspot = buffer;
  for( loop = 1; loop < NUM_BASEPAIRS_VIENNA; loop++ )
    {
      for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA; loop2++ )
        {	
          for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
            {
              cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_2_1_37_dH[loop][loop2][loop3][0][0], NUM_BASES * NUM_BASES);
            }
        }
    }
}

void ViennaEnergyModel::internal_set_interior_2_2_energies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3,loop4,loop5;
  char *cur_bufspot;
  if( ptype == VIENNA )
    {
      fgets( buffer, 2048, fp );
      cur_bufspot = buffer;
      for( loop = 1; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        {
          for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA; loop2++ )
            {	
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

  if( ptype == MFOLD )
    {
      while( buffer[0] == '>')
        fgets( buffer, 2048, fp );


      cur_bufspot = buffer;
      for( loop = 1; loop < NUM_BASEPAIRS_VIENNA-1; loop++ )
        {
          for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA-1; loop2++ )
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
}

void ViennaEnergyModel::internal_set_interior_2_2_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3,loop4,loop5;
  char *cur_bufspot;
  fgets( buffer, 2048, fp );
  cur_bufspot = buffer;
  for( loop = 1; loop < NUM_BASEPAIRS_VIENNA; loop++ )
    {
      for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA; loop2++ )
        {	
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

void ViennaEnergyModel::internal_set_dangle_5_energies( FILE *fp, char *buffer )
{
  if( ptype == VIENNA )
    {
      buffer[8] = '\0';
  
      internal_read_array_data( fp, buffer, buffer, &dangle_5_37_dG[0][0], NUM_BASEPAIRS_VIENNA * NUM_BASES );
    }
  else if( ptype == MFOLD )
    {
      while( buffer[0] == '>' ) fgets(buffer, 2048, fp);

      int loop,loop2;
      for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        for( loop2 = 0; loop2 < NUM_BASES; loop2++)
          dangle_5_37_dG[loop][loop2] = 0;

      for( loop = 1; loop < NUM_BASEPAIRS_VIENNA; loop++)
        internal_read_array_data( fp, buffer, buffer, &dangle_5_37_dG[loop][1], NUM_BASES -1);      
    }
}

void ViennaEnergyModel::internal_set_dangle_3_energies( FILE *fp, char *buffer )
{
  if( ptype == VIENNA )
    {
      buffer[8] = '\0';
  
      internal_read_array_data( fp, buffer, buffer, &dangle_3_37_dG[0][0], NUM_BASEPAIRS_VIENNA * NUM_BASES );
    }
  else if( ptype == MFOLD )
    {
      while( buffer[0] == '>' ) fgets(buffer, 2048, fp);

      int loop,loop2;
      for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        for( loop2 = 0; loop2 < NUM_BASES; loop2++)
          dangle_3_37_dG[loop][loop2] = 0;

      for( loop = 1; loop < NUM_BASEPAIRS_VIENNA; loop++)
        internal_read_array_data( fp, buffer, buffer, &dangle_3_37_dG[loop][1], NUM_BASES -1);      
    }
}

void ViennaEnergyModel::internal_set_dangle_5_enthalpies( FILE *fp, char *buffer )
{
  buffer[8] = '\0';
  internal_read_array_data( fp, buffer, buffer, &dangle_5_37_dH[0][0], NUM_BASEPAIRS_VIENNA * NUM_BASES );
}

void ViennaEnergyModel::internal_set_dangle_3_enthalpies( FILE *fp, char *buffer )
{
  buffer[8] = '\0';
  internal_read_array_data( fp, buffer, buffer, &dangle_3_37_dH[0][0], NUM_BASEPAIRS_VIENNA * NUM_BASES );
}

void ViennaEnergyModel::internal_set_multiloop_parameters( FILE *fp, char *buffer )
{
  int temp[4];
  if( ptype == VIENNA )
    {
      internal_read_array_data( fp, buffer, buffer, temp, 4 );
      multiloop_base = temp[0];
      multiloop_closing = temp[1];
      multiloop_internal = temp[2];
      terminal_AU = temp[3];
    }
  else if( ptype == MFOLD )
    {
      while( buffer[0] == '>') fgets(buffer, 2048, fp);

      internal_read_array_data( fp, buffer, buffer, temp, 3 );
      multiloop_base = temp[2];
      multiloop_closing = temp[0];
      multiloop_internal = temp[1];
    }
}

void ViennaEnergyModel::internal_set_at_penalty( FILE *fp, char *buffer )
{
  if( ptype == MFOLD ) // Vienna parameter set includes this in the multiloop data.
    {
      while( buffer[0] == '>' ) fgets( buffer, 2048, fp);

      internal_read_array_data( fp, buffer, buffer, &terminal_AU, 1 );
    }
}


void ViennaEnergyModel::internal_set_bimolecular_penalty( FILE *fp, char *buffer )
{
  if( ptype == MFOLD ) // Vienna parameter set doesn't have this term.
    {
      while( buffer[0] == '>' ) fgets( buffer, 2048, fp);

      internal_read_array_data( fp, buffer, buffer, &bimolecular_penalty, 1 );
    }
}

void ViennaEnergyModel::internal_set_ninio_parameters( FILE *fp, char *buffer )
{
  int temp[5];
  if( ptype == VIENNA )
    {
      internal_read_array_data( fp, buffer, buffer, temp, 2 );
      ninio_correction_37[2] = temp[0];
      maximum_NINIO = temp[1];
    }
  else if( ptype == MFOLD )
    {
      while( buffer[0] == '>')
        fgets(buffer, 2048, fp);
      internal_read_array_data( fp, buffer, buffer, temp, 5 );
      maximum_NINIO = temp[4];
      ninio_correction_37[0] = temp[0];
      ninio_correction_37[1] = temp[1];
      ninio_correction_37[2] = temp[2];
      ninio_correction_37[3] = temp[3];
    }
}

void ViennaEnergyModel::internal_set_hairpin_tetraloop_parameters( FILE *fp, char *buffer )
{
  int buf_index = 0;
  int lookup_index = 0;

  fgets( buffer, 2048, fp );

  // NOTE:: 4096 = (NUM_BASES-1)^6
  // Initialize the tetraloop parameters to 0. 
  for( int loop = 0; loop < 4096; loop ++ )
    {
      hairpin_tetraloop_37_dG[loop] = 0;
    }

  while( strlen(buffer) > 7  && buffer[0] != '>')
    {
      buf_index = 0;
      while( isspace(buffer[buf_index])) buf_index++;

      lookup_index = ((baseLookup(buffer[buf_index+0])-1) << 10) + 
        ((baseLookup(buffer[buf_index+1])-1) << 8) + 
        ((baseLookup(buffer[buf_index+2])-1) << 6) + 
        ((baseLookup(buffer[buf_index+3])-1) << 4) +
        ((baseLookup(buffer[buf_index+4])-1) << 2) +
        (baseLookup(buffer[buf_index+5])-1);
      hairpin_tetraloop_37_dG[ lookup_index ] = atoi( &buffer[buf_index + 6]);

      fgets(buffer, 2048, fp );
    }
#ifdef DEBUG
  fprintf(stderr,"Tetraloop Paramaters (VIENNA): %d read.\n",tetra_index);
#endif
}

void ViennaEnergyModel::internal_set_hairpin_triloop_parameters( FILE *fp, char *buffer )
{
  int buf_index = 0;
  int lookup_index = 0;
  fgets( buffer, 2048, fp );

  // NOTE:: 1024 = (NUM_BASES-1)^5
  // Initialize the triloop parameters to 0. 
  for( int loop = 0; loop < 1024; loop ++ )
    {
      hairpin_triloop_37_dG[loop] = 0;
    }
  while( strlen(buffer) > 6 && buffer[0] != '>')
    {
      buf_index = 0;
      while( isspace(buffer[buf_index])) buf_index++;

      lookup_index = ((baseLookup(buffer[buf_index+0])-1) << 8) + 
        ((baseLookup(buffer[buf_index+1])-1) << 6) + 
        ((baseLookup(buffer[buf_index+2])-1) << 4) + 
        ((baseLookup(buffer[buf_index+3])-1) << 2) +
        (baseLookup(buffer[buf_index+4])-1);

      hairpin_triloop_37_dG[ lookup_index ] = atoi( &buffer[buf_index + 5]);

      fgets( buffer, 2048, fp );
    }
#ifdef DEBUG
  fprintf(stderr,"Triloop Paramaters (VIENNA): %d read.\n",tri_index);
#endif
}


void ViennaEnergyModel::internal_set_hairpin_mismatch_energies( FILE *fp, char *buffer )
{
  int loop, loop2,loop3;
  char *cur_bufspot;
  
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
	{
	  if( loop == 0 || loop == NUM_BASEPAIRS_VIENNA-1)
	    for( loop2 = 0; loop2 < NUM_BASES * NUM_BASES; loop2++ )
	      hairpin_mismatch_37_dG[loop][0][loop2] = 0;
	  else
	    {
	      cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &hairpin_mismatch_37_dG[loop][0][0], NUM_BASES * NUM_BASES);
	    }
	}

}

void ViennaEnergyModel::internal_set_interior_loop_mismatch_energies( FILE *fp, char *buffer )
{
  int loop, loop2, loop3;
  char *cur_bufspot;

  
  if( ptype == MFOLD )
    {
      int temp[NUM_BASEPAIRS_VIENNA-2];
      while( buffer[0] == '>' ) fgets(buffer, 2048, fp);
      cur_bufspot = buffer;
      for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
          for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
            internal_mismatch_37_dG[loop][loop2][loop3] = 0;

      for( loop = 0; loop < (NUM_BASES-1) * (NUM_BASES-1); loop++)
        {
          cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &temp[0], NUM_BASEPAIRS_VIENNA-2);
          loop3 = (loop-(loop % (NUM_BASES-1))) / (NUM_BASES-1);
          for( loop2 = 1; loop2 < NUM_BASEPAIRS_VIENNA-1; loop2++ )
            internal_mismatch_37_dG[loop2][loop3+1][(loop % (NUM_BASES-1))+1] = temp[loop2-1];

        }
    
    }

  if( ptype == VIENNA )
    {
      cur_bufspot = buffer;
      for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        {
          if( loop == 0 || loop == NUM_BASEPAIRS_VIENNA-1)
            for( loop2 = 0; loop2 < NUM_BASES * NUM_BASES; loop2++ )
              internal_mismatch_37_dG[loop][0][loop2] = 0;
          else
            {
              cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &internal_mismatch_37_dG[loop][0][0], NUM_BASES * NUM_BASES);
            }
        }
    }
}

void ViennaEnergyModel::internal_set_multiloop_mismatch_energies( FILE *fp, char *buffer )
{
  int loop, loop2;
  char *cur_bufspot;
  /*
    if( ptype == MFOLD )
    {
    int temp[NUM_BASEPAIRS_VIENNA-2];
    cur_bufspot = buffer;
    for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
	for( loop2 = 0; loop2 < NUM_BASES; loop2++ )
    for( loop3 = 0; loop3 < NUM_BASES; loop3++ )
    multi

    for( loop = 0; loop < (NUM_BASES-1) * (NUM_BASES-1); loop++ )
	{
    if( loop == 0 || loop == NUM_BASEPAIRS_VIENNA-1)
    for( loop2 = 0; loop2 < NUM_BASES * NUM_BASES; loop2++ )
    multiloop_mismatch_37_dG[loop][0][loop2] = 0;
    else
    {
    cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &multiloop_mismatch_37_dG[loop][0][0], NUM_BASES * NUM_BASES);
    }
	}
    }
  */
  if( ptype == VIENNA )
    {  
      cur_bufspot = buffer;
      for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
        {
          if( loop == 0 || loop == NUM_BASEPAIRS_VIENNA-1)
            for( loop2 = 0; loop2 < NUM_BASES * NUM_BASES; loop2++ )
              multiloop_mismatch_37_dG[loop][0][loop2] = 0;
          else
            {
              cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &multiloop_mismatch_37_dG[loop][0][0], NUM_BASES * NUM_BASES);
            }
        }
    }
}

void ViennaEnergyModel::internal_set_mismatch_enthalpies( FILE *fp, char *buffer )
{
  int loop, loop2;
  char *cur_bufspot;
  
  cur_bufspot = buffer;
  for( loop = 0; loop < NUM_BASEPAIRS_VIENNA; loop++ )
    {
      if( loop == 0 || loop == NUM_BASEPAIRS_VIENNA-1)
        for( loop2 = 0; loop2 < NUM_BASES * NUM_BASES; loop2++ )
          mismatch_37_dH[loop][0][loop2] = 0;
      else
        {
          cur_bufspot = internal_read_array_data( fp, buffer, cur_bufspot, &mismatch_37_dH[loop][0][0], NUM_BASES * NUM_BASES);
        }
    }
}



char *ViennaEnergyModel::internal_read_array_data( FILE *fp, char *buffer, char *start_loc, int *read_loc, int size )
{
  int loop,loop2;
  char *cur_bufspot, *temp_char;
  int temp_int;
  int last = -1;

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
                {
                  read_loc[loop] = read_loc[last] + (int) rint( log_loop_penalty_37 * log( ((double) loop ) / (double) last));
                  
                  // See below note for why this is commented out.
                  //read_loc[loop] = read_loc[loop-1] + (int) rint( log_loop_penalty_37 * log(((double) loop) / ((double) (loop - 1)))) ;
                  //printf("index %d: %d vs %d\n", loop, read_loc[loop], other );

                  // A /very/ very subtle rounding error occurs in the
                  // integer arithmetic between the first version of
                  // this extrapolation, and the second.
                  //
                  // Note that 'last' is the index of the last non-'x' read in value.
                  //
                  // Spot the difference? They're actually equivalent ways of performing
                  // the computation, one requiring a memory and the other does not...
                  //
                  // However, they're equivalent only on the real numbers, not on the integers.
                  // On the dna parameter set, whose last read-in value is 420, at position 8,
                  // we get the following sequence:
                  //         (last)  (loop-1)
                  // index 9: 433 vs 433
                  // ...
                  // index 14: 480 vs 480
                  // index 15: 488 vs 487
                  // index 16: 495 vs 494
                  // index 17: 501 vs 501
                  // ...
                  // index 23: 534 vs 534
                  // index 24: 538 vs 539
                  // index 25: 543 vs 543
                  // ...
                  // This led to some very odd off by .01 errors. May they never show up again,
                  // and may we all use floating point numbers; no more shall we worry about 
                  // relative speeds of int vs double.

                }

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
          last = loop;
        }
    }

  return cur_bufspot;
}





