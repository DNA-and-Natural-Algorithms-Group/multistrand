/*
   Copyright (c) 2007-2008 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/
 
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "energymodel.h"
#include "scomplex.h"
#include "ssystem.h"

#define DEBUG

/* ------------------------------------------------------------------------
   

   Testing Main


   ------------------------------------------------------------------------ */

/*
int main( void )
{
  char insequence[1024],instructure[1024];
  FILE *fp = fopen("dna.par","rt");
  if ( fp == NULL )
    printf("fail\n");
  else
    {
      EnergyModel *temp = NULL;
      temp = new EnergyModel( fp );
      Loop::SetEnergyModel( temp );

      while( !feof(stdin) )
	{
	  fgets( insequence, 1024, stdin );
	  while( insequence[0] == '>' )
	    fgets( insequence, 1024, stdin );
	  fgets( instructure, 1024, stdin );
	  if( feof( stdin)) continue;
	  insequence[strlen(insequence)-1] = '\0';
	  instructure[strlen(instructure)-1] = '\0';
	  StrandComplex *newcomplex = new StrandComplex( insequence, instructure);

	  newcomplex->generateLoops();
	  double energy = newcomplex->getEnergy();
	  
	  printf("%s\n%s (%6.2f)\n",insequence,instructure,energy);

	  delete newcomplex;
	}



      delete temp;
      temp = NULL;
    }
  printf("%d\n",sizeof(class Move));
}
*/


int main( int argc, char **argv )
{
  SimulationSystem *ssystem;
  ssystem = new SimulationSystem( argc, argv );
  ssystem->StartSimulation();
  delete ssystem;

  /*
     // Next block commented out. Pre-Options based code.
     // Works fine though, just not as useful.
  int numsims;
  double simtime;


  ssystem = new SimulationSystem("dna.par");
  ssystem->LoadSystem( stdin );
 
  int iflag = 0;
  numsims = 1;
  simtime = 1000.0;

  for( int loop = 1; loop < argc; loop++ )
    {
      if( !strcmp( argv[loop], "--energy" )) iflag = 1;
      else if( !strcmp( argv[loop], "--single" )) iflag = 2;
      else if( !strcmp( argv[loop], "--mc" )) iflag = 5;
      else if( !strcmp( argv[loop], "--mc1" )) iflag = 4;
      else if( !strcmp( argv[loop], "--silent" )) iflag = 6;
      else if( !strcmp( argv[loop], "--num" )) 
	{
	  if( argc> loop+1)
	    numsims = atoi(argv[loop+1]);
	}
      else if( !strcmp( argv[loop], "--time" )) 
	{
	  if( argc> loop+1)
	    simtime = atof(argv[loop+1]);
	}

    }

  if( iflag == 0 ) iflag = 3;

  
  ssystem->StartSimulation( iflag, numsims, simtime );
  // End pre-options commented code.
  */ 

 /* 
 // Pre-ssystem loading and init code. 

 char insequence[3072],instructure[3072];
  Move *tempmove;
  int iflag;
  FILE *fp = fopen("dna.par","rt");
  if ( fp == NULL )
    printf("fail\n");
  else
    {
      EnergyModel *temp = NULL;
      temp = new EnergyModel( fp );
      Loop::SetEnergyModel( temp );
      
      if( argc > 1 )
	{
	  if( !strcmp( argv[1], "--energy" )) iflag = 1;
	  if( !strcmp( argv[1], "--single" )) iflag = 2;
	  if( !strcmp( argv[1], "--mc" )) iflag = 4;
	  else iflag = 3;
	}
      else iflag = 3;

      srand( time( NULL ) );

      while( !feof(stdin) )
	{
	  fgets( insequence, 3072, stdin );
	  while( insequence[0] == '>' )
	    fgets( insequence, 3072, stdin );
	  fgets( instructure, 3072, stdin );
	  if( feof( stdin)) continue;
	  insequence[strlen(insequence)-1] = '\0';
	  instructure[strlen(instructure)-1] = '\0';
	  StrandComplex *newcomplex = new StrandComplex( insequence, instructure);

	  newcomplex->generateLoops();
	  double energy = newcomplex->getEnergy();
	  newcomplex->moveDisplay();	  
	  double rate = newcomplex->getTotalFlux();
	  printf("%s\n%s (%6.2f) %6.4f\n",insequence,instructure,energy,rate);
	  double stime = 0.0;
	  if( iflag > 1 )
	    {
	      if( iflag == 2 )
		{
		  double rchoice = (rate*rand()/((double)RAND_MAX));

		  printf("%3.3f\n",rchoice);
		  tempmove = newcomplex->getChoice( &rchoice );
		  newcomplex->doChoice(tempmove);
		  char *struc = newcomplex->getStructure();
		  rate = newcomplex->getTotalFlux();
		  energy = newcomplex->getEnergy();
		  printf("%s\n%s (%6.2f) %6.4f\n",insequence,struc,energy,rate);
		}
	      if( iflag == 3 )
		{
		  double rchoice,rsave,moverate;
		  char *struc;
		  int temp1, temp2,type;
		  do {
		    rsave = rchoice = (rate * rand()/((double)RAND_MAX));
		    stime += (log(1. / (double)((rand() + 1)/(double)RAND_MAX)) / rate );
		    tempmove = newcomplex->getChoice( &rchoice );
		    moverate = tempmove->getRate();
		    type = tempmove->getType();
		    newcomplex->doChoice( tempmove );
		    struc = newcomplex->getStructure();
		    rate = newcomplex->getTotalFlux();
		    energy = newcomplex->getEnergy();
		    printf("%s (%6.2f) %6.4f %6.4f pc: %6.4f time: %6.4f type: %d\n",struc,energy,rate,moverate,rsave,stime,type);
		    
		    for( int loop = 0; loop < strlen(struc); loop++)
		      {
			if(struc[loop] == '(') temp1++;
			if(struc[loop] == ')') temp1--;
			assert( temp1 >= 0 );
			if(loop > 0)
			  assert( !(struc[loop-1]=='(' && struc[loop] == ')'));
		      }

		    assert( temp1 == 0 );
		    // printf("%p",tempmove);
		  } while( rate > 0.10 && stime < 1000);
		}
	    }
	  delete newcomplex;
	}



      delete temp;
      temp = NULL;
      }
  // End old code block.
  */
}
