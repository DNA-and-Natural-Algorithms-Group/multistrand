/*
   Copyright (c) 2007-2010 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/
 
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "ssystem.h"

#define DEBUG

/* ------------------------------------------------------------------------
   

   Testing Main


   ------------------------------------------------------------------------ */

int main( int argc, char **argv )
{
  SimulationSystem *ssystem;
  ssystem = new SimulationSystem( argc, argv );
  ssystem->StartSimulation();
  delete ssystem;
}



