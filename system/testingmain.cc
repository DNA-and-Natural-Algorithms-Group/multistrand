/*
  =======================================================================
  Copyright (c) 2007-2015 California Institute of Technology.
  Distributed under the MIT License.
  (See accompanying file LICENSE or copy at
  http://opensource.org/licenses/MIT)
  =======================================================================

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
