/*
   Copyright (c) 2007-2010 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#include "energymodel.h"


EnergyModel::EnergyModel( PyObject *options )
{
  // nothing yet
}

EnergyModel::EnergyModel( void )
{
  // nothing yet either
}

EnergyModel::~EnergyModel( void )
{
  // nothing
}


int pairs[5] = {0,0,0,0,0};
int pairtypes[5][5] = {
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0}
};
int basepair_sw[8] = {0,0,0,0,0,0,0,0};


int lookuphelper[26] = {1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,4,4,0,0,0,0,0};
//                      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z

// // helper function to convert to numerical base format.
int baseLookup( char base )
{
  char temp = toupper(base);
  if( temp < 'A' || temp > 'Z' )
    return base;
  return lookuphelper[temp-'A'];
}

