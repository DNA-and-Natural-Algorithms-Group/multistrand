/*
   Copyright (c) 2007-2008 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/
 
 /* StrandComplex class header. The Complex object contains information about a collection of strands, and has the loop structures contained within it. */
#ifndef __SCOMPLEX_H__
#define __SCOMPLEX_H__


#include "loop.h"
#include <string>


// structure containing information about bases exterior to the complex, IE bases that could pair with other complexes. First incarnation of such.a structure, prolly will change as I work out multiple-complex issues.
struct exterior_bases
{
  int A,T,C,G;
};


#include "strandordering.h"
#include "optionlists.h"

//#include "movetree.h"



class StrandComplex
{
 public:
  // Constructors. Still not sure on exactly how I want to do these. For now, they take a character sequence and structure.
  StrandComplex( char *seq, char *struc );
  StrandComplex( char *seq, char *struc, class identList *id_list );

  StrandComplex( StrandOrdering *newOrdering );

  // Destructor, basic.
  ~StrandComplex( void );
  void cleanup( void );


  // information retreival functions
  double getTotalFlux( void ); // returns total flux for all moves within the complex
  int getStrandCount( void ); // # of strands in the complex.
  double getEnergy( void ); // returns the energy of the complex
  void moveDisplay( void ); // display function to output the dot-paren structure of all moves contained in this complex. Should be preceded by printing the sequence, possibly I should change it to just do that straight out. Used for testing purposes (comparing all moves adjacent and rates).
  char *getSequence( void ); // returns char representation of sequence
  char *getStructure( void );// returns dot-paren notation structure for seq.
  char *getStrandNames( void ); // returns ordered list of strand names
  struct exterior_bases *getExteriorBases( void );
  int checkIDList( class identList *stoplist, int id_count );
  int checkIDBound( char *id );
  
  // published functions to affect the complex, these being a choice being made on the move set inside the complex, usually.
  // Once these are working, will need to add functions to merge complexes and perhaps others. Also, performing a choice will need to be able to pop a disassociation event back to the main system. Maybe do this with exceptions?
  // 11/24 JMS: I think I need to add a strand class in order to cleanly handle disassociation events for single strands. Disassociations that result in two seperate complexes needs to be handled as well, and efficiently. More thought required.
  // 11/25 JMS: Possibly the best thing to do is have the complex which performs the splitting choice return the new complex. If I implement strands, it should be easier to find the splitting point and construct the new complex efficiently.

  Move *getChoice( double *rand_choice ); // get a move chosen stochasticly from all possible moves within the complex. We'll then call perform choice on that move to generate the new setup.
  StrandComplex *doChoice( Move *move );
  int generateLoops( void );

  string toString(void );

  static StrandComplex *performComplexJoin( StrandComplex **complexes, char *types, int *index);

 private:
  StrandOrdering *ordering;
  Loop *beginLoop;
  MoveTree *kineticMoves;
  double totalFlux; // Total flux contained within this complex.

};



#endif
