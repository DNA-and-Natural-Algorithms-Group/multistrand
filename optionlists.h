/*
   Copyright (c) 2007-2008 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/
 
#ifndef __OPTIONLISTS_H__
#define __OPTIONLISTS_H__

#define STOPTYPE_STRUCTURE                     0 
#define STOPTYPE_BOUND                         1
#define STOPTYPE_DISASSOC                      2
#define STOPTYPE_LOOSE_STRUCTURE               3
#define STOPTYPE_PERCENT_OR_COUNT_STRUCTURE    4

#include <stdlib.h>
#include <cstring>

// The bodies for these functions are defined in options.cc, but the
// actual classes get used across a bunch of the files so they have been
// seperated to not cause huge recompiles.


class strandlist {
 public:
  char *id;
  char *seq;
  class strandlist *next;
  char *lookup( char *item_id );
  char *name_lookup( char *item_id);
  strandlist( char *newid, char *newseq, class strandlist *old = NULL );
  ~strandlist( void );
};


class identlist {
 public:
  char *id;
  class identlist *next;
  identlist( char *newid, class identlist *old = NULL);
  void  make_unique( strandlist *strands);
  ~identlist( void );
};

class complex_item {
 public:
  char *structure;
  int type;
  int count; // for use with percentage or count stop types
  class identlist *strand_ids;
  class complex_item *next;
  complex_item( char *struc, class identlist *strands, class complex_item *old = NULL );
  complex_item( char *struc, class identlist *strands, class complex_item *old, int newtype );
  complex_item( char *struc, class identlist *strands, class complex_item *old, int newtype, int newcount );
  ~complex_item( void );
};

class stopcomplexes {
 public:
  char *tag;
  class complex_item *citem;
  class stopcomplexes *next;
  stopcomplexes( char *newtag, class complex_item *newitem, class stopcomplexes *old = NULL);
  ~stopcomplexes( void );
};

#endif
//  __OPTIONLISTS_H__
