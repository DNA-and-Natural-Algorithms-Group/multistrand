/*
   Copyright (c) 2007-2010 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#ifndef __OPTIONLISTS_H__
#define __OPTIONLISTS_H__

#define STOPTYPE_STRUCTURE                     0
#define STOPTYPE_BOUND                         1
#define STOPTYPE_DISASSOC                      2
#define STOPTYPE_LOOSE_STRUCTURE               3
#define STOPTYPE_PERCENT_OR_COUNT_STRUCTURE    4

#include <python2.7/Python.h>
#include <string>
// for PyObject *

// The bodies for these functions are defined in options.cc, but the
// actual classes get used across a bunch of the files so they have been
// seperated to not cause huge recompiles.


class strandList {
 public:
  char *id;
  char *seq;
  class strandList *next;
  char *lookup( char *item_id );
  char *name_lookup( char *item_id);
  strandList( char *newid, char *newseq, class strandList *old = NULL );
  ~strandList( void );
};


class identList {
 public:
  long uid;
  char *id;
  PyObject *pyo_id;  // needed for correct ref counting dealloc.
  class identList *next;
  identList( long newuid, char *newid, class identList *old = NULL);
  void  make_unique( strandList *strands);
  std::string toString(void);
  ~identList( void );
};

class complexItem {
 public:
  char *structure;
  int type;
  int count; // for use with percentage or count stop types
  class identList *strand_ids;
  class complexItem *next;
  complexItem( char *struc, class identList *strands, class complexItem *old = NULL );
  complexItem( char *struc, class identList *strands, class complexItem *old, int newtype );
  complexItem( char *struc, class identList *strands, class complexItem *old, int newtype, int newcount );
  ~complexItem( void );
};

class stopComplexes {
 public:
  char *tag;
  class complexItem *citem;
  class stopComplexes *next;
  stopComplexes( char *newtag, class complexItem *newitem, class stopComplexes *old = NULL);
  ~stopComplexes( void );
};

#endif
//  __OPTIONLISTS_H__
