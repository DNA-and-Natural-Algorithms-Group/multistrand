/*
  =======================================================================
  Copyright (c) 2010-2015 California Institute of Technology.
  Distributed under the MIT License.
  (See accompanying file LICENSE or copy at
  http://opensource.org/licenses/MIT)
  =======================================================================

  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
  Edits by: Chris Berlind    (cberlind@dna.caltech.edu)
*/

 
#include <string.h>
#include "optionlists.h"

/*
  strandlist::strandlist( char *newid, char *newseq, class strandlist *old = NULL )

*/

strandlist::strandlist( char *newid, char *newseq, class strandlist *old )
{
  id = newid;
  seq = newseq;
  next = old;
}

char *strandlist::lookup( char *item_id )
{
  if( strcmp( item_id, id ) == 0) return seq;
  else if( next != NULL ) return next->lookup(item_id);
  else return NULL;
}

/*
  strandlist::~strandlist( void )

  This is a basic linked list destructor. The "id" and "seq" members were allocated in the parser via strdup, and are thus free'd here.

*/

strandlist::~strandlist( void )
{
  if( id != NULL )
    delete id;  
  if( seq != NULL )
    delete seq;
  if( next != NULL )
    delete next;
  next = NULL;
}

/*
  identlist::identlist( char *newid, class identlist *old = NULL )

  This is a very basic linked list constructor.
*/

identlist::identlist( long newuid, char *newid, class identlist *old )
{
  id = new char[strlen(newid) + 1];
  strcpy(id, newid );

  uid = newuid;
  next = old;
}

/*
  identlist::~identlist( void )

  This is a basic linked list destructor. The "id" parameter is assumed to be
  deleted here.

*/

identlist::~identlist( void )
{
  if( id != NULL )
    delete[] id;  
  if( next != NULL )
    delete next;
  next = NULL;
}




/*
  complex_item::complex_item( char *struc, class identlist *strands, class complex_item *old) 

  This is another very basic linked list constructor.
*/


complex_item::complex_item( char *struc, class identlist *strands, class complex_item *old) 
{
  structure = new char[strlen(struc) + 1];
  strcpy(structure, struc );
  
  while( strchr( structure, '_' ) != NULL )
    {
      *strchr(structure,'_') = '+';
    }
  // while( strchr( structure, '+' ) != NULL )
  //{
  //  *strchr(structure,'+') = '.';
  //}
  strand_ids = strands;
  next = old;
  type = STOPTYPE_STRUCTURE;
  count = 0;
}


complex_item::complex_item( char *struc, class identlist *strands, class complex_item *old, int newtype) 
{
  structure = new char[strlen(struc) + 1];
  strcpy(structure, struc );
  
  if( structure != NULL )
    {
      while( strchr( structure, '_' ) != NULL )
	{
	  *strchr(structure,'_') = '+';
	}
      //      while( strchr( structure, '+' ) != NULL )
      //{
      //  *strchr(structure,'+') = '.';
      //}
    }
  strand_ids = strands;
  next = old;
  type = newtype;
  count = 0;
}

complex_item::complex_item( char *struc, class identlist *strands, class complex_item *old, int newtype, int newcount) 
{
  structure = new char[strlen(struc) + 1];
  strcpy(structure, struc );
  
  if( structure != NULL )
    {
      while( strchr( structure, '_' ) != NULL )
	{
	  *strchr(structure,'_') = '+';
	}
      //      while( strchr( structure, '+' ) != NULL )
      //{
      //  *strchr(structure,'+') = '.';
      //}
    }
  strand_ids = strands;
  next = old;
  type = newtype;
  count = newcount;
}

/*
  complex_item::~complex_item( void )

  This is a basic linked list destructor. The "structure" member was strdup'd in the parser, and is deleted here via free.

*/

complex_item::~complex_item( void )
{
  if( structure != NULL )
    delete[] structure;
  if( strand_ids != NULL )
    delete strand_ids;
  if( next != NULL )
    delete next;
  next = NULL;
}




/*
  stopcomplexes::stopcomplexes( char *newtag, class complex_item *newitem, class stopcomplexes *old) 

  This is yet another very basic linked list constructor.
*/

stopcomplexes::stopcomplexes( char *newtag, class complex_item *newitem, class stopcomplexes *old) 
{
  tag = new char[strlen(newtag) + 1];
  strcpy(tag, newtag );
  
  citem = newitem;
  next = old;
}

/*
  stopcomplexes::~stopcomplexes( void )

  This is a basic linked list destructor. The "tag" member was strdup'd in the parser, and is deleted here via free.

*/

stopcomplexes::~stopcomplexes( void )
{
  if( tag != NULL )
    delete[] tag;
  if( citem != NULL )
    delete citem;
  if( next != NULL )
    delete next;
  next = NULL;
}
