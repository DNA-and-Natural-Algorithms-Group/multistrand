/*
 Copyright (c) 2007-2010 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 Edits by: Chris Berlind    (cberlind@dna.caltech.edu)
 */

#include <Python.h>

#include <string>
#include "optionlists.h"
#include "utility.h"

using std::string;

/*
 strandlist::strandlist( char *newid, char *newseq, class strandlist *old = NULL )

 */

strandList::strandList(char *newid, char *newseq, class strandList *old) {
	id = newid;
	seq = newseq;
	next = old;
}

char *strandList::lookup(char *item_id) {
	if (strcmp(item_id, id) == 0)
		return seq;
	else if (next != NULL)
		return next->lookup(item_id);
	else
		return NULL;
}

/*
 strandlist::~strandlist( void )

 This is a basic linked list destructor. The "id" and "seq" members were allocated in the parser via strdup, and are thus free'd here.

 */

strandList::~strandList(void) {
	if (id != NULL)
		delete id;
	if (seq != NULL)
		delete seq;
	if (next != NULL)
		delete next;
	next = NULL;
}

/*
 identlist::identlist( char *newid, class identlist *old = NULL )

 This is a very basic linked list constructor.
 */

identList::identList(long newuid, const char *newid, class identList *old) {
	id = new char[strlen(newid) + 1];
	strcpy(id, newid);

	uid = newuid;
	next = old;
}

std::string identList::toString() {

	if (id != NULL) {
		return utility::copyToString(id);
	} else {
		return string("");
	}

}

/*
 identlist::~identlist( void )

 This is a basic linked list destructor. The "id" parameter is assumed to be
 deleted here.

 */

identList::~identList(void) {
	if (id != NULL)
		delete[] id;
	if (next != NULL)
		delete next;
	next = NULL;
}

/*
 complex_item::complex_item( char *struc, class identlist *strands, class complex_item *old)

 FD: This is another very basic linked list constructor.
 */

complexItem::complexItem(const char *struc, class identList *strands,
		class complexItem *old) {
	structure = new char[strlen(struc) + 1];
	strcpy(structure, struc);

	while (strchr(structure, '_') != NULL) {
		*strchr(structure, '_') = '+';
	}

	strand_ids = strands;
	next = old;
	type = STOPTYPE_STRUCTURE;
	count = 0;
}

complexItem::complexItem(const char *struc, class identList *strands,
		class complexItem *old, int newtype) {
	structure = new char[strlen(struc) + 1];
	strcpy(structure, struc);

	if (structure != NULL) {
		while (strchr(structure, '_') != NULL) {
			*strchr(structure, '_') = '+';
		}

	}
	strand_ids = strands;
	next = old;
	type = newtype;
	count = 0;
}

complexItem::complexItem(const char *struc, class identList *strands,
		class complexItem *old, int newtype, int newcount) {
	structure = new char[strlen(struc) + 1];
	strcpy(structure, struc);

	if (structure != NULL) {
		while (strchr(structure, '_') != NULL) {
			*strchr(structure, '_') = '+';
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

complexItem::~complexItem(void) {
	if (structure != NULL)
		delete[] structure;
	if (strand_ids != NULL)
		delete strand_ids;
	if (next != NULL)
		delete next;
	next = NULL;
}

/*
 stopcomplexes::stopcomplexes( char *newtag, class complex_item *newitem, class stopcomplexes *old)

 This is yet another very basic linked list constructor.
 */

stopComplexes::stopComplexes(const char *newtag, class complexItem *newitem,
		class stopComplexes *old) {
	tag = new char[strlen(newtag) + 1];
	strcpy(tag, newtag);

	citem = newitem;
	next = old;
}

string stopComplexes::toString() {

	string output = "";

	return output;

}

string complexItem::toString() {

	string output = "";

	return output;

}

/*

 stopcomplexes::~stopcomplexes(void)

 This is a basic linked list destructor. The "tag" member was strdup'd in the parser, and is deleted here via free.

 */

stopComplexes::~stopComplexes(void) {
	if (tag != NULL)
		delete[] tag;
	if (citem != NULL)
		delete citem;
	if (next != NULL)
		delete next;
	next = NULL;
}
