/*
   Copyright (c) 2007-2008 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h> 
#include "options.h"
#include <ctype.h>
#include <string>  // for std::string
#include <cstring> // for strcpy on strings.


#define MAX_BUFFER 1024
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

char *strandlist::name_lookup( char *item_id )
{
  char * temp = strchr(id, ':');
  if( strcmp( temp+1, item_id ) == 0) return id;
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

identlist::identlist( char *newid, class identlist *old )
{
  id = newid;
  next = old;
}

void identlist::make_unique( strandlist *strands)
{
  char * temp;
  if (strchr(id, ':') != NULL)
    {
      if (next != NULL)
	next->make_unique(strands);
    }
  else
    {
      temp = strands->name_lookup(id);
      free(id);
      id = temp;
    }
}
/*
  identlist::~identlist( void )

  This is a basic linked list destructor. The "id" parameter is assumed to be
  deleted here.

*/

identlist::~identlist( void )
{
  if( id != NULL )
    free(id);  // id is created by strdup in the parser and so we must use free here.
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
  structure = struc;
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
  structure = struc;
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
  structure = struc;
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
    free(structure);  
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
  tag = newtag;
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
    free(tag);  
  if( citem != NULL )
    delete citem;
  if( next != NULL )
    delete next;
  next = NULL;
}

/*

   Options::Options( void )

   default constructor. Initializes the flag array and data members 
   appropriately. Can't do much else, really.

*/

/*

   Options::Options( char *filename)

   Also initializes the default logfile name.

*/


Options::Options( char *filename )
{

  char *search;
  flagarray = 0;
  inputfilename[0]='\0';
  logfilename[0]='\0';
  logfile = NULL;
  trajectoryfilename[0]='\0';
  trajfile = NULL;
  //  complexcount = 0;
  strandcount = 0;
  strands = NULL;
  start_structure = NULL;
  simulationtime = 10000.0;
  simulationmode = SIMULATION_NORMAL;
  trajectorycount = 1;
  outputinterval = -1;
  currentinterval = 0;
  stopoptions = 0;
  stoplist = NULL;
  stopcount = 0;
  initialseed = 0;
  //  joinrate = 1.0;
  //joinenergy = 0.0;
  //  joinconc = 1000.0;
  trajtype = 0;
  //  dangles = 0;
  gtenable = 0;
  intramolecularscaling = 1.6 * 1000000.0;
  intermolecularscaling =  .5 * 1000000.0;
  logml = -1;
  strcpy(energymodelfilename,"dna.par"); 
  energymodel_type = 0;
  energymodel = 0;
  outputtime = 0.0;
  energymode = 0;
  temperature = 37.0;
  ratemethod = 2; // kawasaki
  search = strstr( filename, ".txt");

  /* python interface data members */
  python_identlist = NULL;
  python_start_structure = NULL;
  python_stop_structure = NULL;
  python_strands = NULL;

  python_current_time = 0.0;
  python_current_seed = 0;
  python_trajectory_completion_flag = 0;
  python_suspend_trajectory_flag = 0; 
  python_halt_trajectory_flag = 0;
  python_trajectory_tag = "";
  python_trajectory_time = 0.0;
  python_k_collision = -1.0;


  
  /*Unique ID*/
   unique_id = 0;

  /* end python interface data members */
  if( search == NULL)
    {
      strcpy(logfilename,filename);
      strcat(logfilename,".log");
    }
  else
    {
      strcpy(logfilename,filename);
      logfilename[ search-filename] = '\0';
      strcat(logfilename,".log");
    }
}


/*

   Options::~Options( void )

   default destructor. Closes the file handle if it is open, and deletes
   the strand/complex/stoplist info.

*/

Options::~Options( void )
{
  if( python_identlist != NULL )
    delete python_identlist;
  if(logfile != NULL)
    fclose(logfile);
  if( strands != NULL )
    delete strands;
  if( start_structure != NULL )
    delete start_structure;
  if( stoplist != NULL )
    delete stoplist;
}


/* 
  int Options::addSequence( char *id, char *seq );
  int Options::addSequence( char *id, char *seq, double conc );
  int Options::addSequence( char *id, char *seq, int count );

*/


int Options::addSequence( char *id, char *seq )
{
 
  /*The following code was added 
   * to allow the insertion of unique IDs
   * for strands */
  short length = 1;

  char *id2 = new char[strlen(id)+12];
  char *seq2 = new char[strlen(seq)+1];
  // printf("%l:%s\n", unique_id, id);
  printf("%s", id);
  printf("%s\n", id2);
  sprintf( id2, "%ld:%s",unique_id, id);
  strcat(id2, id);
  strcpy(seq2, seq);
  
  free( id );
  free( seq );
  unique_id++;
  if( strands == NULL )
    strands = new strandlist( id2, seq2 );
  else
    strands = new strandlist( id2, seq2, strands );
  
  if( strands == NULL || id2 == NULL || seq2 == NULL)
    return 0; // error
  else
    return 1;
}

int Options::addSequence( char *id, char *seq, double conc )
{
  printf("Sequence added\n");
}


int Options::addSequence( char *id, char *seq, int count )
{
  printf("Sequence added\n");
}

int check_bases( char base1, char base2 )
{
  char bases[2];
  bases[0] = tolower(base1);
  bases[1] = tolower(base2);
  if( bases[0] > bases[1] )
    {
      bases[0] = bases[1];
      bases[1] = tolower(base1);
    }
  if( bases[0] == 'a' && (bases[1] == 't' || bases[1] == 'u') )
    return 1;
  if( bases[0] == 'c' && bases[1] == 'g' )
    return 1;
  if( bases[0] == 'g' && (bases[1] == 't' || bases[1] == 'u'))
    return 1;
  
  return 0;
  
}

int Options::addStopStructures( class stopcomplexes *stopstructs, int errortype )
{
  class stopcomplexes *traverse = stopstructs;
  if (stoplist != NULL )
    {
      if( errortype == 0 )
	fprintf(stderr,"WARNING: multiple stop structure lists given in input.\n          Only the final set of stop structures will be used.\n");
      delete stoplist;
    }
  stoplist = stopstructs;
  while( traverse != NULL )
    {
      stopcount++;
      traverse = traverse->next;
    }

  // error checking time!
  // For each stop structure given that is of the normal type, 
  //  1. Make sure total length is correct (strlen(strands)+strandcount-1)
  //  2. Make sure there are seperators between each strand
  //  3. Balanced parens check
  //  4. Maybe check that the bases are viable for balanced parens.
  //  5. Completely connected structure.
  //  6. Make sure all strand ids used match up with valid strands.
  //
  // Note that this is the same error checking for start structures. 

  traverse = stoplist;
  class complex_item *ctrav = stoplist->citem;
  int struct_len;
  int strand_count;
  class identlist *strav;
  char *copy_seq;
  int length_total, cur_strand, error_flag = 0, parencount, *pair_array;
  char *composed_sequence;
  for( int loop = 0; loop < stopcount; loop++, traverse = traverse->next )
    {
      ctrav = traverse->citem;
      if( ctrav->type != STOPTYPE_STRUCTURE )
	continue;   // TODO: error checking for bound/disassoc type stops.
      while( ctrav != NULL )
	{
	  strav = ctrav->strand_ids;
	  strav->make_unique(strands);
	  strand_count = 0;
	  while( strav != NULL )
	    {
	      strav = strav->next;
	      strand_count++;
	    }
	  //	  fprintf(stderr,"Strand Count: %d\n",strand_count);

	  struct_len = strlen(ctrav->structure);
	  
	  composed_sequence = new char[struct_len+1];
	  composed_sequence[0] = '\0';

	  pair_array = new int[struct_len];
	  for( int loop4 = 0; loop4 < struct_len; loop4++ )
	    pair_array[loop4]=0;

	  strav = ctrav->strand_ids;
	  length_total = 0;
	  while( strav != NULL )
	    {
	      copy_seq = strands->lookup(strav->id);
	      if( copy_seq == NULL ) // Check #6
		{
		  if( errortype == 0 )
		    {
		      fprintf(stderr,"ERROR: In StopStruc(%s), Strand ID %s not found in list of strands.\n",traverse->tag,strav->id);
		      fprintf(stderr,"ERROR: Can't continue checking for errors, halting.\n");
		      exit(0);
		    }
		  else if (errortype == 1)
		    return -1; // STRAND ID NOT FOUND
		  error_flag = 1;
		}
	      length_total += strlen(copy_seq);
	      //	      fprintf(stderr,"Len: %d Seq: %s\n",strlen(copy_seq),copy_seq);
	      if( length_total > struct_len )
		break;
	      strcat( composed_sequence, copy_seq);
	      if( strav->next != NULL )
		{
		  length_total ++;
		  if( length_total > struct_len )
		    break;
		  strcat( composed_sequence, "+");
		}
	      strav = strav->next;
	    }

	  // we now have, for this particular structure, the base sequence
	  // corresponding to the strand ids listed, with _ for gaps, and
	  // the lengths of each strand. 


	  if( struct_len != length_total ) // Check #1
	    {
	      if( errortype == 0)
		fprintf(stderr,"ERROR: StopStruc(%s): Strand length (%d) and structure length (%d) mismatched.\nERROR: Structure: %s\n",traverse->tag,length_total,struct_len,ctrav->structure);
	      else if( errortype == 1)
		return -2; // strand length and structure length mismatch
	      error_flag = 1;
	    }
	  else
	    {
	  
	  // Checks #2-5.
	  parencount = 0;
	  cur_strand = 0;
	  for( int loop2 = 0; loop2 < length_total; loop2++ )
	    {
	      if( composed_sequence[loop2] == '+' )
		{
		  if( ctrav->structure[loop2] == '(' || ctrav->structure[loop2] == ')') // Check #2: valid seperators, otherwise our paren counting could be off.
		    {
		      if( errortype == 0)
			fprintf(stderr,"ERROR: StopStruc(%s): Structure does not have a valid strand break character at strand break in position %d.\nERROR: Sequence : %s\nERROR: Structure: %s\n",traverse->tag,loop2+1,composed_sequence,ctrav->structure);
		      else if( errortype == 1)
			return -3; // structure does not have valid break character
		      error_flag = 1;
		      loop2 = length_total;
		      continue;
		    }
		  else if( parencount == 0 ) // check #5: connected structure
		    {
		      if( errortype == 0)
			fprintf(stderr,"ERROR: StopStruc(%s): Structure is not connected, break at position %d.\nERROR: Sequence : %s\nERROR: Structure: %s\n",traverse->tag,loop2+1,composed_sequence,ctrav->structure);
		      else if (errortype == 1)
			return -4;
		      error_flag = 1;
		      loop2 = length_total;
		      continue;
		  }
		}
	      if( ctrav->structure[loop2] == '(' )
		{
		  parencount++;
		  pair_array[loop2] = parencount;
		}
	      if( ctrav->structure[loop2] == ')' )
		{int loop3;
		  for( loop3 = loop2-1; loop3 >= 0; loop3-- )
		    if( pair_array[loop3] == parencount )
		      break;
		  parencount--;

		  if( parencount >= 0 && check_bases(composed_sequence[loop3], composed_sequence[loop2]) != 1)
		  {
		    if( errortype == 0 )
		      fprintf(stderr,"ERROR: StopStruc(%s): Base %c at position %d and Base %c at position %d cannot pair.\nERROR: Sequence : %s\nERROR: Structure: %s\n",traverse->tag, composed_sequence[loop3], loop3+1, composed_sequence[loop2], loop2+1,composed_sequence,ctrav->structure);
		    else if (errortype == 1 )
		      return -5;
		    error_flag = 1;
		  }
		}
	      if( parencount < 0 ) // check #3, part 1
		{
		  if( errortype == 0)
		    fprintf(stderr,"ERROR: StopStruc(%s): Mismatched ) paren at position %d.\nERROR: Sequence : %s\nERROR: Structure: %s\n",traverse->tag,loop2+1,composed_sequence,ctrav->structure);
		  else if ( errortype == 1 )
		    return -6; // mismatched paren
		  error_flag = 1;
		  loop2 = length_total;
		  continue;
		}
	    }
	  if( parencount > 0 ) // check #3, part 2
	    {
	      if( errortype == 0 )
		fprintf(stderr,"ERROR: StopStruc(%s): Mismatched parens, too many (.\nERROR: Sequence : %s\nERROR: Structure: %s\n",traverse->tag,composed_sequence,ctrav->structure);
	      else if( errortype == 1)
		return -7; // mismatched paren
	      error_flag = 1;
	    }
	    }

	  ctrav = ctrav->next;
	  delete[] composed_sequence;
	  delete[] pair_array;
	}
    }
  if( error_flag == 1 )
    {
      if( errortype == 0)
	{
	  fprintf(stderr,"ERROR: Too many errors in stop structures, halting.\n");
	  exit(0);
	}
      if( errortype == 1)
	return 0;
    }
  return 1;
}

int Options::addStartStructure( class complex_item *startstruct, int errortype)
{
  if (start_structure != NULL )
    {
      if( errortype == 0)
	fprintf(stderr,"WARNING: multiple start structures given in input.\n          Only the final start structure will be used.\n");
      delete start_structure;
    }
  start_structure = startstruct;

  // error checking time!
  // For each stop structure given, 
  //  1. Make sure total length is correct (strlen(strands)+strandcount-1)
  //  2. Make sure there are seperators between each strand
  //  3. Balanced parens check
  //  4. Maybe check that the bases are viable for balanced parens.
  //  5. Completely connected structure.
  //  6. Make sure all strand ids used match up with valid strands.
  //
  // Note that this is the same error checking for start structures. 

  class complex_item *ctrav = start_structure;
  int struct_len;
  int strand_count;
  class identlist *strav;
  char *copy_seq;
  int length_total, cur_strand, error_flag = 0, parencount, *pair_array;
  char *composed_sequence;
  ctrav = start_structure;
  while( ctrav != NULL )
    {
      strav = ctrav->strand_ids;
      strav->make_unique(strands);
      strand_count = 0;
      while( strav != NULL )
	{
	  strav = strav->next;
	  strand_count++;
	}
      //	  fprintf(stderr,"Strand Count: %d\n",strand_count);
      
      struct_len = strlen(ctrav->structure);
      
      composed_sequence = new char[struct_len+1];
      composed_sequence[0] = '\0';
      
      pair_array = new int[struct_len];
      for( int loop4 = 0; loop4 < struct_len; loop4++ )
	pair_array[loop4]=0;
      
      strav = ctrav->strand_ids;
      length_total = 0;
      while( strav != NULL )
	{
	  copy_seq = strands->lookup(strav->id);
	  if( copy_seq == NULL ) // Check #6
	    {
	      if( errortype == 0)
		{
		  fprintf(stderr,"ERROR: In Start Structure, Strand ID %s not found in list of strands.\n",strav->id);
		  fprintf(stderr,"ERROR: Can't continue checking for errors, halting.\n");
		}
	      else if( errortype == 1 )
		return -1;
	      exit(0);
	      error_flag = 1;
	    }
	  length_total += strlen(copy_seq);
	  //	      fprintf(stderr,"Len: %d Seq: %s\n",strlen(copy_seq),copy_seq);
	  if( length_total > struct_len )
	    break;
	  strcat( composed_sequence, copy_seq);
	  if( strav->next != NULL )
	    {
	      length_total ++;
	      if( length_total > struct_len )
		break;
	      strcat( composed_sequence, "+");
	    }
	  strav = strav->next;
	}
      
      // we now have, for this particular structure, the base sequence
      // corresponding to the strand ids listed, with _ for gaps, and
      // the lengths of each strand. 
      
      
      if( struct_len != length_total ) // Check #1
	{
	  if( errortype == 0)
	    fprintf(stderr,"ERROR: StartStruct: Strand length (%d) and structure length (%d) mismatched.\nERROR: Structure: %s\n",length_total,struct_len,ctrav->structure);
	  else if( errortype == 1)
	    return -2;
	  error_flag = 1;
	}
      else
	{
	  
	  // Checks #2-5.
	  parencount = 0;
	  cur_strand = 0;
	  for( int loop2 = 0; loop2 < length_total; loop2++ )
	    {
	      if( composed_sequence[loop2] == '+' )
		{
		  if( ctrav->structure[loop2] == '(' || ctrav->structure[loop2] == ')') // Check #2: valid seperators, otherwise our paren counting could be off.
		    {
		      if( errortype == 0)
			fprintf(stderr,"ERROR: StartStruct: Structure does not have a valid strand break character at strand break in position %d.\nERROR: Sequence : %s\nERROR: Structure: %s\n",loop2+1,composed_sequence,ctrav->structure);
		      else if( errortype == 1)
			return -3;
		      error_flag = 1;
		      loop2 = length_total;
		      continue;
		    }
		  else if( parencount == 0 ) // check #5: connected structure
		    {
		      if( errortype == 0)
			fprintf(stderr,"ERROR: StartStruct: Structure is not connected, break at position %d.\nERROR: Sequence : %s\nERROR: Structure: %s\n",loop2+1,composed_sequence,ctrav->structure);
		      else if( errortype == 1)
			return -4;
		      error_flag = 1;
		      loop2 = length_total;
		      continue;
		  }
		}
	      if( ctrav->structure[loop2] == '(' )
		{
		  parencount++;
		  pair_array[loop2] = parencount;
		}
	      if( ctrav->structure[loop2] == ')' )
		{int loop3;
		  for( loop3 = loop2-1; loop3 >= 0; loop3-- )
		    if( pair_array[loop3] == parencount )
		      break;
		  parencount--;

		  if( parencount >= 0 && check_bases(composed_sequence[loop3], composed_sequence[loop2]) != 1)
		  {
		    if( errortype == 0)
		      fprintf(stderr,"ERROR: StartStruct: Base %c at position %d and Base %c at position %d cannot pair.\nERROR: Sequence : %s\nERROR: Structure: %s\n", composed_sequence[loop3], loop3+1, composed_sequence[loop2], loop2+1,composed_sequence,ctrav->structure);
		    else if( errortype == 1)
		      return -5;
		    error_flag = 1;
		  }
		}
	      if( parencount < 0 ) // check #3, part 1
		{
		  if( errortype == 0)
		    fprintf(stderr,"ERROR: StartStruct: Mismatched ) paren at position %d.\nERROR: Sequence : %s\nERROR: Structure: %s\n",loop2+1,composed_sequence,ctrav->structure);
		  else if( errortype == 1)
		    return -6;
		  error_flag = 1;
		  loop2 = length_total;
		  continue;
		}
	    }
	  if( parencount > 0 ) // check #3, part 2
	    {
	      if( errortype == 0)
		fprintf(stderr,"ERROR: StartStruct(%s): Mismatched parens, too many (.\nERROR: Sequence : %s\nERROR: Structure: %s\n",composed_sequence,ctrav->structure);
	      else if( errortype == 1)
		return -7;
	      error_flag = 1;
	    }
	    }

	  ctrav = ctrav->next;
	  delete[] composed_sequence;
	  delete[] pair_array;
    }
  if( error_flag == 1 )
    {
      if( errortype == 0)
	{
	  fprintf(stderr,"ERROR: Too many errors in start structure, halting.\n");
	  exit(0);
	}
      // unknown error
      if( errortype == 1)
	return 0;
    }
  return 1;
}

/*  Simple mutators for various private members, used by the parser.

void Options::setSimulationTime( double newtime )
void Options::setSimulationMode( int mode )
void Options::setStopOptions( int newoptions )
void Options::setLogfile( char *logfilename )
void Options::setNumSimulations( int newnumsims )
void Options::setOutputInterval( int newoutputinterval )
void Options::setInitialSeed( long seedval )
void Options::setParameterFile( char *newparamfile )
void Options::setParameterType( int type )
void Options::setIntramolecularScaling( double newscale );
void Options::setIntermolecularScaling( double newscale );
void Options::setEnergyMode( int mode );
*/ 

void Options::setSimulationMode( int mode )
{
  if ( mode < SIMULATION_NORMAL || mode > SIMULATION_PYTHON_FIRST_BI )
    fprintf(stderr, "WARNING: Simulation Mode %d not recognized. Using normal simulation mode. \n", mode );
  simulationmode = mode;
}

void Options::incrementOutputState( void )
{
  if( currentinterval == outputinterval ) currentinterval = 0;
  else currentinterval++;

}

void Options::setSimulationTime( double newtime )
{
  if( newtime <= 0.0 )
    fprintf(stderr,"WARNING: Simulation Time entered as a non-positive value. Using defaults.\n");
  else
    simulationtime = newtime;
}

void Options::setStopOptions( int newoptions )
{
  // TODO: some error checking here once options are finalized.
  stopoptions = newoptions;
}

void Options::setLogfile( char *newlogfilename )
{
  if( strlen(newlogfilename) > 79 )
    fprintf(stderr,"WARNING: Logfile entry too long. Please use a shorter name.\n");
  else
    {
      strcpy(logfilename,newlogfilename);
      free(newlogfilename); // strdup'd in parser.

      // TODO: open the log file here?
    }
}


void Options::setNumSimulations( int newnumsims )
{
  if( newnumsims <= 0 )
    fprintf(stderr,"WARNING: Number of simulations must be positive. Using default values.\n");
  else
    trajectorycount = newnumsims;
}


void Options::setOutputInterval( int newoutputinterval )
{
  // TODO: error checking here? it can be negative, but not very.
  outputinterval = newoutputinterval;
  outputtime = 0.0;
}

void Options::setOutputTime( double newoutputtime )
{
  outputinterval = -1;
  outputtime = newoutputtime;
}

void Options::setGTWobble( int gtstate )
{
  gtenable = gtstate;
}

void Options::setLogML( int newlogml )
{
  logml = newlogml;
}

void Options::setEnergyMode( int mode )
{
  energymode = mode;
}

void Options::setInitialSeed( long seedval )
{
  initialseed = seedval;
}

void Options::setDangles( int dangleopt )
{
  dangles = dangleopt;
}

void Options::setTemperature( double temp )
{
  temperature = temp;
}

void Options::setParameterFile( char *newparamfile )
{
  if( strlen(newparamfile) > 79 )
    fprintf(stderr,"WARNING: Logfile entry too long. Please use a shorter name.\n");
  else
    {
      strcpy(energymodelfilename,newparamfile);
      free(newparamfile); // strdup'd in parser.
    }
  if( energymodel != 0 )
    fprintf(stderr,"WARNING: You used both energy model parameter options (#paramfiletype) and the automated energy model searching options (#energymodel).\n");
}

void Options::setTrajectoryfile( char *newtrajfilename )
{
  if( strlen(newtrajfilename) > 79 )
    fprintf(stderr,"WARNING: Trajectory filename entry too long. Please use a shorter name.\n");
  else
    {
      strcpy(trajectoryfilename,newtrajfilename);
      free(newtrajfilename); // strdup'd in parser.
    }
}

void Options::setParameterType( int type ) 
{
  energymodel_type = type;
  if( energymodel != 0 )
    fprintf(stderr,"WARNING: You used both energy model parameter options (#paramfiletype) and the automated energy model searching options (#energymodel).\n");
}

void Options::setIntramolecularScaling( double newscale )
{
  intramolecularscaling = newscale;
}

void Options::setIntermolecularScaling( double newscale )
{
  intermolecularscaling = newscale;
}


void Options::setEnergyModel( int modeltype )
{
  energymodel = modeltype;
  if(modeltype == VIENNADNA )
    energymodel_type = 0;
  else
    energymodel_type = 1;
}

void Options::setTrajectoryType( int ttype )
{
  trajtype = ttype;
}

/* 
 Not so simple mutators for interacting with the join rate between complexes.

 void Options::setJoinRate( double newjoinrate ) 
 void Options::setJoinRateByConcentration( double concentration )
 void Options::setJoinRateByVolume( double volume ) 
 void Options::setJoinRateByEnergy( double energy )
*/

void Options::setJoinRate( double newjoinrate ) 
{
  joinconc = newjoinrate * 55.6 / .001;
}

void Options::setJoinRateByConcentration( double concentration )
{
  // input concentration is in milli-molar (mM) units
  // future versions of the parser will convert to these units.
  //
  // To get the energy difference (and thus the join rate)
  // we must find the number of solvent molecules M_s in the volume
  // which would contain a single molecule at this concentration.
  // the energy term is then kT ln M_s + join term.

  // Roughly, this conversion is as follows:
  // Convert concentration to moles per cubic meter.
  // Convert concentration to molecules per cubic meter
  // invert to get cubic meters per molecule
  // multiply by water density to get grams water per molecule
  // divide by water molecular weight (g/mol) to get mols water per (strand) molecule
  // multiply by avogadro's number to get molecules water per (strand) molecule
  // This result is our M_s, the number of solvent molecules in our volume.

  // Shortened form, to try and avoid floating point error:
  // joinenergy = kT log( W / (.001 * C)) 
  // W = 55.6 mol/L (molarity of water)

  // joinrate = exp( -dG / kT )
  //          = exp( - kT log( W / (.001*C)) / kT)
  //          = exp( -log( W / (.001*C)))
  //          = .001 C / W

  //  joinenergy = log( 55.6 / (.001 * concentration));
  //joinrate = .001 * concentration / 55.6;
  joinconc = concentration;

  // Temperature (C) - water density (g/L)
  // 0  : .99984
  // 10 : .9970
  // 20 : .99821
  // 30 : .99565
  // 40 : .99222
  // 50 : .98803
  // 60 : .98320
  // 70 : .97778
  // 80 : .97182
  // 90 : .96535
  // 100: .95840
  
}

void Options::setJoinRateByVolume( double volume ) 
{

}


void Options::setJoinRateByEnergy( double energy )
{


}


/*
  Accessors for basic info. 

  double Options::getSimulationTime( void )
  int Options::getSimulationMode( void )
  double Options::getNumSimulations( void )
  double Options::getJoinRate( void )
 int Options::getStopCount( void );
  long Options::getInitialSeed( void );
  class stopcomplexes *Options::getStopComplexList( int index );
  double Options::getIntramolecularScaling( void );
  double Options::getIntermolecularScaling( void );
  double Options::getJoinEnergyExtra( void );
  int Options::getEnergyMode( void )

*/

int Options::getEnergyMode( void )
{
  return energymode;
}

int Options::getSimulationMode( void )
{
  return simulationmode;
}

/*double Options::getJoinEnergyExtra( void )
{
  return joinenergy;
}*/

double Options::getJoinConcentration( void )
{
  return joinconc;
}

double Options::getIntramolecularScaling( void )
{
  return intramolecularscaling;
}

double Options::getIntermolecularScaling( void )
{
  return intermolecularscaling;
}

class stopcomplexes *Options::getStopComplexList( int index )
{
  class stopcomplexes *traverse = stoplist;
  while( index > 0 && traverse != NULL )
    {
      index--;
      traverse = traverse->next;
    }
  if( traverse == NULL ) // ERROR!
    {
      fprintf(stderr,"ERROR: Options::getStopComplexList indexed past end of list.\n");
      exit(1);
    }
  return traverse;
}

int Options::getStopCount( void )
{
  return stopcount;
}

int Options::getStopOptions( void )
{
  return stopoptions;
}

long Options::getInitialSeed( void )
{
  return initialseed;
}

double Options::getSimulationTime( void )
{
  return simulationtime;
}

double Options::getTemperature( void  )
{
  return temperature;
}

/*double Options::getJoinRate( void )
{
  return joinrate;
}*/

int Options::getDangles( void )
{
  return dangles;
}

int Options::getGTWobble( void )
{
  return gtenable;
}

int Options::getLogML( void )
{
  return logml;
}

int Options::getNumSimulations( void )
{
  return trajectorycount;
}


int Options::getOutputInterval( void )
{
  return outputinterval;
}

double Options::getOutputTime( void )
{
  return outputtime;
}

int Options::getOutputState( void )
{
  if( outputinterval < 0 ) return 0;

  if( currentinterval == outputinterval )
    return 1;
  return 0;
}

int Options::getTrajectoryType( void )
{
  return trajtype;
}

/* 

   Output Related Functions

   void Options::printStatusLine( long r_seed, char *stop_id, double end_time )
   void Options::printStatusLine( long r_seed )
   void Options::printTrajLine( char *trigger_id, double cur_time )
   FILE *Options::getOutputDescriptor( void )
   void Options::closeOutput( void )

*/

void Options::printStatusLine( long r_seed, char *stop_id, double end_time )
{
  if( simulationmode & SIMULATION_PYTHON )
    {
      python_trajectory_time = end_time;
      python_trajectory_tag = stop_id;
      python_trajectory_completion_flag = 1;
    }
  else
    {
      if( logfile == NULL )
	{
	  logfile = fopen( logfilename, "at");
	  if( logfile == NULL ) // something bad happened
	    {
	      fprintf(stderr,"ERROR: Could not open logfile (%s) for writing.\n",logfilename);
	      exit(1);
	    }
	  
	}


      fprintf(logfile, "(0x%07x) %s %lf\n",r_seed, stop_id, end_time );
      fflush(logfile);
    }

}

void Options::printStatusLine_First_Bimolecular( long r_seed, int completiontype, double completiontime, double forwardrate, char *tag )
{
  if( simulationmode & SIMULATION_PYTHON )
    {
      python_trajectory_time = completiontime;
      python_k_collision = forwardrate;
      python_current_seed = r_seed;
     
      if( completiontype == STOPCONDITION_FORWARD)
	python_trajectory_tag = "FORWARD:"+ (std::string) (const char *)tag;
      else if ( completiontype == STOPCONDITION_REVERSE )
	python_trajectory_tag = "REVERSE";
      else if (completiontype == STOPCONDITION_TIME )
	python_trajectory_tag = "TIME";
      else if ( completiontype == STOPCONDITION_ERROR )
	python_trajectory_tag = "ERROR";

      python_trajectory_completion_flag = 1;
    }
  else
    {
      if( logfile == NULL )
	{
	  logfile = fopen( logfilename, "at");
	  if( logfile == NULL ) // something bad happened
	    {
	      fprintf(stderr,"ERROR: Could not open logfile (%s) for writing.\n",logfilename);
	      exit(1);
	    }
	}
  

      if( completiontype == STOPCONDITION_FORWARD)
	fprintf(logfile, "(0x%07x) Rate: %lf FORWARD:%s %le\n",r_seed, forwardrate, tag, completiontime );
      else if ( completiontype == STOPCONDITION_REVERSE )
	fprintf(logfile, "(0x%07x) Rate: %lf REVERSE %le\n",r_seed, forwardrate, completiontime );
      else if (completiontype == STOPCONDITION_TIME )
	fprintf(logfile, "(0x%07x) Rate: %lf TIME %le\n",r_seed, forwardrate, completiontime );
      else if( completiontype == STOPCONDITION_ERROR )
	fprintf(logfile, "(0x%07x) Rate: %lf ERROR %le\n", r_seed, forwardrate, completiontime );
      
      
      fflush(logfile);
    }
}

void Options::printStatusLine_Warning( int warningtype, int data )
{
  if( logfile == NULL )
    {
      logfile = fopen( logfilename, "at");
      if( logfile == NULL ) // something bad happened
	{
	  fprintf(stderr,"ERROR: Could not open logfile (%s) for writing.\n",logfilename);
	  exit(1);
	}

    }


  if( warningtype == 0 )
    fprintf(logfile, "WARNING: A total of %d trajectories reached the maximum allowed time! Forward rate will be measured incorrectly!\n", data ); 
    

}

void Options::printStatusLine_Final_First_Bimolecular( double total_rate, double *total_times, long *total_types, long total_count, double *computed_rate_means, double *computed_rate_mean_diff_squared )
{
  if( logfile == NULL )
    {
      logfile = fopen( logfilename, "at");
      if( logfile == NULL ) // something bad happened
	{
	  fprintf(stderr,"ERROR: Could not open logfile (%s) for writing.\n",logfilename);
	  exit(1);
	}

    }

  fprintf(logfile, "Simulation Complete: %ld trajectories total.\n", total_count);


  fprintf(logfile, "Collision Rate (k_coll): %lf (/M/s)\n", (total_rate / (double) total_count  ) );
  fprintf(logfile, "Collision Rate (k_coll) mean, std: %lf, %lf (/M/s)\n", computed_rate_means[2], sqrt(computed_rate_mean_diff_squared[2]));

  fprintf(logfile, "Reactive Collision Rate (k_1): %lf (/M/s)\n", (total_types[0] / (double) total_count ) * (total_rate / (double) total_count  ) );
  fprintf(logfile, "Reactive Collision Rate (k_1) mean, std: %lf, %lf (/M/s)\n", computed_rate_means[2] * (total_types[0] / (double) total_count ), sqrt(computed_rate_mean_diff_squared[2] * (total_types[0] / (double) total_count ) * (total_types[0] / (double) total_count)));

  fprintf(logfile, "Failed Collision Rate (k_1'): %lf (/M/s)\n", (total_types[1] / (double) total_count ) * (total_rate / (double) total_count  ) );
  fprintf(logfile, "Failed Collision Rate (k_1') mean, std: %lf, %lf (/M/s)\n", computed_rate_means[2] * (total_types[1] / (double) total_count ), sqrt(computed_rate_mean_diff_squared[2] * (total_types[1] / (double) total_count ) * (total_types[1] / (double) total_count)));

  /// (double) total_types[0])  );
  if( total_types[0] != 0 && total_types[0] != 1 )
    {
      fprintf(logfile, "Reactive Trajectory Rate (k_2) mean, std: %lf, %lf (/s)\n", computed_rate_means[0], sqrt(computed_rate_mean_diff_squared[0] / (double) (total_types[0] - 1 )));
      fprintf(logfile, "Reactive Trajectories: %ld, Average Time: %le\n", total_types[0], total_times[0] / total_types[0] );
    }
  else if( total_types[0] != 0 )
    {
      fprintf(logfile, "Reactive Trajectory Rate (k_2) mean: %lf (/s)\n", computed_rate_means[0]);
      fprintf(logfile, "Reactive Trajectories: %ld, Average Time: %le\n", total_types[0], total_times[0] / total_types[0] );
    }
  else
    {
        fprintf(logfile, "Reactive Trajectory Rate (k_2) mean: N/A (/s)\n");
	fprintf(logfile, "Reactive Trajectories: %ld, Average Time: N/A\n", total_types[0] );
    }

  if( total_types[1] != 0 && total_types[1] != 1)
    {
      fprintf(logfile, "Failed Trajectory Rate (k_2') mean, std: %lf, %lf (/s)\n", computed_rate_means[1], sqrt(computed_rate_mean_diff_squared[1] / (double) (total_types[1] - 1)));
      fprintf(logfile, "Failed Trajectories: %ld, Average Time: %le\n", total_types[1], total_times[1] / total_types[1] );
    }
  else
    {
      fprintf(logfile, "Failed Trajectory Rate (k_2') mean: N/A (/s)\n");
      fprintf(logfile, "Failed Trajectories: %ld, Average Time: N/A\n", total_types[1] );
    }

  fflush(logfile);
}


void Options::printStatusLine( long r_seed )
{
  if( logfile == NULL )
    {
      logfile = fopen( logfilename, "at");
      if( logfile == NULL ) // something bad happened
	{
	  fprintf(stderr,"ERROR: Could not open logfile (%s) for writing.\n",logfilename);
	  exit(1);
	}

    }

  fprintf(logfile, "(0x%0x) No Stop Found\n",r_seed );
  fflush(logfile);

}


FILE *Options::getOutputDescriptor( void )
{
  return logfile;
}

void Options::closeOutput( void )
{
  if( logfile != NULL )
    fclose( logfile );
  if( trajfile != NULL )
    fclose( trajfile );
}

void Options::printTrajLine( char *trigger_id, double cur_time )
{
  // WARNING: static member!!
  static char lasttag[80]="\0";
  if( trajfile == NULL )
    {
      trajfile = fopen( trajectoryfilename, "at");
      if( trajfile == NULL ) // something bad happened
	{
	  fprintf(stderr,"ERROR: Could not open trajectory file (%s) for writing.\n",trajectoryfilename);
	  exit(1);
	}

    }

  if( trigger_id == NULL )
    {
      fprintf(trajfile, "Trajectory %d:\n",(int) cur_time);
      strcpy(lasttag,"");
    }
  else if( strcmp(trigger_id,lasttag) != 0 )
    {
      fprintf(trajfile, "%s : %lf\n", trigger_id, cur_time );
      strcpy(lasttag,trigger_id);
    }
  fflush(trajfile);
}


/*
  char *Options::getParameterFile( void )
  int Options::getParameterType( void )
*/

char *Options::getParameterFile( void )
{
  return energymodelfilename;
}

int Options::getParameterType( void )
{
  return energymodel_type;
}

int Options::getEnergyModel( void )
{
  return energymodel;
}

/*
  char *Options::getSequence( int complex_id );
  char *Options::getStructure( int complex_id );
  char *Options::getSequenceByID( char *id );
  char *Options::getID( int complex_id );
*/

char *Options::getSequence( int complex_id )
{
  char *sequence = NULL;
  class identlist *traverse = NULL;

  assert( complex_id >= 0 && complex_id <= start_count );

  if( complex_id == start_count) return NULL;

  traverse = start_list[complex_id]->strand_ids;

  sequence = new char[strlen(start_list[complex_id]->structure)+1];
  strcpy( sequence, getSequenceByID( traverse->id ) );
  traverse = traverse->next;

  while( traverse != NULL )
    {
      strcat( sequence, "_" );
      strcat( sequence, getSequenceByID( traverse->id ));
      traverse = traverse->next;  
    }
  return sequence;
}

char *Options::getBoltzmannStructure( int complex_id )
{
  return getStructure( complex_id ); // for now.
}

char *Options::getStructure( int complex_id )
{
  char *structure = NULL;

  assert( complex_id >= 0 && complex_id <= start_count );

  if( complex_id == start_count) return NULL;

  structure = new char[strlen(start_list[complex_id]->structure)+1];
  strcpy( structure, start_list[complex_id]->structure );

  return structure;
}

class identlist *Options::getID_list( int complex_id )
{
  assert( complex_id >= 0 && complex_id <= start_count );
  if( complex_id == start_count) return NULL;
  return start_list[complex_id]->strand_ids;
}

char *Options::getSequenceByID( char *id )
{
  class strandlist *traverse = strands;
  while( traverse != NULL )
    {
      if( strcmp(id, traverse->id) == 0)
	return traverse->seq;
      traverse = traverse->next;
    }

  fprintf(stderr, "WARNING: Strand ID used did not match any listed as strands.\n");
  return NULL;
}

void printError( int location, char **argv );


/*

   int Options::readCommandLine( int argc, char **argv )


   Handles the command line options, most importantly the input file.
   If others are missing, that's ok, but input file is the killer.

   Currently only handles input file, will add the appropriate overrides
   later.

*/
int Options::readCommandLine( int argc, char **argv )
{
  int index;
  char *temp;

  for( index = 1; index <= argc; index ++)  // loop over all arguments.
    {
      if( argv[index][0] != '-' || argv[index][1] != '-')
	{
	  printError( index, argv );
	  return -1;
	}

      /* Case blocks */
      if( strstr(argv[index],"inputfile") == &argv[index][2] )
	{
	  // Input filename.
	  temp = strchr( argv[index], '=' );
	  if( temp == NULL )
	    {
	      printError( index, argv );
	      return -1;
	    }
	  strcpy(inputfilename, temp+1);
	}
      else
	{
	  fprintf(stderr,"Unknown argument: %s\n",argv[index]);
	}
    }
  return 1;
}


/*

   printError( int location, char **argv )
   
   Auxilary function for printing an error message on malformed input.

*/
void printError( int location, char **argv)
{
  fprintf(stderr,"Argument #%d invalid: %s.\n", location, argv[location]);
  fprintf(stderr,"Usage info would go here once I type it.\n");

}


/*

   int Options::loadInputFile( void )

   Read in the options from the input file (which should be in inputfilename).

   TODO: for this function, add sanity checks on parameter sizes. IE, no negative or 0 entries for simulation time, etc.
*/

int Options::loadInputFile( void )
{
  char buffer[MAX_BUFFER];
  int buf_len;
  FILE *fp = NULL;
  char *temp;

  fp = fopen(inputfilename,"rt");
  if( fp == NULL )
    {
      fprintf(stderr,"Could not open input file: %s\n",inputfilename);
      return -1;
    }

  while( !feof(fp) )
    {
      fgets( buffer, MAX_BUFFER-1, fp);
      buffer[1023]='\0';
      buf_len = strlen(buffer);
      if( buffer[0] != '#') continue;                 // malformed line
      if( buf_len > 1 && buffer[1] == '#' ) continue; // this line is a comment

      // Sanitize the input buffer. Removes whitespace, drops caps.
      sanitize(buffer,&buf_len);

      // Parse single-line option block.

      if( !(flagarray & OPTION_ENERGYMODEL)
	  && strstr(buffer,"energymodel") == &buffer[1] )
	{
	  flagarray |= OPTION_ENERGYMODEL;
	  // Input filename.
	  temp = strchr( buffer, '=' );
	  if( temp == NULL )
	    {
	      fprintf(stderr,"Option #energymodel requires a parameter.\n");
	      return -1;
	    }
	  strcpy(energymodelfilename, temp+1);
	}
      else if( !(flagarray & OPTION_LOGFILE)
	       && strstr(buffer,"logfile") == &buffer[1] )
	{
	  flagarray |= OPTION_LOGFILE;
	  // Input filename.
	  temp = strchr( buffer, '=' );
	  if( temp == NULL )
	    {
	      fprintf(stderr,"Option #logfile requires a parameter.\n");
	      return -1;
	    }
	  strcpy(logfilename, temp+1);
	}

      else if( !(flagarray & OPTION_STRANDCOUNT)
	       && strstr(buffer,"strandcount") == &buffer[1] )
	{
	  flagarray |= OPTION_STRANDCOUNT;
	  // Input filename.
	  temp = strchr( buffer, '=' );
	  if( temp == NULL )
	    {
	      fprintf(stderr,"Option #strandcount requires a parameter.\n");
	      return -1;
	    }
	  strandcount = atoi(temp+1);
	}
      else if( !(flagarray & OPTION_STRANDLIST)
	       && strstr(buffer,"strandlist") == &buffer[1] )
	{
	  flagarray |= OPTION_STRANDLIST;
	  if( !(flagarray & OPTION_STRANDCOUNT ))
	    {
	      fprintf(stderr, "Must specify number of strands before giving the list of strands!\n");
	      return -1;
	    }
	  readStrandList( fp, buffer );
	}


      else if( !(flagarray & OPTION_COMPLEXCOUNT)
	       && strstr(buffer,"complexcount") == &buffer[1] )
	{
	  flagarray |= OPTION_COMPLEXCOUNT;
	  // Input filename.
	  temp = strchr( buffer, '=' );
	  if( temp == NULL )
	    {
	      fprintf(stderr,"Option #complexcount requires a parameter.\n");
	      return -1;
	    }
	  //
	}
      else if( !(flagarray & OPTION_COMPLEXLIST)
	       && strstr(buffer,"complexlist") == &buffer[1] )
	{
	  flagarray |= OPTION_COMPLEXLIST;
	  if( !(flagarray & OPTION_COMPLEXCOUNT ))
	    {
	      fprintf(stderr, "Must specify number of complexes before giving the list of complexes!\n");
	      return -1;
	    }
	  readComplexList( fp, buffer );
	}
      else if( !(flagarray & OPTION_SIMULATIONTIME)
	       && ((strstr(buffer,"simtime") == &buffer[1] ) ||  (strstr(buffer,"simulationtime") == &buffer[1])))
	{
	  flagarray |= OPTION_SIMULATIONTIME;
	  // Input filename.
	  temp = strchr( buffer, '=' );
	  if( temp == NULL )
	    {
	      fprintf(stderr,"Option #simtime requires a parameter.\n");
	      return -1;
	    }
	  simulationtime = atof(temp+1);
	}
      else if( !(flagarray & OPTION_TRAJECTORYCOUNT)
	       && ((strstr(buffer,"numsims") == &buffer[1] ) ||  (strstr(buffer,"trajectorycount") == &buffer[1])))
	{
	  flagarray |= OPTION_TRAJECTORYCOUNT;
	  // Input filename.
	  temp = strchr( buffer, '=' );
	  if( temp == NULL )
	    {
	      fprintf(stderr,"Option #numsims requires a parameter.\n");
	      return -1;
	    }
	  trajectorycount = atoi(temp+1);
	}



      
    }

}

void Options::finalizeInput( void )
{
  class complex_item *traverse = start_structure;
  int location = 0;

  if( traverse == NULL ) 
    {
      fprintf(stderr,"ERROR: no start structure given.\n");
      exit(1);
    }

  start_count = 1;
  while( traverse != NULL )
    {
      traverse = traverse->next;
      if( traverse != NULL )
	start_count++;
    }
  start_list = new complex_item *[start_count];
  traverse = start_structure;
  while( traverse != NULL )
    {
      start_list[location] = traverse;
      traverse = traverse->next;
      location++;
    }

}


void Options::sanitize( char *buf, int *length )
{
  ;
}


void Options::readStrandList( FILE *fp, char *buffer )
{
  ;
}


void Options::readComplexList( FILE *fp, char *buffer )
{
  ;
}



/*

  Options
      Python Interface functions

 
  int addSequence_Python(  std::string id, std::string seq );
    //... plus other versions for later implementation
  
  // list of strands interface
  int beginStrandList_Python( void );
  int addStrandtoStrandList_Python( type name);  // define order as right to left
  int clearStrandList_Python( void );

  int beginStartStructure_Python( void ); 
  int beginStartStructureComplex_Python( void );
    // uses implied strand list already created
  int addStartStructureComplexStructure_Python( type struc );
  int addStartStructureComplexType_Python( type type );
  int completeStartStructureComplex_Python( void );
  int completeStartStructure_Python( void );

  int beginStopStructure_Python( type tag ); 
  int beginStopStructureComplex_Python( void );
    // uses implied strand list already created
  int addStopStructureComplexStructure_Python( type struc );
  int addStopStructureComplexType_Python( type type );
  int completeStopStructureComplex_Python( void );
  int completeStopStructure_Python( void );

  void setLogfile_Python( type newlogfilename );
  void setTrajectoryfile_Python( type newtrajfilename );
  void setParameterFile_Python( type newparamfile );

*/

/* 
   int Options::addSequence_Python( std::string id, std::string seq )

   Adds a sequence to our list of sequences, with id and seq as given in parameters. Returns 1 on success, 0 on failure.

   Only failure mode is allocation failure.
*/

int Options::addSequence_Python( std::string id, std::string seq)
{
  char *id2 = new char[id.length()+1];
  char *seq2 = new char[seq.length()+1];

  if( id2 == NULL || seq2 == NULL )
    return 0;
  strcpy( id2, id.c_str());
  strcpy(seq2, seq.c_str());

 if( strands == NULL )
    strands = new strandlist( id2, seq2 );
  else
    strands = new strandlist( id2, seq2, strands );
  
 if( strands == NULL )
    return 0; // memory allocation error
  else
    return 1;
}


/*
  Options
       Python Interface

  // list of idents interface
  int addStrandtoIdentList_Python( std::string name);  // define order as right to left
  int clearIdentList_Python( void );


*/


/*
  int Options::addStrandToIdentList_Python( std::string name );

  adds name to the current python_identlist, at the left edge (beginning edge)

  Returns 1 on success, 0 on failure (memory allocation failure only).

*/

int Options::addStrandToIdentList_Python( std::string name )
{
  char *name2 = new char[ name.length() + 1 ];
  if( name2 == NULL )
    return 0;
  strcpy( name2, name.c_str() );

  if( python_identlist == NULL )
    python_identlist = new identlist( name2 );
  else
    python_identlist = new identlist( name2, python_identlist );

  if( python_identlist == NULL )
    return 0;
  else
    return 1;
}


/*
  int Options::clearIdentList_Python( void );
  
  returns 1 if python_identlist was cleared, 0 if there was no current python_identlist.

*/

int Options::clearIdentList_Python( void )
{
  if( python_identlist == NULL )
    return 0;
  else
    {
      // delete python_identlist;
      python_identlist = NULL;
      return 1;
    }
}

/* 
   int Options::beginStartStructure_Python( void )

   returns 1 always.

*/


int Options::beginStartStructure_Python( void )
{
  if( python_start_structure != NULL )
    delete python_start_structure;
  return 1;
}

/*
  int addStartStructureComplexStructure_Python( std::string struc )

*/

int Options::addStartStructureComplex_Python( std::string struc, int type )
{
  char *struc2 = new char[ struc.length() + 1 ];
  if( struc2 == NULL )
    return 0;
  strcpy( struc2, struc.c_str() );

  python_start_structure = new complex_item( struc2, python_identlist, python_start_structure );

  if( python_start_structure->structure == NULL )
    return 0;
  if( python_start_structure->strand_ids == NULL )
    return 0;
  else
    return 1;  
}

/*
  int Options::finalizeStartStructure_Python( void )

  uses the current start struc to be the system's start structure.

*/

int Options::finalizeStartStructure_Python( void )
{
  if( start_structure != NULL )
    return 0;
 
  return addStartStructure( python_start_structure,1 );
}



/* 
   int Options::beginStopStructure_Python( void )

   returns 0 on allocation failure; 1 otherwise

*/


int Options::beginStopStructure_Python( std::string tag )
{
  char *tag2 = new char[ tag.length() + 1];
  if( tag2 == NULL )
    return 0;
  strcpy( tag2, tag.c_str() );

  python_stop_structure = new stopcomplexes( tag2, NULL, python_stop_structure );
  if( python_stop_structure == NULL )
    return 0;

  return 1;
}

/* 
   int Options::beginStopStructureComplex_Python

   returns 0 for memory allocation failure, 1 otherwise (success)

*/

int Options::addStopStructureComplex_Python( std::string struc, int type, int count )
{
  char *struc2 = new char[ struc.length() + 1 ];
  if( struc2 == NULL )
    return 0;
  strcpy( struc2, struc.c_str() );

  if( python_stop_structure == NULL )
    return 0;

  python_stop_structure->citem = new complex_item( struc2, python_identlist, python_stop_structure->citem, type, count );

  if( python_stop_structure == NULL || python_stop_structure->citem == NULL )
    return 0;
  if( python_stop_structure->citem->structure == NULL )
    return 0;
  if( python_stop_structure->citem->strand_ids == NULL )
    return 0;
  else
    return 1;

}


/*
  int Options::finalizeStopStructures_Python( void )

  sets the stop structures to be the previously completed list

*/

int Options::finalizeStopStructures_Python( void )
{
  if( python_stop_structure == NULL )
    return 0;
  return addStopStructures( python_stop_structure, 1 );
}



/*
  void Options::setLogfile_Python( std::string newlogfilename )
  

*/

void Options::setLogfile_Python( std::string newlogfilename )
{
  strncpy( logfilename, newlogfilename.c_str(), 79 );
  logfilename[79] = '\0';
}

/*
  void Options::setTrajectoryfile_Python( std::string newtrajfilename )

*/

void Options::setTrajectoryfile_Python( std::string newtrajfilename )
{
  strncpy( trajectoryfilename, newtrajfilename.c_str(), 79 );
  trajectoryfilename[79] = '\0';

}

/*
    void setParameterFile_Python( std::string newparamfile );

*/

void Options::setParameterFile_Python( std::string newparamfile )
{
  strncpy( energymodelfilename, newparamfile.c_str(), 79 );
  energymodelfilename[79] = '\0';
}


/* Python Interface

   simulation time and output methods
   
   setCurSimTime
  double checkCurrentTime_Python(void ); 
  int checkCompleted_Python( void );
  std::string getTrajectoryTag_Python( void );
  double getTrajectoryTime_Python( void );


*/


/*
  void Options::setCurSimTime( double cur_time )

*/

void Options::setCurSimTime( double cur_time )
{
  python_current_time = cur_time;
}

/*
  void Options::setRateMethod( int method)

*/

void Options::setRateMethod( int method )
{
  ratemethod = method;
}

/*
  int Options::getRateMethod( void )

*/

int Options::getRateMethod( void )
{
  return  ratemethod;
}

/*
  double Options::checkCurrentTime_Python( void )
*/

double Options::checkCurrentTime_Python( void )
{
  return python_current_time;
}

/* 
   int Options::checkCompleted_Python( void )
*/

int Options::checkCompleted_Python( void )
{
  return python_trajectory_completion_flag;
}

/* 
   int Options::resetCompleted( void )
*/

void Options::resetCompleted_Python( void )
{
  python_trajectory_completion_flag = 0;
}

/*
  std::string Options::getTrajectoryTag_Python( void )
*/

std::string Options::getTrajectoryTag_Python( void )
{
  return python_trajectory_tag;
}

/*
  double Options::getTrajectoryTime_Python( void )
*/

double Options::getTrajectoryTime_Python( void )
{
  return python_trajectory_time;
}

/* 
   void Options::haltCurrentTrajectory_Python( void )
*/

void Options::haltCurrentTrajectory_Python( void )
{
  python_halt_trajectory_flag = 1;
  python_suspend_trajectory_flag = 0;
}

/*
  int Options::getCurrentTrajectoryHaltFlag( void )
*/

int Options::getCurrentTrajectoryHaltFlag( void )
{
  if( python_halt_trajectory_flag)
    {
      python_halt_trajectory_flag = 0;
      return 1;
    }
  else
    return 0;
}

/* 
   void Options::suspendCurrentTrajectory_Python( void )
*/

void Options::suspendCurrentTrajectory_Python( void )
{
  python_suspend_trajectory_flag = 1;
}

/* 
   void Options::resumeCurrentTrajectory_Python( void )
*/

void Options::resumeCurrentTrajectory_Python( void )
{
  python_suspend_trajectory_flag = 0;
}

/*
  int Options::getCurrentTrajectorySuspendFlag( void )
*/

int Options::getCurrentTrajectorySuspendFlag( void )
{
  return python_suspend_trajectory_flag;
}

/*
  double Options::getCollisionRate_Python( void )
*/

double Options::getCollisionRate_Python( void )
{
  return python_k_collision;
}

/*
  void Options::setCollisionRate_Python( double newrate )
*/

void Options::setCollisionRate_Python( double newrate )
{
  python_k_collision = newrate;
}


/*
  long Options::getCurrentSeed_Python( void )
*/

long Options::getCurrentSeed_Python( void )
{
  return python_current_seed;
}

