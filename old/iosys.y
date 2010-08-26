
%error-verbose
%{

/*
   Copyright (c) 2007-2008 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../options.h"

int yydebug = 1;
int stopstartflag = 0;
 int curLine = 1;
//LogicNode *rootNode;
extern Options *GlobalOptions;

 extern int yywrap();
 extern int yylex( void );
 int yyparse( void );

 //#define YYERROR_VERBOSE


int yywrap()
{
	return 1;
}

void yyerror( const char *str)
{
  fprintf(stderr,"error: %s near line %d\n",str,curLine );//, yygetlineno());//yylloc.first_line);
}
 
//main()
//{
//	yyparse();
//}

int isloosestructure( char *struc )
{
  if( strchr( struc , (int) '*') != NULL )
    return 1;
  else return 0;
}

int countpercent( char *struc, double perc )
{
  int count = 0;
  int loop;
  for( loop = 0; loop < strlen(struc); loop++ )
    {
      if( struc[loop] == '(' || struc[loop] == ')' || struc[loop] == '.' )
	count++;
    }
  return (int) floor( (double) count * perc / 100);
}

%}

%union
{
  int num;
  double fnum;
  char *buffer;
  class identlist *listid;
  class complex_item *citem;
  class stopcomplexes *stoplist;
  //  class LogicNode *node;  
}

%token <num> NUMBER 
%token <fnum> FLOAT
%token <buffer> ID SEQUENCE STRUCTURE
%token TOKPERCENT TOKSTARTSTRUCTURE TOKSTOPSTRUCTURE TOKCOMMA TOKSTRANDS TOKNEWLINE TOKTAGS TOKSIMTIME TOKNUMSIMS TOKOUTPUTINTERVAL TOKSTOPOPTIONS TOKLOGFILE TOKEQUALS TOKJOINENERGY TOKCONCENTRATION TOKVOLUME TOKJOINRATE TOKPARAMFILE TOKPARAMTYPE TOKSTARTSEED TOKTRAJFILE TOKTRAJTYPE TOKENERGYMODEL TOKVIENNADNA TOKNUPACKDNA23 TOKNUPACKRNA23 TOKTEMPERATURE TOKOUTPUTTIME TOKGTWOBBLE TOKLOGML STRING TOKINTRA TOKINTER TOKDANGLES TOKBOUND TOKDISASSOC TOKSIMMODE TOKRATEMETHOD
%type <listid> idlist
%type <buffer> tagdata
%type <citem> complex complexlist
%type <stoplist> stopstructurelist


%%
infile:  
        statement extra_newlines
        |
        statement 
        |
        statement infile
        |
        statement extra_newlines infile
        ;

extra_newlines:
        TOKNEWLINE
        |
        extra_newlines TOKNEWLINE
        ;

statement:
        stranddefs
	|
        stopstructure
        |
	startstructure
	|
	singlelineoptions
	|
	error TOKNEWLINE
{ fprintf(stderr,"ERROR: Bad statement name on line %d.\n",curLine - 1);}
;

singlelineoptions:
	simulationtime
	|
	simulationmode
	|
	numsimulations
	|
	outputinterval
        |
        outputtime
	|
	stopoptions
	|
	logfilename
	|
	paramfile
	|
	paramtype
	|
	startseed
	|
	joinrateoption
	|
	trajfilename
	|
	trajtype
	|
	temperature
        |
	energymodel
	|
	dangles
	|
	gtwobble
	|
	logml
	|
	intrascale
	|
	interscale
	|
	ratemethod
	;

stranddefs:
	TOKSTRANDS strandinfo
        |
        TOKSTRANDS error TOKNEWLINE
{fprintf(stderr,"ERROR: malformed #Strands block, near line %d\n",curLine);}
	;


strandinfo:
	strand
	|
	strandinfo strand
	;

strand:
	ID TOKCOMMA SEQUENCE TOKNEWLINE
{ /*printf("Sequence added: %s:%s\n",$1,$3);*/ GlobalOptions->addSequence( $1, $3 );}
	|
	ID TOKCOMMA SEQUENCE TOKCOMMA FLOAT TOKNEWLINE
{/* printf("Sequence added: %s:%s %lf\n",$1,$3,$5);*/ GlobalOptions->addSequence( $1,$3,$5 ); }
	|
	ID TOKCOMMA SEQUENCE TOKCOMMA NUMBER TOKNEWLINE
{/* printf("Sequence added: %s:%s %d\n",$1,$3,$5); */GlobalOptions->addSequence( $1,$3,$5 ); }
        |
	SEQUENCE TOKCOMMA SEQUENCE TOKNEWLINE
{/* printf("Sequence added: %s:%s\n",$1,$3);*/ GlobalOptions->addSequence( $1, $3);}
	;

startstructure:
	TOKSTARTSTRUCTURE complexlist
{ GlobalOptions->addStartStructure( $2 ); /*printf("Startstructure Parsed\n");*/ }
        |
        TOKSTARTSTRUCTURE error TOKNEWLINE
{fprintf(stderr,"ERROR: malformed #StartStructure block, near line %d\n",curLine);}
        ;

stopstructure:
	TOKSTOPSTRUCTURE stopstructurelist
{ GlobalOptions->addStopStructures( $2 ); /*printf("Stopstructure Parsed\n"); */ }
        |
        TOKSTOPSTRUCTURE error TOKNEWLINE
{fprintf(stderr,"ERROR: malformed #StopStructure block, near line %d\n",curLine-1);}
	;

complexlist:
	complex complexlist
{ $1->next = $2; $$ = $1; }
	|
	complex
{ $$ = $1; }
	;

stopstructurelist:
        complexlist TOKTAGS tagdata TOKNEWLINE stopstructurelist
{ $$ = new stopcomplexes($3, $1, $5); }
	|
	complexlist TOKTAGS tagdata TOKNEWLINE
{ $$ = new stopcomplexes($3, $1);  }
;


complex:
        idlist TOKNEWLINE STRUCTURE TOKNEWLINE
	{ if(isloosestructure( $3))
	    $$ = new complex_item( $3, $1, NULL, STOPTYPE_LOOSE_STRUCTURE);
	  else
	    $$ = new complex_item( $3, $1, NULL, STOPTYPE_STRUCTURE); 
	}
        |
        idlist TOKNEWLINE STRUCTURE TOKCOMMA FLOAT TOKPERCENT TOKNEWLINE
{ 
  $$ = new complex_item( $3, $1, NULL, STOPTYPE_PERCENT_OR_COUNT_STRUCTURE, countpercent( $3, $5 ));
}
        |
        idlist TOKNEWLINE STRUCTURE TOKCOMMA NUMBER TOKPERCENT TOKNEWLINE
{ 
  $$ = new complex_item( $3, $1, NULL, STOPTYPE_PERCENT_OR_COUNT_STRUCTURE, countpercent($3, (double) $5 ));
}
        |
        idlist TOKNEWLINE STRUCTURE TOKCOMMA FLOAT TOKNEWLINE
{ 
  $$ = new complex_item( $3, $1, NULL, STOPTYPE_PERCENT_OR_COUNT_STRUCTURE, (int) floor($5));
}
        |
        idlist TOKNEWLINE STRUCTURE TOKCOMMA NUMBER TOKNEWLINE
{ 
  $$ = new complex_item( $3, $1, NULL, STOPTYPE_PERCENT_OR_COUNT_STRUCTURE, $5 * 2);
}
        |
	idlist TOKNEWLINE TOKBOUND TOKNEWLINE
{ $$ = new complex_item( NULL, $1, NULL, STOPTYPE_BOUND); }
        |
        idlist TOKNEWLINE TOKDISASSOC TOKNEWLINE
{ $$ = new complex_item( NULL, $1, NULL, STOPTYPE_DISASSOC);}
        ;

idlist:
	ID TOKCOMMA idlist
{$$ = new identlist( $1, $3 ); /* printf("ID: %s\n",$1);*/}
	|
	ID
	{$$ = new identlist($1); /*printf("ID: %s\n",$1);*/}
	;

tagdata:
	ID
{ $$ = $1; /*printf("ID: %s\n",$1);*/}
;

simulationmode:
        TOKSIMMODE TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setSimulationMode( (int) $3 ); }
        |
	TOKSIMMODE TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setSimulationMode( (int) $3 ); }
        ;

simulationtime:
	TOKSIMTIME TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setSimulationTime( (double) $3 ); } 
	|
	TOKSIMTIME TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setSimulationTime( $3 );} 
	;

intrascale:
	TOKINTRA TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setIntramolecularScaling( $3 ); }
        |
	TOKINTRA TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setIntramolecularScaling( (double) $3 ); }
        ;

interscale:
	TOKINTER TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setIntermolecularScaling( $3 ); }
        |
	TOKINTER TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setIntermolecularScaling( (double) $3 ); }
        ;

temperature:
	TOKTEMPERATURE TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setTemperature( $3 );}
        |
	TOKTEMPERATURE TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setTemperature((double) $3);}
        ;

dangles:
	TOKDANGLES TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setDangles( (int) $3 ); }
        |
	TOKDANGLES TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setDangles( $3 ); }
        ;

gtwobble:
        TOKGTWOBBLE TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setGTWobble( $3 );}
        |
        TOKGTWOBBLE TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setGTWobble( (int) $3 );}
        ;

logml:
	TOKLOGML TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setLogML( $3 ); }
        |
	TOKLOGML TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setLogML( (int) $3); }
        ;

startseed:
	TOKSTARTSEED TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setInitialSeed( $3 );}
 
numsimulations:
	TOKNUMSIMS TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setNumSimulations( $3 );} 
	|
	TOKNUMSIMS TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setNumSimulations( (int) $3 );} 
	;

 
outputinterval:
	TOKOUTPUTINTERVAL TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setOutputInterval( $3 );} 
	|
	TOKOUTPUTINTERVAL TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setOutputInterval( (int) $3 );} 
	;

outputtime:
	TOKOUTPUTTIME TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setOutputTime( $3 );} 
	|
	TOKOUTPUTTIME TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setOutputTime( (int) $3 );} 
	;

ratemethod:
        TOKRATEMETHOD TOKEQUALS NUMBER TOKNEWLINE
	{GlobalOptions->setRateMethod( $3);}
        |
        TOKRATEMETHOD TOKEQUALS FLOAT TOKNEWLINE
	{GlobalOptions->setRateMethod( (int) $3 );}
	;

stopoptions:
	TOKSTOPOPTIONS TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setStopOptions( $3 );} 
	|
	TOKSTOPOPTIONS TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setStopOptions( (int) $3 );} 
	;
     
logfilename:
	TOKLOGFILE TOKEQUALS ID TOKNEWLINE
{ GlobalOptions->setLogfile( $3 );} 
        ;

trajfilename:
	TOKTRAJFILE TOKEQUALS ID TOKNEWLINE
{ GlobalOptions->setTrajectoryfile( $3 );} 
        ;

trajtype:
	TOKTRAJTYPE TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setTrajectoryType( $3 );}

paramfile:
	TOKPARAMFILE TOKEQUALS ID TOKNEWLINE
{ GlobalOptions->setParameterFile( $3 ); }
        ;

paramtype:
	TOKPARAMTYPE TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setParameterType( $3 ); }
        ;

energymodel:
	TOKENERGYMODEL TOKEQUALS TOKVIENNADNA TOKNEWLINE
{ GlobalOptions->setEnergyModel( VIENNADNA );}
        |
	TOKENERGYMODEL TOKEQUALS TOKNUPACKDNA23 TOKNEWLINE
{ GlobalOptions->setEnergyModel( NUPACKDNA23 );}
        |
	TOKENERGYMODEL TOKEQUALS TOKNUPACKRNA23 TOKNEWLINE
{ GlobalOptions->setEnergyModel( NUPACKRNA23 );}
	;

joinrateoption:
	TOKJOINRATE TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setJoinRate( (double) $3); }
        |
	TOKJOINRATE TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setJoinRate( $3 ); }
        |
	TOKCONCENTRATION TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setJoinRateByConcentration( (double) $3 ); }
        |
	TOKCONCENTRATION TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setJoinRateByConcentration( $3 ); }
        |
	TOKVOLUME TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setJoinRateByVolume( (double) $3 ); }
        |
	TOKVOLUME TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setJoinRateByVolume( $3 ); }
        |
	TOKJOINENERGY TOKEQUALS NUMBER TOKNEWLINE
{ GlobalOptions->setJoinRateByEnergy( (double) $3 ); }
        |
	TOKJOINENERGY TOKEQUALS FLOAT TOKNEWLINE
{ GlobalOptions->setJoinRateByEnergy( $3 ); }
        ;
