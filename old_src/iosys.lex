%{
/*
   Copyright (c) 2007-2008 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/
 
#include <stdio.h>
  //#include "logic.h"
#include "iosys.tab.h"
#include <string.h>
  extern int curLine;

%}

%x comment statement strands startstructure stopstructure numentry stringentry hexentry energymodel

%option noyywrap
%%

^"##"                             BEGIN(comment);
^"//"                             BEGIN(comment);
^"%"                              BEGIN(comment);

<comment>[^\n]*                  /* destroy all characters up til a \n */
<comment>\n                      BEGIN(INITIAL); curLine++;

<INITIAL>"#"                     BEGIN(statement); //return TOKPERCENT;

<statement>[ \t]*                /* eat up whitespace */
<statement>("Strands"|"strands")[ \t]*\n      BEGIN(strands); curLine++; return TOKSTRANDS;
<statement>("StartStructure"|"startstructure")[ \t]*\n  BEGIN(startstructure); curLine++;  return TOKSTARTSTRUCTURE;
<statement>("stopstructure"|"StopStructure"|"StopStructures"|"stopstructures")[ \t]*\n    BEGIN(stopstructure);  curLine++; return TOKSTOPSTRUCTURE;
<statement>("SimTime"|"simtime")[ \t]*                BEGIN(numentry);  return TOKSIMTIME;
<statement>("numsims"|"NumSims")[ \t]*                BEGIN(numentry);  return TOKNUMSIMS;
<statement>("OutputInterval"|"outputinterval")[ \t]*  BEGIN(numentry);  return TOKOUTPUTINTERVAL;
<statement>("OutputTime"|"outputtime")[ \t]*  BEGIN(numentry);  return TOKOUTPUTTIME;
<statement>("StopOptions"|"stopoptions")[ \t]*        BEGIN(numentry);  return TOKSTOPOPTIONS;
<statement>("SimulationMode"|"simulationmode")[ \t]*  BEGIN(numentry); return TOKSIMMODE;
<statement>("Logfile"|"logfile")[ \t]*                BEGIN(stringentry);  return TOKLOGFILE;
<statement>"dangles"[ \t]*         BEGIN(numentry); return TOKDANGLES;
<statement>("Trajfile"|"trajfile")[ \t]*              BEGIN(stringentry);  return TOKTRAJFILE;
<statement>("Trajtype"|"Trajtype")[ \t]*              BEGIN(numentry);  return TOKTRAJTYPE;
<statement>("Energymodelfile"|"energymodelfile")[ \t]*      BEGIN(stringentry);   return TOKPARAMFILE;
<statement>("Energymodeltype"|"energymodeltype")[ \t]*      BEGIN(numentry);  return TOKPARAMTYPE;
<statement>("Energymodel"|"energymodel")[ \t]*                      BEGIN(energymodel);  return TOKENERGYMODEL;
<statement>("NupackEnableGTWobble"|"nupackenablegtwobble")[ \t]*        BEGIN(numentry);  return TOKGTWOBBLE; 
<statement>("LogML"|"logml")[ \t]*        BEGIN(numentry); return TOKLOGML;
<statement>("StartSeed"|"Startseed"|"startseed")[ \t]*            BEGIN(hexentry); return TOKSTARTSEED;
<statement>("JoinEnergy"|"joinenergy"|"Joinenergy")[ \t]* BEGIN(numentry); return TOKJOINENERGY;
<statement>("Concentration"|"concentration")[ \t]*    BEGIN(numentry); return TOKCONCENTRATION;
<statement>("Ratemethod"|"ratemethod")[ \t]*    BEGIN(numentry); return TOKRATEMETHOD;
<statement>("Intra"|"intra"|"IntramolecularScale"|"intramolecularscale")[ \t]*    BEGIN(numentry); return TOKINTRA;
<statement>("Inter"|"inter"|"IntermolecularScale"|"intermolecularscale")[ \t]*    BEGIN(numentry); return TOKINTER;
<statement>("Temperature"|"temperature")[ \t]*        BEGIN(numentry); return TOKTEMPERATURE;
<statement>("Volume"|"volume")[ \t]*                  BEGIN(numentry); return TOKVOLUME;
<statement>("JoinRate"|"joinrate")[ \t]*              BEGIN(numentry); return TOKJOINRATE;
<statement>{
^"##"                             BEGIN(comment);
^"//"                             BEGIN(comment);
^"%"                              BEGIN(comment);
    [^\n]
      \n             BEGIN(INITIAL); curLine++; return TOKNEWLINE;
}

<strands>{
[AGCTUagctu]+                     yylval.buffer = strdup(yytext); yylval.buffer[yyleng]='\0'; return SEQUENCE;
[[:alpha:]_\-][[:alnum:]_\-]*     yylval.buffer = strdup(yytext); yylval.buffer[yyleng]='\0'; return ID;
 ","                              return TOKCOMMA;
[0-9]*"."[0-9]+                   yylval.fnum = atof(yytext); return FLOAT;
[0-9]+                            yylval.num = atoi(yytext); return NUMBER;
[ \t]*
  \n                                curLine++;return TOKNEWLINE;
^"##"                             BEGIN(comment);
^"//"                             BEGIN(comment);
^"%"                              BEGIN(comment);
^"#"                              BEGIN(statement); //return TOKPERCENT;
}

<stopstructure>{
  "TAG:"             return TOKTAGS;
  "BOUND"            return TOKBOUND;
  "DISASSOC"         return TOKDISASSOC;
}

<startstructure,stopstructure>{
  [[:alpha:]_\-][[:alnum:]_\-]*   yylval.buffer = strdup(yytext); yylval.buffer[yyleng] = '\0'; return ID;
  [\(\)\._\+\*]+                       yylval.buffer = strdup(yytext); yylval.buffer[yyleng] = '\0'; return STRUCTURE;
  [0-9]*"."[0-9]+                   yylval.fnum = atof(yytext); return FLOAT;
  "-"?[0-9]+                        yylval.num = atoi(yytext); return NUMBER;
  ","                             return TOKCOMMA;

  ^"%"                              BEGIN(comment);

  "%"                             return TOKPERCENT;
  \n                              curLine++;return TOKNEWLINE;
  [ \t]*
  ^"##"                             BEGIN(comment);
  ^"//"                             BEGIN(comment);
  ^"#"                            BEGIN(statement);
}

<numentry>{
  [0-9]*"."[0-9]+                   yylval.fnum = atof(yytext); return FLOAT;
  "-"?[0-9]+                        yylval.num = atoi(yytext); return NUMBER;
}

<hexentry>{
  "0x"?[0-9a-fA-F]+                 yylval.num = strtol( yytext, NULL, 0 ); return NUMBER;
}

<numentry,stringentry,hexentry,energymodel>{
  "="                             return TOKEQUALS;
  [ \t]*
  \n                              curLine++;return TOKNEWLINE;
  ^"##"                           BEGIN(comment);
  ^"//"                             BEGIN(comment);
  ^"%"                              BEGIN(comment);
  ^"#"                            BEGIN(statement);
}

<stringentry>{
  [[[:alnum:]_\+\.\-/]*     yylval.buffer = strdup(yytext); yylval.buffer[yyleng]='\0'; return ID;
}

<energymodel>{
    "VIENNA_DNA"              return TOKVIENNADNA;
    "NUPACK_DNA_2_3"          return TOKNUPACKDNA23;
    "NUPACK_RNA_2_3"          return TOKNUPACKRNA23;
}

[0-9]+           yylval.num=atoi(yytext); return NUMBER;
","              return TOKCOMMA;
\n               curLine++; return TOKNEWLINE; /* ignore EOL */
[ \t]+           /* ignore whitespace */
.                return STRING;/*[^\n]+           yylval.buffer = strdup(yytext); return ID;/*return (int) yytext[0];*/
%%
