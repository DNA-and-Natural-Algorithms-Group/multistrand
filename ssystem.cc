/*
  Copyright (c) 2007-2010 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#include "ssystem.h"
#include "python_options.h"
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>

/*Options *GlobalOptions = NULL;
  extern int yyparse( void );
  extern int yyrestart ( FILE *);
*/
#ifndef PYTHON_THREADS
SimulationSystem::SimulationSystem( int argc, char **argv )
{
  assert(0);  // entire function is fail now.
  
  // FILE *fp = NULL;
  // if( GlobalOptions == NULL )
  //   {
  //     //      GlobalOptions = new Options();

  //     if( argc > 1 ) // some command line arguments.
  //       {
  //         if( argc > 2 ) // could be a filename arg
  //           {
  //             if( strcmp( argv[1], "--inputfile") == 0 )
  //               {
  //                 fp = fopen( argv[2], "rt");
  //                 if( fp == NULL )
  //                   fprintf(stderr,"ERROR: Could not open input file (%s) for reading.\n",argv[2]);
  //                 else
  //                   GlobalOptions = new Options( argv[2] );
  //               }
  //             else if( strcmp( argv[1], "--energy") == 0 )
  //               {
  //                 fp = fopen( argv[2], "rt");
  //                 if( fp == NULL )
  //                   {
  //                     fprintf(stderr,"ERROR: Could not open input file (%s) for reading.\n",argv[2]);
  //                     exit(1);
  //                   }
  //                 else
  //                   GlobalOptions = new Options( argv[2] );
          
  //                 GlobalOptions->setEnergyMode(1); // activate the energy-only version.
  //               }
  //           }
  //       }
  //     else
  //       GlobalOptions = new Options();

  //     if( fp != NULL )
  //       yyrestart( fp );
  //     yyparse();
  //     if( fp != NULL )
  //       fclose(fp);
  //   }
  
  // system_options = GlobalOptions;

  // if( Loop::GetEnergyModel() == NULL)
  //   {
  //     dnaEnergyModel = NULL;
  //     if( getLongAttr(system_options, parameter_type) == 0 ) // VIENNA = 0
  //       dnaEnergyModel = new ViennaEnergyModel( system_options );
  //     else
  //       dnaEnergyModel = new NupackEnergyModel( system_options );
  //     Loop::SetEnergyModel( dnaEnergyModel );
  //   }
  // else
  //   {
  //     dnaEnergyModel = Loop::GetEnergyModel();
  //   }
  
  // startState = NULL;
  // complexList = NULL; // new SComplexList( dnaEnergyModel );

  // //  GlobalOptions->finalizeInput();

  // getBoolAttr( system_options, boltzmann_sample, &boltzmann_sampling );
}
#endif
// end #ifndef PYTHON_THREADS - we do not want to have this constructor as it uses yyparse.

SimulationSystem::SimulationSystem( PyObject *system_options )
{
    
  if( Loop::GetEnergyModel() == NULL)
    {
      dnaEnergyModel = NULL;

      if(   getLongAttr(system_options, parameter_type) == 0 )
        dnaEnergyModel = new ViennaEnergyModel( system_options );
      else
        dnaEnergyModel = new NupackEnergyModel( system_options );
      Loop::SetEnergyModel( dnaEnergyModel );
    }
  else
    {
      dnaEnergyModel = Loop::GetEnergyModel();
    }
  
  startState = NULL;
  complexList = NULL; // new SComplexList( dnaEnergyModel );
  
  //  system_options->finalizeInput();
  getBoolAttr( system_options, boltzmann_sample, &boltzmann_sampling );
}

SimulationSystem::~SimulationSystem( void )
{
  if( complexList != NULL )
    delete complexList;
#ifndef PYTHON_THREADS
  if( dnaEnergyModel != NULL )
    delete dnaEnergyModel;
#endif

}


//#define SRANDOMDEV
// uncomment above line if you have srandomdev() available for better
// random number generation. 

void SimulationSystem::StartSimulation( void )
{
  int curcount = 0;
  long random_seed = 0;
  int ointerval = getLongAttr(system_options, output_interval);
  int initial_seed = getLongAttr(system_options, initial_seed);
  FILE *fp; // used only for initial random number generation.
  //#ifndef SRANDOMDEV
  //  srandom( time( NULL) );

  if(getLongAttr(system_options, simulation_mode) & SIMULATION_PYTHON)
    {
      callFunc_NoArgsToNone(system_options, reset_completed_python);
      setDoubleAttr(system_options, python_collision_rate, -1.0);
    }

  if((fp = fopen("/dev/urandom","r")) != NULL )
    {  // if urandom exists, use it to provide a seed
      long deviceseed;
      fread(&deviceseed, sizeof(long), 1, fp);
      
      srand48( deviceseed );
      fclose(fp);
      //printf("device seeded: %ld\n",deviceseed);
    }
  else // use the possibly flawed time as a seed.
    {
      srand48( time(NULL) );
      //printf("time seeded\n");
    }
  //#endif

  if( getLongAttr(system_options, energy_mode) == 1 ) // need energy only.
    {
      InitializeSystem();
      complexList->initializeList();
      complexList->printComplexList(1); // NUPACK energy output : bimolecular penalty, no Volume term.
      return;
    }
  
  if( getLongAttr(system_options, simulation_mode) & SIMULATION_FIRST )
    {
      StartSimulation_First_Bimolecular();
      return;
    }

  
  while( curcount < getLongAttr(system_options, num_simulations) )
    {
      
      //#ifdef SRANDOMDEV
      //srandomdev(); // initialize random() generator from /dev/random device.
      //#endif
      //      random_seed = random(); // get a seed from that generator
      random_seed = lrand48();

      if( curcount == 0 && initial_seed > 0)
        random_seed = initial_seed;
       
      //srandom( random_seed ); // re-init generator using that seed
      // Sort of a hack, but I need to have
      // the sequences repeatable.
      srand48( random_seed );


      if( ointerval < 0 && !(getLongAttr(system_options, simulation_mode) & SIMULATION_PYTHON))
        printf("Seed: 0x%lx\n",random_seed);
      InitializeSystem();

      if( ointerval < 0 && !(getLongAttr(system_options, simulation_mode) & SIMULATION_PYTHON))
        {
          printf("Size: %ld, %ld, %ld, %ld\n",sizeof(time_t),sizeof(long),RAND_MAX,1<<31 -1);
          printf("System %d Initialized\n",curcount);
        }
      if( getLongAttr(system_options, trajectory_type) > 0 )
        system_options->printTrajLine( NULL, curcount );

    
      SimulationLoop( random_seed );
      curcount++;
    }
}

#ifdef PYTHON_THREADS
void SimulationSystem::StartSimulation_threads( void )
{
  using namespace boost::python;
  Py_BEGIN_ALLOW_THREADS
    StartSimulation();
  Py_END_ALLOW_THREADS
    }
#endif


void SimulationSystem::SimulationLoop( long r_seed )
{
  double rchoice,rate,stime=0.0;
  double maxsimtime = getDoubleAttr(system_options, simulation_time);
  int stopcount = getLongAttr(system_options, stop_count);
  int stopoptions = getLongAttr(system_options, stop_options);
  int curcount = 0;
  int checkresult = 0;
  class stopcomplexes *traverse;
  int ointerval = getLongAttr(system_options, output_interval);
  int sMode = getLongAttr(system_options, simulation_mode); 
  double otime = getDoubleAttr(system_options, output_time);
  double ctime = 0.0;
  //  long temp_r=0;
  // double temp_deltat=0.0;
  complexList->initializeList();

  rate = complexList->getTotalFlux();
 
  if( getLongAttr(system_options, trajectory_type) == 0)
    {
      do {
        //assert( rate > -0.0 );
        //  rchoice = (rate * random()/((double)RAND_MAX+1));
        rchoice = rate * drand48();
        //  temp_r = random();
        //temp_deltat = (log(1. / (double)((temp_r + 1)/(double)RAND_MAX)) / rate);
        //if( isnan( temp_deltat ) )
        /// printf("Stime = NAN\n");
        //  stime += (log( 1. / ( ((double)random() + 1.0) /((double)RAND_MAX+1.0))) / rate );
        stime += (log( 1. / (1.0 - drand48()) ) / rate ); // 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.

        //temp_deltat; //
        //  assert( stime > -0.0 );

        if( sMode & SIMULATION_PYTHON )
          setDoubleAttr(system_options, python_current_time, stime);

        complexList->doBasicChoice( rchoice, stime );
        rate = complexList->getTotalFlux();

        if( stopcount > 0 && stopoptions==2)
          {
            curcount = 0;
            checkresult = 0;
            while( curcount < stopcount && checkresult == 0 )
              {
                traverse = getStopComplexList( system_options, curcount );
                checkresult = complexList->checkStopComplexList( traverse->citem );
                curcount++;
              }
          }
        if( checkresult == 0 && (sMode & SIMULATION_PYTHON) )
          {
            // if external python interface has asked us to complete.
            if( getLongAttr(system_options, python_trajectory_suspend_flag))
              {
                while( getLongAttr(system_options, python_trajectory_suspend_flag) )
                  sleep(1);
              }

            if( getLongAttr(system_options, python_trajectory_halt_flag) )
              checkresult = -1;

          }
    
        // trajectory output via outputtime option
        if( otime > 0.0 )
          {
            if( stime - ctime > otime )
              {
                ctime += otime;
                printf("Current State: Choice: %6.4f, Time: %6.6f\n",rchoice, ctime);
                complexList->printComplexList( 0 );
              }

          }

        // trajectory output via outputinterval option
        if( ointerval >= 0 )
          {
            if( getLongAttr(system_options, output_state) )
              {
                printf("Current State: Choice: %6.4f, Time: %6.6f\n",rchoice, stime);
                complexList->printComplexList( 0 );
              }
            callFunc_NoArgsToNone(system_options, increment_output_state);
          }

      } while( /*rate > 0.01 &&*/ stime < maxsimtime && checkresult == 0);
      
      if( otime > 0.0 )
        printf("Final state reached: Time: %6.6f\n",stime);
      if( ! (sMode & SIMULATION_PYTHON ) )
        if( ointerval < 0 || getLongAttr(system_options, output_state) == 0)
          complexList->printComplexList( 0 );
      
      if( stime == NAN )
        system_options->printStatusLine( r_seed, "ERROR", stime );
      else if ( checkresult > 0 )
        system_options->printStatusLine( r_seed, traverse->tag, stime );
      else
        system_options->printStatusLine( r_seed, "INCOMPLETE", stime );
      
      if( ! (sMode & SIMULATION_PYTHON) )
        printf("Trajectory Completed\n");
    }
  else if( getLongAttr(system_options, trajectory_type) == 1 )
    {
      int stopindex=-1;
      if( stopcount > 0 )
        {
          curcount = 0;
          while( curcount < stopcount && stopindex < 0 )
            {
              traverse = getStopComplexList( system_options, curcount );
              if( strstr( traverse->tag, "stop") != NULL )
                stopindex = curcount;
              curcount++;
            }
        }

      system_options->printTrajLine("Start",0);
      do {
        rchoice = (rate * random()/((double)RAND_MAX));
        stime += (log(1. / (double)((random() + 1)/(double)RAND_MAX)) / rate );
        complexList->doBasicChoice( rchoice, stime );
        rate = complexList->getTotalFlux();
        if( ointerval >= 0 )
          {
            if( getLongAttr(system_options, output_state) )
              complexList->printComplexList( 0 );
            callFunc_NoArgsToNone(system_options, increment_output_state);
          }
        if( stopcount > 0 )
          {
            curcount = 0;
            checkresult = 0;
            while( curcount < stopcount && checkresult == 0 )
              {
                traverse = getStopComplexList( system_options, curcount );
                checkresult = complexList->checkStopComplexList( traverse->citem );
                curcount++;
              }
          }
        if( checkresult > 0 )
          {
            system_options->printTrajLine( traverse->tag, stime );
          }
        else
          system_options->printTrajLine("NOSTATE", stime );
      } while( /*rate > 0.01 && */ stime < maxsimtime && !(checkresult > 0 && stopindex == curcount-1));
      
      if( ointerval < 0 || getLongAttr(system_options, output_state) == 0)
        complexList->printComplexList( 0 );

      if( stime == NAN )
        system_options->printStatusLine( r_seed, "ERROR", stime );
      else
        system_options->printStatusLine( r_seed, "INCOMPLETE", stime );
      
      printf("Trajectory Completed\n");
    }
}


void SimulationSystem::StartSimulation_First_Bimolecular( void )
{
  int curcount = 0;
  long random_seed = 0;
  int ointerval = getLongAttr(system_options, output_interval);
  int initial_seed = getLongAttr(system_options, initial_seed);

  /* these are used for compiling statistics on the runs */
  double completiontime;
  int completiontype;
  double forwardrate;

  char *tag; // tracking stop state tag for output reasons

  /* accumulators for our total types, rates and times. */
  double total_rate = 0.0 ;
  double total_time[2] = {0.0, 0.0};
  long total_types[3] = {0, 0,0};
  double computed_rate_means[3] = {0.0, 0.0,0.0};
  double computed_rate_mean_diff_squared[3] = {0.0, 0.0,0.0}; // difference from current mean, squared, used for online variance algorithm.
  double delta = 0.0; // temporary variable for computing  delta from current mean.

  int sMode = getLongAttr(system_options, simulation_mode) & SIMULATION_PYTHON;
  assert( getLongAttr(system_options, simulation_mode) & SIMULATION_FIRST );
  // random number seed has already been generated.

  while( curcount < getLongAttr(system_options, num_simulations) )
    {
      random_seed = lrand48();

      if( curcount == 0 && initial_seed > 0)
        random_seed = initial_seed;

      srand48( random_seed );

      InitializeSystem();

      //      if( getLongAttr(system_options, trajectory_type) > 0 )
      //    system_options->printTrajLine( NULL, curcount );
      //      if( sMode() )
      //        callFunc_NoArgsToNone(system_options, reset_completed);


      SimulationLoop_First_Bimolecular( random_seed, &completiontime, &completiontype, &forwardrate, &tag );

      // now we need to process the rate, time and type information.
      total_rate = total_rate + forwardrate;
      delta = forwardrate - computed_rate_means[2];
      computed_rate_means[2] += delta / (double) (curcount+1);
      computed_rate_mean_diff_squared[2] += delta * ( forwardrate - computed_rate_means[2]);
      
      if( completiontype == STOPCONDITION_TIME || completiontype == STOPCONDITION_FORWARD)
        {
          total_time[0] = total_time[0] + completiontime;
          total_types[0]++;
          if( completiontype == STOPCONDITION_TIME )
            total_types[2]++;
      
          delta = 1.0 / completiontime - computed_rate_means[0];
          computed_rate_means[0] += delta / (double ) total_types[0];
          computed_rate_mean_diff_squared[0] += delta * ( 1.0 / completiontime - computed_rate_means[0]);

          completiontype = STOPCONDITION_FORWARD;
        }
      if( completiontype == STOPCONDITION_REVERSE )
        {
          total_time[1] = total_time[1] + completiontime;
          total_types[1]++;

          delta = 1.0 / completiontime - computed_rate_means[1];
          computed_rate_means[1] += delta / (double ) total_types[1];
          computed_rate_mean_diff_squared[1] += delta * (1.0 / completiontime - computed_rate_means[1]);

        }

      system_options->printStatusLine_First_Bimolecular( random_seed, completiontype, completiontime, forwardrate, tag );
      curcount++;
    }
  if( !sMode )
    system_options->printStatusLine_Final_First_Bimolecular( total_rate, total_time, total_types, curcount, computed_rate_means, computed_rate_mean_diff_squared );
  if( total_types[2] > 0 && !sMode )
    system_options->printStatusLine_Warning( 0, total_types[2] );

}


void SimulationSystem::SimulationLoop_First_Bimolecular( long r_seed, double *completiontime, int *completiontype, double *frate, char **tag )
{
  double rchoice,rate,stime=0.0;
  double maxsimtime = getDoubleAttr(system_options, simulation_time);
  int stopcount = getLongAttr(system_options, stop_count);
  int stopoptions = getLongAttr(system_options, stop_options);
  int curcount = 0;
  int checkresult = 0;
  class stopcomplexes *traverse;
  int ointerval = getLongAttr(system_options, output_interval);
  int sMode = getLongAttr(system_options, simulation_mode) & SIMULATION_PYTHON;
  int trajMode = getLongAttr(system_options, trajectory_type);
  double otime = getDoubleAttr(system_options, output_time);
  double ctime = 0.0;
  //  long temp_r=0;
  // double temp_deltat=0.0;
  complexList->initializeList();

  rate = complexList->getJoinFlux();

  if ( rate < 0.0 )
    { // no initial moves
      *completiontime = 0.0;
      *completiontype = STOPCONDITION_ERROR;
      return;
    }

  rchoice = rate * drand48();

  if( ointerval >= 0 || getDoubleAttr(system_options, output_time) > 0.0 )
    complexList->printComplexList( 0 );
  if( ointerval >= 0 )
    printf("Initial State (before join).\n",rchoice, ctime);

  complexList->doJoinChoice( rchoice );

  *frate = rate * dnaEnergyModel->getJoinRate_NoVolumeTerm() / dnaEnergyModel->getJoinRate() ; // store the forward rate used for the initial step.

  if( ointerval >= 0 || getDoubleAttr(system_options, output_time) > 0.0 )
    {
      complexList->printComplexList( 0 );
      printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, ctime);
    }

  if( sMode )
    setDoubleAttr(system_options, python_collision_rate, *frate);
 
  // Begin normal steps.
  rate = complexList->getTotalFlux();
  do {

    rchoice = rate * drand48();
    stime += (log( 1. / (1.0 - drand48()) ) / rate ); // 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.

    if( !sMode && otime > 0.0 )
      {
        if( stime - ctime > otime )
          {
            ctime += otime;
            printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, ctime);
            complexList->printComplexList( 0 );
          }

      }

    if( sMode )
      setDoubleAttr(system_options, python_current_time, stime);
      

    complexList->doBasicChoice( rchoice, stime );
    rate = complexList->getTotalFlux();
    
    if( !sMode && ointerval >= 0 )
      {
        if( getLongAttr(system_options, output_state) )
          {
            complexList->printComplexList( 0 );
            printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, stime);
          }
        callFunc_NoArgsToNone(system_options, increment_output_state);
      }

    if( stopcount > 0 && stopoptions==2)
      {
        curcount = 0;
        checkresult = 0;
        while( curcount < stopcount && checkresult == 0 )
          {
            traverse = getStopComplexList( system_options, curcount );
            checkresult = complexList->checkStopComplexList( traverse->citem );
            curcount++;
          }
      }
    if( checkresult == 0 && sMode )
      {
        // if external python interface has asked us to complete.
        if( getLongAttr(system_options, python_trajectory_suspend_flag))
          {
            while( getLongAttr(system_options, python_trajectory_suspend_flag) )
              sleep(1);
          }

        if( getLongAttr(system_options, python_trajectory_halt_flag) )
          checkresult = -1;

      }
  } while( stime < maxsimtime && checkresult == 0);
      
  if( !sMode && otime > 0.0 )
    {
      complexList->printComplexList(0);
      printf("Final state reached: Time: %6.6e\n",stime);
    }
  // if( ointerval >= 0 )
  //printf("Final state reached: Time: %6.6e\n",stime);

    
  /*   if( stime == NAN )
       system_options->printStatusLine( r_seed, "ERROR", stime );
       else if ( checkresult > 0 )
       system_options->printStatusLine( r_seed, traverse->tag, stime );
       else
       system_options->printStatusLine( r_seed, "INCOMPLETE", stime );
  */
  // printing is handled at the upper level for the status information. We do, however, need to return the info.

  if( checkresult > 0 )
    {
      if( strcmp( traverse->tag, "REVERSE") == 0 )
        *completiontype = STOPCONDITION_REVERSE;
      else 
        {
          *tag = traverse->tag;
          *completiontype = STOPCONDITION_FORWARD;
        }
      *completiontime = stime;
    }
  else
    {
      *completiontime = stime;
      *completiontype = STOPCONDITION_TIME;
    }
  if( ! sMode )
    printf("Trajectory Completed\n");
}


///////////////////////////////////////////////////
// This has now been ref count checked, etc etc. //
///////////////////////////////////////////////////

        
void SimulationSystem::InitializeSystem( void )
{                             
  class StrandComplex *tempcomplex;
  char *sequence, *structure;
  class identlist *id;
  int start_count;
  PyObject *py_start_state, *py_complex;
  PyObject *py_seq, *py_struc;

  if( complexList != NULL )
    delete complexList;

  complexList = new SComplexList( dnaEnergyModel );
  
  py_start_state = getListAttr(system_options, start_state);
  start_count = PyList_GET_SIZE(py_start_state);  
  // doesn't need reference counting
  
  for( int index = 0; index < start_count; index++ )
    {
      py_complex = PyList_GET_ITEM(py_start_state, index);
      // does need to decref when no longer in use.

      sequence = getStringAttr(py_complex, sequence, py_seq);
      
      //  TODO: simmode is now stored internally to ssystem.
      if ( getLongAttr(system_options, simulation_mode) == SIMULATION_FIRST_BIMOLECULAR )
        structure = getStringAttr(py_complex, boltzmann_structure, py_struc);
      else
        structure = getStringAttr(py_complex, structure, py_seq);
      
      id = getID_list( system_options, index );
      
      tempcomplex = new StrandComplex( sequence, structure, id );
      // StrandComplex does make its own copy of the seq/structure, so we can now decref.
      
      Py_DECREF( py_seq );
      Py_DECREF( py_struc );
      Py_DECREF( py_complex );
      startState = tempcomplex;
      complexList->addComplex( tempcomplex );
      tempcomplex = NULL;
    }
  return;
}


int SimulationSystem::ResetSystem( void )
{
  char insequence[3072],instructure[3072];
  char *newseq, *newstruc;
  int totlength, curstart, curindex, curcount;

  strcpy( insequence , systemsequence  );
  strcpy( instructure , systemstructure);

  if( complexList != NULL ) delete complexList;
  complexList = new SComplexList( dnaEnergyModel );
  startState = NULL;
  
  totlength = strlen(insequence);

  curindex = 0; curstart = 0; curcount = 0;
  while( curindex < totlength )
    {
      while( (curstart < totlength) && ((curcount > 0) || (insequence[curstart] != '_')))
        {
          if( instructure[curstart] == '(' ) curcount++;
          if( instructure[curstart] == ')' ) curcount--;
          if( curcount < 0 ) printf("Mismatched Parenthesis in input.");
          curstart++;
        }
      newseq = new char[curstart-curindex+1];
      newstruc = new char[curstart-curindex+1];
      insequence[curstart] = '\0';
      instructure[curstart] = '\0';
      strcpy( newseq, &insequence[curindex] );
      strcpy( newstruc, &instructure[curindex] );
      StrandComplex *tempcomplex = new StrandComplex( newseq, newstruc);
      delete[] newseq;
      delete[] newstruc;
      while( insequence[curstart] == '_' || (curstart < totlength && insequence[curstart] == '\0')) curstart++;
      curindex = curstart;
      curcount = 0;
      startState = tempcomplex;
      complexList->addComplex( tempcomplex );
    }

}


int SimulationSystem::StartSimulation( int input_flags, int num_sims, double simtime )
{
  // DEPRECATED: simulations always start via StartSimulation( void ) now,
  // and args are passed via the constructor. 
  Move *tempmove;
  StrandComplex *newcomplex = startState;
  srand( time( NULL ) );

  if( input_flags > 1 && input_flags < 4)
    {
      newcomplex->generateLoops();
      double energy = newcomplex->getEnergy();
      newcomplex->moveDisplay();      
      double rate = newcomplex->getTotalFlux();
      printf("%s\n%s (%6.2f) %6.4f\n",newcomplex->getSequence(),newcomplex->getStructure(),energy,rate);
      double stime = 0.0;
      
      if( input_flags == 2 )
        {
          double rchoice = (rate*rand()/((double)RAND_MAX));
      
          printf("%3.3f\n",rchoice);
          tempmove = newcomplex->getChoice( &rchoice );
          newcomplex->doChoice(tempmove);
          char *struc = newcomplex->getStructure();
          rate = newcomplex->getTotalFlux();
          energy = newcomplex->getEnergy();
          printf("%s\n%s (%6.2f) %6.4f\n",newcomplex->getSequence(),struc,energy,rate);
        }
      if( input_flags == 3 )
        {
          double rchoice,rsave,moverate;
          char *struc;
          int temp1, type;
          do {
            rsave = rchoice = (rate * rand()/((double)RAND_MAX));
            stime += (log(1. / (double)((rand() + 1)/(double)RAND_MAX)) / rate );
            tempmove = newcomplex->getChoice( &rchoice );
            moverate = tempmove->getRate();
            type = tempmove->getType();
            newcomplex->doChoice( tempmove );
            struc = newcomplex->getStructure();
            rate = newcomplex->getTotalFlux();
            energy = newcomplex->getEnergy();
            printf("%s (%6.2f) %6.4f %6.4f pc: %6.4f time: %6.4f type: %d\n",struc,energy,rate,moverate,rsave,stime,type);
            temp1=0;
            for( int loop = 0; loop < strlen(struc); loop++)
              {
                if(struc[loop] == '(') temp1++;
                if(struc[loop] == ')') temp1--;
                assert( temp1 >= 0 );
                if(loop > 0)
                  assert( !(struc[loop-1]=='(' && struc[loop] == ')'));
              }
        
            assert( temp1 == 0 );
            // printf("%p",tempmove);
          } while( rate > 0.10 && stime < 1000);
        }
    }
  if( input_flags == 4 ) // multi-complex operation mode - print start
    {
      complexList->initializeList();
      complexList->printComplexList( 0 );
      printf("That's it!\n");
    }
  if( input_flags == 5 ) // multi-complex operation mode - run basic moves - verbose output
    {
      int curcount = 0;
      
      while( curcount < num_sims )
        {
          complexList->initializeList();
          complexList->printComplexList( 0 );

          double rchoice,rate,stime=0.0;
      
          rate = complexList->getTotalFlux();
          do {
            rchoice = (rate * rand()/((double)RAND_MAX));
            stime += (log(1. / (double)((rand() + 1)/(double)RAND_MAX)) / rate );
            complexList->doBasicChoice( rchoice, stime );
            rate = complexList->getTotalFlux();
            complexList->printComplexList( 0 );
          } while( rate > 0.10 && stime < simtime);
          complexList->printComplexList( 0 );
          curcount++;
          if( curcount != num_sims)
            ResetSystem();
        }
    }
  
  if( input_flags == 6 ) // multi-complex operation mode - run basic moves - silent output
    {
      int curcount = 0;
      
      while( curcount < num_sims )
        {
          complexList->initializeList();
          //      complexList->printComplexList();

          double rchoice,rate,stime=0.0;
      
          rate = complexList->getTotalFlux();
          do {
            rchoice = (rate * rand()/((double)RAND_MAX));
            stime += (log(1. / (double)((rand() + 1)/(double)RAND_MAX)) / rate );
            complexList->doBasicChoice( rchoice, stime );
            rate = complexList->getTotalFlux();
            //      complexList->printComplexList();
          } while( rate > 0.10 && stime < simtime);
          complexList->printComplexList( 0 );
          printf("Trajectory Complete: Time: %.2lf\n",stime);
          curcount++;
          if( curcount != num_sims)
            ResetSystem();
        }
    }
  
}
