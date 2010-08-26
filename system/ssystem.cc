/*
  Copyright (c) 2007-2010 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#include "python_options.h"
#include "ssystem.h"

#include <string.h>
#include <time.h>
#include <stdlib.h>

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
  //     if( testLongAttr(system_options, parameter_type ,=, 0 ) ) // VIENNA = 0
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

SimulationSystem::SimulationSystem( PyObject *system_o )
{
  system_options = system_o;
  
  getLongAttr(system_options, simulation_mode, &simulation_mode );
  getLongAttr(system_options, num_simulations, &simulation_count_remaining);
  if( Loop::GetEnergyModel() == NULL)
    {
      dnaEnergyModel = NULL;

      if(  testLongAttr(system_options, parameter_type,=,0) )
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
  //  getBoolAttr( system_options, boltzmann_sample, &boltzmann_sampling );
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

void SimulationSystem::StartSimulation( void )
{
  FILE *fp; // used only for initial random number generation.
  int curcount = 0;
  long ointerval;

  getLongAttr(system_options, output_interval,&ointerval);


  if( simulation_mode & SIMULATION_MODE_FLAG_PYTHON )
    {
      pingAttr( system_options, interface_reset_completion_flag );
    }

  if( simulation_mode & SIMULATION_MODE_ENERGY_ONLY) // need energy only.
    {
      InitializeSystem();
      complexList->initializeList();
      complexList->printComplexList(1); // NUPACK energy output : bimolecular penalty, no Volume term.
      return;
    }
  
  if( simulation_mode & SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR )
    {
      StartSimulation_First_Bimolecular();
      return;
    }

  InitializeRNG();

  while( simulation_count_remaining > 0)
    {
      InitializeSystem();

      if( ointerval < 0 && !(simulation_mode & SIMULATION_MODE_FLAG_PYTHON))
        {
          printf("System %d Initialized\n",curcount);
        }
      //if( getLongAttr(system_options, trajectory_type) > 0 )
      //_m_printTrajLine(system_options, NULL, curcount );
      // Currently deactivated til prints are ready.
    
      SimulationLoop();
      simulation_count_remaining--;
      pingAttr( system_options, increment_trajectory_count );
      generateNextRandom();
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


void SimulationSystem::SimulationLoop( void )
{
  double rchoice,rate,stime=0.0;
  class stopcomplexes *traverse;
  int curcount = 0;
  int checkresult = 0;
  double ctime = 0.0;
  double maxsimtime;
  long stopcount;
  long stopoptions;
  long ointerval;
  long sMode;
  double otime;

  getLongAttr(system_options, simulation_mode,&sMode); 
  getLongAttr(system_options, output_interval,&ointerval);
  getLongAttr(system_options, use_stop_conditions,&stopoptions);
  getLongAttr(system_options, stop_count,&stopcount);
  getDoubleAttr(system_options, simulation_time,&maxsimtime);
  getDoubleAttr(system_options, output_time,&otime);

  complexList->initializeList();

  rate = complexList->getTotalFlux();
 
  if( testLongAttr(system_options, trajectory_type ,=, 0 ))
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

        if( sMode & SIMULATION_MODE_FLAG_PYTHON )
          setDoubleAttr( system_options, interface_current_time, stime);
        // replaces callFunc_DoubleToNone(...).

        complexList->doBasicChoice( rchoice, stime );
        rate = complexList->getTotalFlux();

        if( stopcount > 0 && stopoptions)
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
        if( checkresult == 0 && (sMode & SIMULATION_MODE_FLAG_PYTHON) )
          {
            // previously: callFunc_NoArgsToLong(system_options, get_python_trajectory_suspend_flag)

            // if external python interface has asked us to complete.
            if( testBoolAttr(system_options, interface_suspend_flag) ) 
              {
                while( testBoolAttr(system_options, interface_suspend_flag) ) 

                  sleep(1);
              }

            if( testBoolAttr( system_options, interface_halt_flag ))
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
            if( testBoolAttr(system_options, output_state) )
              {
                printf("Current State: Choice: %6.4f, Time: %6.6f\n",rchoice, stime);
                complexList->printComplexList( 0 );
              }
            pingAttr( system_options, increment_output_state );
          }

      } while( /*rate > 0.01 &&*/ stime < maxsimtime && checkresult == 0);
      
      if( otime > 0.0 )
        printf("Final state reached: Time: %6.6f\n",stime);
      if( ! (sMode & SIMULATION_MODE_FLAG_PYTHON ) )
        if( ointerval < 0 || testLongAttr(system_options, output_state ,=, 0 ))
          complexList->printComplexList( 0 );
      
      if( stime == NAN )
        printStatusLine(system_options, current_seed, STOPRESULT_NAN, 0.0, NULL);
      else if ( checkresult > 0 )
        printStatusLine(system_options, current_seed, STOPRESULT_NORMAL, stime, traverse->tag );
      else
        printStatusLine(system_options, current_seed, STOPRESULT_TIME, stime, NULL );
      
      if( ! (sMode & SIMULATION_MODE_FLAG_PYTHON) )
        printf("Trajectory Completed\n");
    }
  else if( testLongAttr(system_options, trajectory_type ,=, 1 ) )
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

      //      printTrajLine(system_options,"Start",0);
      do {
        rchoice = (rate * random()/((double)RAND_MAX));
        stime += (log(1. / (double)((random() + 1)/(double)RAND_MAX)) / rate );
        complexList->doBasicChoice( rchoice, stime );
        rate = complexList->getTotalFlux();
        if( ointerval >= 0 )
          {
            if( testBoolAttr(system_options, output_state) )
              complexList->printComplexList( 0 );
            pingAttr( system_options, increment_output_state );
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
            ;// printTrajLine(system_options, traverse->tag, stime );
          }
        else
          ;
          //          printTrajLine(system_options,"NOSTATE", stime );
      } while( /*rate > 0.01 && */ stime < maxsimtime && !(checkresult > 0 && stopindex == curcount-1));
      
      if( ointerval < 0 || testLongAttr(system_options, output_state ,=, 0 ))
        complexList->printComplexList( 0 );

      if( stime == NAN )
        printStatusLine(system_options, current_seed, STOPRESULT_NAN, 0.0, NULL);
      else
        printStatusLine(system_options, current_seed, STOPRESULT_TIME, stime, NULL );
      
      //      printf("Trajectory Completed\n");
    }
}


void SimulationSystem::StartSimulation_First_Bimolecular( void )
{
  assert((simulation_mode & SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR)  &&  (simulation_mode & SIMULATION_MODE_FLAG_PYTHON));

  InitializeRNG();
  
  while( simulation_count_remaining > 0 )
    {
      setLongAttr( system_options, interface_current_seed, current_seed );
      InitializeSystem();

      SimulationLoop_First_Bimolecular();
      generateNextRandom();
      pingAttr( system_options, increment_trajectory_count );
      simulation_count_remaining--;
    }
}



void SimulationSystem::SimulationLoop_First_Bimolecular( void )
{
  double rchoice,rate,stime=0.0, ctime=0.0;
  int curcount = 0;
  int checkresult = 0;

  double maxsimtime;
  long stopcount;
  long stopoptions;
  class stopcomplexes *traverse;
  long ointerval;
  int sMode = simulation_mode & SIMULATION_MODE_FLAG_PYTHON;
  long trajMode;
  double otime;
  double otime_interval;
  double frate = 0.0;

  getLongAttr(system_options, trajectory_type,&trajMode);
  getLongAttr(system_options, output_interval,&ointerval);
  getDoubleAttr(system_options, output_time,&otime);
  getLongAttr(system_options, use_stop_conditions,&stopoptions);
  getLongAttr(system_options, stop_count,&stopcount);
  getDoubleAttr(system_options, simulation_time,&maxsimtime);
  getDoubleAttr(system_options, output_time, &otime_interval );

  complexList->initializeList();

  rate = complexList->getJoinFlux();

  // scomplexlist returns a 0.0 rate if there was a single complex in
  // the system, and a -1.0 rate if there are exactly 0 join moves. So
  // the 0.0 rate should probably be caught, though if you use a
  // single complex system for a starting state it's probably
  // deserved.

  if ( rate < 0.0 )
    { // no initial moves
      printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_NOMOVES, 0.0, 0.0, NULL );
      return;
    }

  rchoice = rate * drand48();

  if( ointerval >= 0 || otime_interval >= 0.0 )
    complexList->printComplexList( 0 );
  /*  if( ointerval >= 0 )
    printf("Initial State (before join).\n",rchoice, ctime);
  */

  complexList->doJoinChoice( rchoice );

  // store the forward rate used for the initial step so we can record it.
  frate = rate * dnaEnergyModel->getJoinRate_NoVolumeTerm() / dnaEnergyModel->getJoinRate() ; 

  // rate is the total flux across all join moves - this is exactly equal to total_move_count *
  // dnaEnergyModel->getJoinRate()

  // This join rate is the dG_volume * bimolecular scaling constant
  // used for forward transitions.  What we actually need is the
  // bimolecular scaling constant * total move count. (dG volume is
  // the volume dependent term that is not actually related to the
  // 'collision' rate, but rather the volume we are simulating.

  if( ointerval >= 0 || otime_interval > 0.0 )
    {
      complexList->printComplexList( 0 );
      printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, ctime);
    }

  // Begin normal steps.
  rate = complexList->getTotalFlux();
  do {
    
    rchoice = rate * drand48();
    stime += (log( 1. / (1.0 - drand48()) ) / rate ); 
    
    // 1.0 - drand as drand returns in the [0.0, 1.0) range, 
    //  we need a (0.0,1.0] range for this to be the correct
    // distribution - log 1/U => U in (0.0,1.0] => 1/U => (+Inf,0.0] => 
    //    natural log:   (log(+Inf), 1.0]

    if( !sMode && otime > 0.0 )
      {
        if( stime - ctime > otime )
          {
            ctime += otime;
            printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, ctime);
            complexList->printComplexList( 0 );
          }

      }

    //    if( sMode )
    //    setDoubleAttr(system_options, interface_current_time, stime);
      

    complexList->doBasicChoice( rchoice, stime );
    rate = complexList->getTotalFlux();
    
    if( !sMode && ointerval >= 0 )
      {
        if( testBoolAttr(system_options, output_state) )
          {
            complexList->printComplexList( 0 );
            printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, stime);
          }
        pingAttr( system_options, increment_output_state );
        
      }

    if( stopcount > 0 && stopoptions)
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
        // if external python interface has asked us to pause.
        if( testBoolAttr(system_options, interface_suspend_flag) ) 
          {
            while( testBoolAttr(system_options, interface_suspend_flag) ) 
              sleep(1);
          }
        
        if( testBoolAttr( system_options, interface_halt_flag ))
          checkresult = -1;
      }
  } while( stime < maxsimtime && checkresult == 0);
      
  if( !sMode && otime > 0.0 )
    {
      complexList->printComplexList(0);
      printf("Final state reached: Time: %6.6e\n",stime);
    }

  if( checkresult > 0 )
    {
      if( strcmp( traverse->tag, "REVERSE") == 0 )
        printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_REVERSE, stime, frate, traverse->tag );
      else 
        printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_FORWARD, stime, frate, traverse->tag );
    }
  else
    {
      printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_FTIME, stime, frate, NULL );
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

  startState = NULL;
  if( complexList != NULL )
    delete complexList;

  complexList = new SComplexList( dnaEnergyModel );
  
  py_start_state = getListAttr(system_options, start_state);
  start_count = PyList_GET_SIZE(py_start_state);  
  // doesn't need reference counting for this size call.
  // the getlistattr call we decref later.
  
  for( int index = 0; index < start_count; index++ )
    {
      py_complex = PyList_GET_ITEM(py_start_state, index);
      // does need to decref when no longer in use.

      sequence = getStringAttr(py_complex, sequence, py_seq);
      
      // boltzmann sampling is not our problem - the options object
      // should be providing structures according to the boltzmann or
      // not as it decides, not our job now. :)

      structure = getStringAttr(py_complex, structure, py_struc);
      
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
  Py_DECREF( py_start_state );
  return;
}

void SimulationSystem::InitializeRNG( void )
{
  bool use_fixed_random_seed = false;
  FILE *fp = NULL;
  getBoolAttr( system_options, initial_seed_flag, &use_fixed_random_seed);

  if( use_fixed_random_seed )
    getLongAttr(system_options, initial_seed,&current_seed);
  else
    {
      if((fp = fopen("/dev/urandom","r")) != NULL )
        {  // if urandom exists, use it to provide a seed
          long deviceseed;
          fread(&deviceseed, sizeof(long), 1, fp);

          current_seed = deviceseed;
          fclose(fp);
        }
      else // use the possibly flawed time as a seed.
        {
          current_seed = time(NULL);
        }
    }
  // now initialize this generator using our random seed, so that we can reproduce as necessary.
  srand48( current_seed );
}

void SimulationSystem::generateNextRandom( void )
{
  current_seed = lrand48();
  srand48( current_seed );
}
