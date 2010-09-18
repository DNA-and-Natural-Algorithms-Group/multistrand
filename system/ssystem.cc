/*
  Copyright (c) 2007-2010 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#include "options.h"
#include "ssystem.h"

#include <string.h>
#include <time.h>
#include <stdlib.h>

#ifdef PROFILING
#include "google/profiler.h"
#endif

SimulationSystem::SimulationSystem( int argc, char **argv )
{
  return;
  // entire function is FAIL. Needs replacing.
}

SimulationSystem::SimulationSystem( PyObject *system_o )
{
#ifdef PROFILING
  ProfilerStart("ssystem_init_profile.prof");
#endif

  system_options = system_o;
  Py_INCREF( system_options );

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
  complexList = NULL; 
#ifdef PROFILING
  ProfilerStop();
#endif
}

SimulationSystem::~SimulationSystem( void )
{
  if( complexList != NULL )
    delete complexList;  
  complexList = NULL;

  // the remaining members are not our responsibility, we null them out
  // just in case something thread-unsafe happens.
  dnaEnergyModel = NULL;

  if( system_options != NULL )
    Py_DECREF( system_options );

  system_options = NULL;
  startState = NULL;
}

void SimulationSystem::StartSimulation( void )
{
  FILE *fp; // used only for initial random number generation.
  int curcount = 0;
  long ointerval;
#ifdef PROFILING
  ProfilerStart("ssystem_run_profile.prof");
#endif

  getLongAttr(system_options, output_interval,&ointerval);

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

      // if( ointerval < 0 && !(simulation_mode & SIMULATION_MODE_FLAG_PYTHON))
      //   {
      //     printf("System %d Initialized\n",curcount);
      //   }
      //if( getLongAttr(system_options, trajectory_type) > 0 )
      //_m_printTrajLine(system_options, NULL, curcount );
      // Currently deactivated til prints are ready.
    
      SimulationLoop();
      simulation_count_remaining--;
      pingAttr( system_options, increment_trajectory_count );
      generateNextRandom();
    }
#ifdef PROFILING
  ProfilerStop();
#endif
}

void SimulationSystem::SimulationLoop( void )
{
  double rchoice,rate,stime,ctime;
  // Could really use some commenting on these local vars.
  rchoice = rate = stime = ctime = 0.0;

  double maxsimtime, otime;
  maxsimtime = otime = -1.0;
  
  int curcount = 0, checkresult = 0;
  long stopcount = 0, stopoptions = 0, sMode = 0;
  long ointerval = -1;
  class stopcomplexes *traverse = NULL, *first=NULL;

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
      // Note: trajectory type 0 is normal, trajectory type 1 is transition time data.
      // TODO: wrap this into simulation modes. 

      do {

        rchoice = rate * drand48();

        stime += (log( 1. / (1.0 - drand48()) ) / rate ); 
        // 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.
        // see notes below in First Step mode.

        complexList->doBasicChoice( rchoice, stime );
        rate = complexList->getTotalFlux();

        if( stopcount > 0 && stopoptions)
          {
            checkresult = 0;
            first = getStopComplexList( system_options, 0 );
            checkresult = complexList->checkStopComplexList( first->citem );
            traverse = first;
            while( traverse->next != NULL && checkresult == 0 )
              {
                traverse = traverse->next;
                checkresult = complexList->checkStopComplexList( traverse->citem );
              }
            // Note: we cannot delete first here if checkresult != 0,
            // as traverse->tag may be needed. It will get checked at
            // that point and deleted once traverse->tag is used,
            // later.
            if( checkresult == 0 )
              delete first;
          }
    
        // trajectory output via outputtime option
        // if( otime > 0.0 )
        //   {
        //     if( stime - ctime > otime )
        //       {
        //         ctime += otime;
        //         printf("Current State: Choice: %6.4f, Time: %6.6f\n",rchoice, ctime);
        //         complexList->printComplexList( 0 );
        //       }

        //   }

        // trajectory output via outputinterval option
        // if( ointerval >= 0 )
        //   {
        //     if( testBoolAttr(system_options, output_state) )
        //       {
        //         printf("Current State: Choice: %6.4f, Time: %6.6f\n",rchoice, stime);
        //         complexList->printComplexList( 0 );
        //       }
        //     pingAttr( system_options, increment_output_state );
        //   }

      } while( /*rate > 0.01 &&*/ stime < maxsimtime && checkresult == 0);
      
      // if( otime > 0.0 )
      //   printf("Final state reached: Time: %6.6f\n",stime);
      // if( ! (sMode & SIMULATION_MODE_FLAG_PYTHON ) )
      //   if( ointerval < 0 || testLongAttr(system_options, output_state ,=, 0 ))
      //     complexList->printComplexList( 0 );
      
      if( stime == NAN )
        printStatusLine(system_options, current_seed, STOPRESULT_NAN, 0.0, NULL);
      else if ( checkresult > 0 )
        {
          printStatusLine(system_options, current_seed, STOPRESULT_NORMAL, stime, traverse->tag );
          delete first;
        }
      else
        printStatusLine(system_options, current_seed, STOPRESULT_TIME, stime, NULL );
      
      // if( ! (sMode & SIMULATION_MODE_FLAG_PYTHON) )
      //   printf("Trajectory Completed\n");
    }
  else if( testLongAttr(system_options, trajectory_type ,=, 1 ) )
    {
      // begin transition times mode case
      int stopindex=-1;
      if( stopcount > 0 )
        {
          curcount = 0;
          first = getStopComplexList( system_options, 0 );
          traverse = first;
          while( curcount < stopcount && stopindex < 0 )
            {
              if( strstr( traverse->tag, "stop") != NULL )
                stopindex = curcount;
              curcount++;
              traverse = traverse->next;
            }
          delete first;
        }

      do {
        rchoice = (rate * random()/((double)RAND_MAX));
        stime += (log(1. / (double)((random() + 1)/(double)RAND_MAX)) / rate );
        complexList->doBasicChoice( rchoice, stime );
        rate = complexList->getTotalFlux();
        // if( ointerval >= 0 )
        //   {
        //     if( testBoolAttr(system_options, output_state) )
        //       complexList->printComplexList( 0 );
        //     pingAttr( system_options, increment_output_state );
        //   }
        if( stopcount > 0 )
          {
            curcount = 1;
            checkresult = 0;
            first = getStopComplexList( system_options, 0 );
            checkresult = complexList->checkStopComplexList( first->citem );
            traverse = first;
            while( curcount < stopcount && checkresult == 0 )
              {
                traverse = traverse->next;
                checkresult = complexList->checkStopComplexList( traverse->citem );
                curcount++;
              }
            if( checkresult == 0 )
              {
                delete first;
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
      if( checkresult > 0 )
        delete first;
      // if( ointerval < 0 || testLongAttr(system_options, output_state ,=, 0 ))
      //   complexList->printComplexList( 0 );

      if( stime == NAN )
        printStatusLine(system_options, current_seed, STOPRESULT_NAN, 0.0, NULL);
      else
        printStatusLine(system_options, current_seed, STOPRESULT_TIME, stime, NULL );
      
    }
}


void SimulationSystem::StartSimulation_First_Bimolecular( void )
{
  assert((simulation_mode & SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR)  &&  (simulation_mode & SIMULATION_MODE_FLAG_PYTHON));
  // Second part of this assert probably needs changing eventually.

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
  int checkresult = 0;

  double maxsimtime;
  long stopcount;
  long stopoptions;
  class stopcomplexes *traverse = NULL, *first = NULL;
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

  // if( ointerval >= 0 || otime_interval >= 0.0 )
  //   complexList->printComplexList( 0 );

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

  // if( ointerval >= 0 || otime_interval > 0.0 )
  //   {
  //     complexList->printComplexList( 0 );
  //     printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, ctime);
  //   }

  // Begin normal steps.
  rate = complexList->getTotalFlux();
  do {
    
    rchoice = rate * drand48();
    stime += (log( 1. / (1.0 - drand48()) ) / rate ); 
    
    // 1.0 - drand as drand returns in the [0.0, 1.0) range, 
    //  we need a (0.0,1.0] range for this to be the correct
    // distribution - log 1/U => U in (0.0,1.0] => 1/U => (+Inf,0.0] => 
    //    natural log:   (log(+Inf), 1.0]

    // if( !sMode && otime > 0.0 )
    //   {
    //     if( stime - ctime > otime )
    //       {
    //         ctime += otime;
    //         printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, ctime);
    //         complexList->printComplexList( 0 );
    //       }

    //   }

    complexList->doBasicChoice( rchoice, stime );
    rate = complexList->getTotalFlux();
    
    // if( !sMode && ointerval >= 0 )
    //   {
    //     if( testBoolAttr(system_options, output_state) )
    //       {
    //         complexList->printComplexList( 0 );
    //         printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, stime);
    //       }
    //     pingAttr( system_options, increment_output_state );
        
    //   }
    if( stopcount > 0 && stopoptions)
      {
        checkresult = 0;
        first = getStopComplexList( system_options, 0 );
        traverse = first;
        checkresult = complexList->checkStopComplexList( traverse->citem );
        while( traverse->next != NULL && checkresult == 0 )
          {
            traverse = traverse->next;
            checkresult = complexList->checkStopComplexList( traverse->citem );
          }
        if( checkresult == 0 && first != NULL)
          delete first;
      }
  } while( stime < maxsimtime && checkresult == 0);
      
  // if( !sMode && otime > 0.0 )
  //   {
  //     complexList->printComplexList(0);
  //     printf("Final state reached: Time: %6.6e\n",stime);
  //   }

  if( checkresult > 0 )
    {
      if( strcmp( traverse->tag, "REVERSE") == 0 )
        printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_REVERSE, stime, frate, traverse->tag );
      else 
        printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_FORWARD, stime, frate, traverse->tag );
      delete first;
    }
  else
    {
      printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_FTIME, stime, frate, NULL );
    }
  // if( ! sMode )
  //   printf("Trajectory Completed\n");
  
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
  PyObject *py_start_state = NULL, *py_complex = NULL;
  PyObject *py_seq = NULL, *py_struc = NULL;

  startState = NULL;
  if( complexList != NULL )
    delete complexList;

  complexList = new SComplexList( dnaEnergyModel );
  
  py_start_state = getListAttr(system_options, start_state);
  // new reference

  start_count = PyList_GET_SIZE(py_start_state);  
  // doesn't need reference counting for this size call.
  // the getlistattr call we decref later.
  
  for( int index = 0; index < start_count; index++ )
    {
      // #ifndef DEBUG_MACROS
      py_complex = PyList_GET_ITEM(py_start_state, index);
      // Borrowed reference, we do NOT decref it at end of loop.
            
      // #else
      //       py_complex = PyList_GetItem(py_start_state, index);
      // #endif

#ifdef DEBUG_MACROS
      printPyError_withLineNumber();
#endif

      sequence = getStringAttr(py_complex, sequence, py_seq);
      // new reference
      
      structure = getStringAttr(py_complex, structure, py_struc);
      // new reference
      
      id = getID_list( system_options, index );
      
      tempcomplex = new StrandComplex( sequence, structure, id );
      // StrandComplex does make its own copy of the seq/structure, so we can now decref.
      
      Py_DECREF( py_seq );
      Py_DECREF( py_struc );
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
