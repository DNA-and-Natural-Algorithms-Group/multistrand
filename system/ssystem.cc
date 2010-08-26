/*
  Copyright (c) 2007-2010 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#include "python_options.h"
#include "ssystem.h"

#include <string.h>
#include <time.h>
#include <stdlib.h>

//#include <math.h>
//#include <assert.h>
//#include <unistd.h>


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


//#define SRANDOMDEV
// uncomment above line if you have srandomdev() available for better
// random number generation. 

void SimulationSystem::StartSimulation( void )
{
  FILE *fp; // used only for initial random number generation.
  int curcount = 0;
  long ointerval;

  getLongAttr(system_options, output_interval,&ointerval);

  //#ifndef SRANDOMDEV
  //  srandom( time( NULL) );

  if( simulation_mode & SIMULATION_MODE_FLAG_PYTHON )
    {
      pingAttr( system_options, interface_reset_completion_flag );
      // replaces: 
      //   callFunc_NoArgsToNone(system_options, reset_completed_python);
      // by making reset_completion_flag be a @property descriptor on options object.
      
      // no longer needed: reset_completion_flag clears the rate value.
      //callFunc_DoubleToNone(system_options, set_python_collision_rate, -1.0);
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
      // if( ointerval < 0 && !(simulation_mode & SIMULATION_MODE_FLAG_PYTHON))
      //   printf("Seed: 0x%lx\n",current_seed);
      InitializeSystem();

      if( ointerval < 0 && !(simulation_mode & SIMULATION_MODE_FLAG_PYTHON))
        {
          printf("System %d Initialized\n",curcount);
        }
      //if( getLongAttr(system_options, trajectory_type) > 0 )
      //m_printTrajLine(system_options, NULL, curcount );
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
        m_printStatusLine(system_options, current_seed, "ERROR", stime );
      else if ( checkresult > 0 )
        m_printStatusLine(system_options, current_seed, traverse->tag, stime );
      else
        m_printStatusLine(system_options, current_seed, "INCOMPLETE", stime );
      
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

      m_printTrajLine(system_options,"Start",0);
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
            m_printTrajLine(system_options, traverse->tag, stime );
          }
        else
          m_printTrajLine(system_options,"NOSTATE", stime );
      } while( /*rate > 0.01 && */ stime < maxsimtime && !(checkresult > 0 && stopindex == curcount-1));
      
      if( ointerval < 0 || testLongAttr(system_options, output_state ,=, 0 ))
        complexList->printComplexList( 0 );

      if( stime == NAN )
        m_printStatusLine(system_options, current_seed, "ERROR", stime );
      else
        m_printStatusLine(system_options, current_seed, "INCOMPLETE", stime );
      
      printf("Trajectory Completed\n");
    }
}


void SimulationSystem::StartSimulation_First_Bimolecular( void )
{
  int curcount = 0;
  long random_seed = 0;
  long ointerval = -1;
  long initial_seed = 0;

  getLongAttr(system_options, output_interval,&ointerval);

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

  int sMode = simulation_mode & SIMULATION_MODE_FLAG_PYTHON;

  assert( simulation_mode & SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR );

  InitializeRNG();

  while( simulation_count_remaining > 0 )
    {
      setLongAttr( system_options, interface_current_seed, current_seed );
      InitializeSystem();

      //      if( getLongAttr(system_options, trajectory_type) > 0 )
      //    m_printTrajLine(system_options, NULL, curcount );
      //      if( sMode() )
      //        callFunc_NoArgsToNone(system_options, reset_completed);


      SimulationLoop_First_Bimolecular( &completiontime, &completiontype, &forwardrate, &tag );
      if (tag != NULL )
        setStringAttr( system_options, interface_current_tag, tag );

      setLongAttr( system_options, interface_current_completion_type, completiontype );
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

      m_printStatusLine_First_Bimolecular(system_options, current_seed, completiontype, completiontime, forwardrate, tag );
      generateNextRandom();
      pingAttr( system_options, increment_trajectory_count );
      simulation_count_remaining--;
    }
  if( !sMode )
    m_printStatusLine_Final_First_Bimolecular(system_options, total_rate, total_time, total_types, curcount, computed_rate_means, computed_rate_mean_diff_squared );
  if( total_types[2] > 0 && !sMode )
    m_printStatusLine_Warning(system_options, 0, total_types[2] );

}

/*




 */


void SimulationSystem::SimulationLoop_First_Bimolecular( double *completiontime, int *completiontype, double *frate, char **tag )
{
  double rchoice,rate,stime=0.0;
  int curcount = 0;
  int checkresult = 0;
  double ctime = 0.0;

  double maxsimtime;
  long stopcount;
  long stopoptions;
  class stopcomplexes *traverse;
  long ointerval;
  int sMode = simulation_mode & SIMULATION_MODE_FLAG_PYTHON;
  long trajMode;
  double otime;
  double otime_interval;

  getLongAttr(system_options, trajectory_type,&trajMode);
  getLongAttr(system_options, output_interval,&ointerval);
  getDoubleAttr(system_options, output_time,&otime);
  getLongAttr(system_options, use_stop_conditions,&stopoptions);
  getLongAttr(system_options, stop_count,&stopcount);
  getDoubleAttr(system_options, simulation_time,&maxsimtime);
  getDoubleAttr(system_options, output_time, &otime_interval );
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

  if( ointerval >= 0 || otime_interval >= 0.0 )
    complexList->printComplexList( 0 );
  if( ointerval >= 0 )
    printf("Initial State (before join).\n",rchoice, ctime);

  complexList->doJoinChoice( rchoice );

  *frate = rate * dnaEnergyModel->getJoinRate_NoVolumeTerm() / dnaEnergyModel->getJoinRate() ; // store the forward rate used for the initial step.

  if( ointerval >= 0 || otime_interval > 0.0 )
    {
      complexList->printComplexList( 0 );
      printf("Current State: Choice: %6.4f, Time: %6.6e\n",rchoice, ctime);
    }

  if( sMode )
    setDoubleAttr( system_options, interface_collision_rate, *frate );
 
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
      setDoubleAttr(system_options, interface_current_time, stime);
      

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
        // if external python interface has asked us to complete.
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
  // if( ointerval >= 0 )
  //printf("Final state reached: Time: %6.6e\n",stime);

    
  /*   if( stime == NAN )
       m_printStatusLine(system_options, current_seed, "ERROR", stime );
       else if ( checkresult > 0 )
       m_printStatusLine(system_options, current_seed, traverse->tag, stime );
       else
       m_printStatusLine(system_options, current_seed, "INCOMPLETE", stime );
  */
  // printing is handled at the upper level for the status information. We do, however, need to return the info.

  if( checkresult > 0 )
    {
      if( strcmp( traverse->tag, "REVERSE") == 0 )
        *completiontype = STOPCONDITION_REVERSE;
      else 
        {
          *completiontype = STOPCONDITION_FORWARD;
        }
      *tag = traverse->tag;
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
      
      //simulation_mode == SIMULATION_MODE_FIRST_BIMOLECULAR )
      //      if ( boltzmann_sampling )
      //        structure = getStringAttr(py_complex, boltzmann_structure, py_struc);
      //      else
      // boltzmann sampling is not our problem - the options object
      // should be providing structures according to the boltzmann or
      // not as it decides, not our job now. :)

      structure = getStringAttr(py_complex, structure, py_struc);
      
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
      //      current_seed = lrand48(); // grab a seed to use as initial seed.
    }
  // now initialize this generator using our random seed, so that we can reproduce as necessary.
  srand48( current_seed );
}

void SimulationSystem::generateNextRandom( void )
{
  current_seed = lrand48();
  srand48( current_seed );
}
