/*
   Copyright (c) 2007-2008 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/
 
#ifndef __OPTIONS_H__n
#define __OPTIONS_H__

#include "optionlists.h"
#include <string>

#define OPTION_ENERGYMODEL     0x0001
#define OPTION_LOGFILE         0x0002
#define OPTION_STRANDCOUNT     0x0004
#define OPTION_STRANDLIST      0x0008
#define OPTION_COMPLEXCOUNT    0x0010
#define OPTION_COMPLEXLIST     0x0020
#define OPTION_SIMULATIONTIME  0x0040
#define OPTION_TRAJECTORYCOUNT 0x0080

// Used for setEnergyModel
#define VIENNADNA     1
#define NUPACKDNA23   2
#define NUPACKRNA23   3

#define STOPTYPE_STRUCTURE                     0 
#define STOPTYPE_BOUND                         1
#define STOPTYPE_DISASSOC                      2
#define STOPTYPE_LOOSE_STRUCTURE               3
#define STOPTYPE_PERCENT_OR_COUNT_STRUCTURE    4




#define STOPCONDITION_NORMAL           1
#define STOPCONDITION_REVERSE          2
#define STOPCONDITION_TIME            -1
#define STOPCONDITION_FORWARD          3
#define STOPCONDITION_ERROR           -2



/* 

   class Options

   This class contains all of the appropriate input/output options.
   These options are read in from a file, or set via the command line (which
   overrides any settings occuring in the input file). This class is used
   by the Simulation System to initialize new trajectories and to handle
   output events. Most of the data members are accessible by standard
   accessor functions, though a few have a more specialized interface.

   There are a few auxilary structures defined for use in the class, as well
   as for data passing back to the simulation system.

*/



class Options
{
 public:
  /* The default constructor is the standard one which we will use, as 
     most of the input takes place after reading the input file name off
     the command line. */
  Options( void );
  Options( char *filename );

  
  /* Default destructor. */
  ~Options( void );

  /* Initialization functions */
  int readCommandLine( int argc, char ** argv );
  int loadInputFile( void );
  void finalizeInput( void );

  /* Output Function */

  void printStatusLine( long r_seed, char *stop_id, double end_time );
  void printStatusLine( long r_seed );
  void printStatusLine_First_Bimolecular( long r_seed, int completiontype, double completiontime, double forwardrate, char *tag );
  void printStatusLine_Final_First_Bimolecular( double total_rate, double *total_times, long *total_types, long total_count, double *computed_rate_means, double *computed_rate_mean_diff_squared );
  void printStatusLine_Warning( int warnid, int data );
  void printTrajLine( char *trigger_id, double cur_time );
  FILE *getOutputDescriptor( void );
  void closeOutput( void );

  /* Accessors */  
  int addSequence( char *id, char *seq );
  int addSequence( char *id, char *seq, double conc );
  int addSequence( char *id, char *seq, int count   );

  int addStopStructures( class stopcomplexes *stopstructs, int errortype = 0 );
  int addStartStructure( class complex_item *startstruct, int errortype = 0);

  /* Accessors for Python Interface */
  
  int addSequence_Python(  std::string id, std::string seq );

  //... plus other versions for later implementation
  
  // list of strands interface
  //int beginStrandList_Python( void );
  int addStrandToIdentList_Python( std::string name);  // define order as right to left
  int clearIdentList_Python( void );

  int beginStartStructure_Python( void ); 

    // uses implied strand list already created
  int addStartStructureComplex_Python( std::string struc, int type = STOPTYPE_STRUCTURE );

  int finalizeStartStructure_Python( void );

  int beginStopStructure_Python( std::string tag ); 

    // uses implied strand list already created
  int addStopStructureComplex_Python( std::string struc , int type = STOPTYPE_STRUCTURE, int count = 0);


  //int completeStopStructure_Python( void );
  int finalizeStopStructures_Python( void );
  // ERROR CODES FOR FINALIZE STEPS:

  // 0  : no stop/start list created.
  // -1 : Strand ID not found.
  // -2 : strand length and structure length mismatch
  // -3 : strand break not found/not valid
  // -4 : structure not connected
  // -5 : invalid base pairing in structure
  // -6 : mismatched paren, too many )'s
  // -7 : mismatched paren, too many ('s


  void setLogfile_Python( std::string newlogfilename );
  void setTrajectoryfile_Python( std::string newtrajfilename );
  void setParameterFile_Python( std::string newparamfile );

  // simulation time and output methods:

  void setCurSimTime( double cur_time ); // used internally by SimulationSystem object
  double checkCurrentTime_Python(void ); 
  int checkCompleted_Python( void );
  void resetCompleted_Python( void ); // used in start simulation to reset completion in the options object.
  std::string getTrajectoryTag_Python( void );
  double getTrajectoryTime_Python( void );

  void haltCurrentTrajectory_Python( void );
  void suspendCurrentTrajectory_Python( void );
  void resumeCurrentTrajectory_Python( void );

  // directly related to the above function and the python simulation mode, but
  // not part of the direct interface as it's used internally.

  int getCurrentTrajectoryHaltFlag( void );
  int getCurrentTrajectorySuspendFlag( void );

  void setCollisionRate_Python( double newrate );
  double getCollisionRate_Python( void );
  long getCurrentSeed_Python( void );
  
  /* End Python Interface */


  void setSimulationTime( double newtime );
  void setDangles( int dangleopt );
  void setStopOptions( int newoptions );
  void setLogfile( char *newlogfilename );
  void setTrajectoryfile( char *newtrajfilename );
  void setTrajectoryType( int ttype );
  void setSimulationMode( int mode ); /* See mode listing above */
  void setTemperature( double temp );
  void setNumSimulations( int newnumsims );
  void setOutputInterval( int newoutputinterval );
  void setOutputTime( double newoutputtime );
  void setOutputComplexEnergyOptions( int complexoutputoptions );
  void setInitialSeed( long seedval );
  void setParameterFile( char *newparamfile );
  void setParameterType( int type );
  void setEnergyModel( int modeltype );
  void setEnergyMode( int mode );
  void setIntermolecularScaling( double rate );
  void setIntramolecularScaling( double rate );
  void setJoinRate( double newjoinrate );
  void setJoinRateByConcentration( double concentration );
  void setJoinRateByVolume( double volume ); 
  void setJoinRateByEnergy( double energy );
  void setJoinConcentration( double concentration );
  void setRateMethod( int method);
  void setGTWobble( int gtstate );
  void setLogML( int newlogml );
  void incrementOutputState( void );

  long getInitialSeed( void );
  char *getSequenceByID( char *id );
  char *getStructure( int complex_id );
  char *getBoltzmannStructure( int complex_id );
  char *getSequence( int complex_id );
  class identlist *getID_list( int complex_id );

  char *getParameterFile( void );
  int getParameterType( void );
  int getEnergyModel( void );
  double getTemperature( void );

  int getOutputInterval( void );
  double getOutputTime( void );
  int getOutputState( void );
  int getOutputComplexEnergyOptions( void );
  int getNumSimulations( void );
  double getSimulationTime( void );
  //  double getJoinRate( void );
  // double getJoinEnergyExtra( void );
  double getJoinConcentration( void );
  double getIntramolecularScaling( void );
  double getIntermolecularScaling( void );
  int getDangles( void );
  int getGTWobble( void );
  int getLogML( void );
  int getTrajectoryType( void );
  int getSimulationMode( void );
  int getStopCount( void );
  int getStopOptions( void );
  int getEnergyMode( void );  
  int getRateMethod( void );
  class stopcomplexes *getStopComplexList( int index );

 private:
  /* Data members */
  /* Control data */
  int flagarray;  // Contains bit positions for each option type to indicate 
                  // whether those options have been set via command line, and
                  // thus overriding the same option in the input file.

  /* Options Data */
  char energymodelfilename[80]; // Energy model paramater file. Defaults to
                                // "dna.par", the vienna parameter file.

  int energymodel_type;   // type of the parameter file. 0 for dna.par (Vienna)
                          // 1 for mfold parameter dG filename. (Mfold)

  int energymodel;        // flag for automated searching for particular energy
                          // files. 

  char inputfilename[80]; // this is the only piece of options data which
                          // cannot be set from the input file. Duh.

  char logfilename[80];   // Filename to which log entries should be written.
                          // Stored for posterity, as the next option is the
                          // one used for most log file entries.


  FILE *logfile;          // Logfile file pointer.

  char trajectoryfilename[80]; // Filename for trajectory output.
  FILE *trajfile;
  int trajtype;

  int simulationmode;     // See simulation mode defines above.
                          // Currently covers normal simulation, and the 
                          // fixed first step simulation mode.

  int energymode;         // Set to 1 for energy only, 0 otherwise.

  int strandcount;        // Total number of strands in the system.
  class strandlist *strands; // A list of all strands, containing sequence
                              // and naming information.

  int complexcount;       // Total number of complexes in the initial state
                          // of the system.

  int gtenable;           // For nupack parameter files, enables GT pairings (by default, penalizes 100000 per GT pair inside an multiloop or open loop).

  int logml;              // Use logarithmic multiloop size penalties.
                          // Nupack defaults to no, vienna to yes. 
                          // This parameter should be -1 for default, 
                          // and 0 or 1 to override the defaults.

  class complex_item *start_structure;  // A list of all initial complexes, contains
                                  // naming, sequence and structure information
  class complex_item **start_list; // flat list of initial complexes.
  int start_count;

  double simulationtime;  // Max time units to run the simulation for.
  
  //  double joinenergy;      // Join energy. (without kT multiplier)
  //double joinrate;        // Join rate for complex moves.
  double joinconc;        // effective strand conc for deriving join rate
  double temperature;

  double intramolecularscaling; // Scaling factor for reactions between bases within the same complex.
  
  double intermolecularscaling; // Scaling factor for reactions between multiple complexes.

  int trajectorycount;    // Number of trajectories to run.
  int ratemethod;         // rate model to use
  int dangles;            // dangle option for energy model
  int outputinterval;     // Controls the frequency of individual state output.
  int currentinterval;
  double outputtime;      // alternate way to control trajectory output
  int stopoptions;        // Flag variable for which type of stopping condition
                          // should be used with this simulation.

  int stopcount;          // Number of Stop States in the system.

  class stopcomplexes *stoplist;   // A list of all complex names in the
                                    // stop state.


  long initialseed;     // initial random number seed to use. 0 = none.
  
  long unique_id;       // Stores the ID of the next strand to be input into the system. Increments from 0.

  /* Auxilary Functions */

  void sanitize( char *buf, int *length ); // Cleans up log file lines so that
                                          // parsing is easier. 


  void readStrandList( FILE *fp, char *buffer );
  void readComplexList( FILE *fp, char *buffer );


  /* Private Data Members for Python Interface functions */
  class identlist *python_identlist;

  class complex_item *python_start_structure;
  class stopcomplexes *python_stop_structure;
  class strandlist *python_strands;

  double python_current_time;
  int python_trajectory_completion_flag;
  std::string python_trajectory_tag;
  double python_trajectory_time;
  int python_halt_trajectory_flag;
  int python_suspend_trajectory_flag;
  double python_k_collision;
  long python_current_seed;


};

#endif
