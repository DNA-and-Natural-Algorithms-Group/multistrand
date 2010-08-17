################################################################################
#                                                                              #
# Python implementation of the options object.                                 #
# Copyright blah blah blah 2010 Caltech                                        #
# Written by:  Joseph Schaeffer.                                               #
#                                                                              #
# (others add your names as you modify files)                                  #
#                                                                              #
################################################################################

class MultistrandOptions( object ):
    def __init__(self):
        """ Creates an options object with default (presumably useful) values """


        ####################
        #
        # BEGIN energymodel
        #
        ####################
		self.gtenable = FILLIN
		""" GT allowed or not?"""
		##int
        ###  needs to be used in energy model files.

		self.logml = FILLIN
		""" log ml for loops """
		##int
		('LogML', 'log_ml')

		self.joinconc = FILLIN
		""" conc for V calcs """
		##double
        ('JoinConcentration', 'join_concentration')

		self.temperature = FILLIN
		""" Temperature, in degrees C """
		##double
		('Temperature', 'temperature')

		self.intramolecularscaling = FILLIN
		""" scale factors, rename to:  uni-molecular-scaling """
		##double
        ('IntramolecularScaling','intramolecular_scaling')

		self.intermolecularscaling = FILLIN
		""" rename to: bi-molecular-scaling """
		##double
        ('IntermolecularScaling','intermolecular_scaling')

		self.ratemethod = FILLIN
		""" kawasaki / metropolis """
		##int
		('RateMethod', 'rate_method')        

		self.dangles = FILLIN
		""" dangle options """
		##int
		('Dangles', 'dangles')        

		self.energymodel_type = FILLIN
		""" nupack / vienna / etc enumeration """
		##int
		('ParameterType', 'parameter_type')

		self.energymodel = FILLIN
		""" automated search path or specific file  flag"""
		##int
		('EnergyModel', 'energy_model')        

		self.energymodelfilename = FILLIN
		""" energy model parameter file path/name """
		##char
		('ParameterFile', 'parameter_file')        


        ####################
        #
        # BEGIN simmode
        #
        ####################
        
		self.inputfilename = FILLIN
		""" deprecated """
		##char


		self.logfilename = FILLIN
		""" filename to send log entries to """
		##char


		self.trajectoryfilename = FILLIN
		""" filename to send trajectory entries to """
		##char


		self.trajtype = FILLIN
		""" 0 == normal, 1 == transition time data  """
		##int
		('TrajectoryType', 'trajectory_type')

		self.simulationmode = FILLIN
		""" normal / first step mode / python """
		##int
		('SimulationMode', 'simulation_mode')

		self.energymode = FILLIN
		""" deprecated, energy info """
		##int
		('EnergyMode', 'energy_mode')

		self.strandcount = FILLIN
		""" total # of strands """
		##int


		self.simulationtime = FILLIN
		""" max sim time """
		##double
		('SimulationTime', 'simulation_time')

		self.trajectorycount = FILLIN
		""" total # of trajectories to run. """
		##int
		('NumSimulations', 'num_simulations')

		self.initialseed = FILLIN
		""" initial random # seed to use """
		##long
		('InitialSeed', 'initial_seed')

		self.unique_id = FILLIN
		""" unique identifier, incremented for each strand added. """
		##long


		self.flagarray = FILLIN
		""" deprecated, bit position info """
		##int


		self.complexcount = FILLIN
		""" count of complexes in initial state """
		##int

        ####################
        #
        # BEGIN startstop
        #
        ####################
		self.start_count = FILLIN
		""" random? """
		##int

		self.stopoptions = FILLIN
		""" 2 = use stop structures """
		##int
		('StopOptions', 'stop_options')

		self.stopcount = FILLIN
		""" random? """
		##int
		('StopCount', 'stop_count')
        
		self.outputinterval = FILLIN
		""" number of states between output counts  """
		##int
        ('OutputInterval', 'output_interval')

		self.currentinterval = FILLIN
		""" current state interval """
		##int
        #### vaguely related to getOutputState and stuff like that.
        ('OutputState', 'output_state')

		self.outputtime = FILLIN
		""" time between output counts """
		##double
		('OutputTime', 'output_time')

		self.complex_item = FILLIN
		None
		##class

		self.stopcomplexes = FILLIN
		None
		##class
        ('StopComplexList', 'stop_complex_list')

		self.strandlist = FILLIN
		None
		##class

        ####################
        #
        # BEGIN pythondata
        #
        ####################
		self.identlist = FILLIN
		None
		##class


		self.complex_item = FILLIN
		None
		##class


		self.stopcomplexes = FILLIN
		None
		##class


		self.strandlist = FILLIN
		None
		##class


		self.python_current_time = FILLIN
		None
		##double


		self.python_trajectory_completion_flag = FILLIN
		None
		##int


		self.python_trajectory_time = FILLIN
		None
		##double


		self.python_halt_trajectory_flag = FILLIN
		None
		##int
        ('CurrentTrajectoryHaltFlag', 'current_trajectory_halt_flag')


		self.python_suspend_trajectory_flag = FILLIN
		None
		##int
        ('CurrentTrajectorySuspendFlag', 'current_trajectory_suspend_flag')


		self.python_k_collision = FILLIN
		None
		##double


		self.python_current_seed = FILLIN
		None
		##long

        ##############################
        #
        # functions used
        #
        ##############################
		('BoltzmannStructure', 'boltzmann_structure')
		('Sequence', 'sequence')
		('Structure', 'structure')


		('finalizeInput', 'finalize_input')
		('incrementOutputState', 'increment_output_state')
		('printStatusLine', 'print_status_line')
		('printStatusLine_Final_First_Bimolecular', 'print_status_line__final__first__bimolecular')
		('printStatusLine_First_Bimolecular', 'print_status_line__first__bimolecular')
		('printStatusLine_Warning', 'print_status_line__warning')
		('printTrajLine', 'print_traj_line')
		('resetCompleted', 'reset_completed')
		('resetCompleted_Python', 'reset_completed__python')
		('setCollisionRate_Python', 'set_collision_rate__python')
		('setCurSimTime', 'set_cur_sim_time')






