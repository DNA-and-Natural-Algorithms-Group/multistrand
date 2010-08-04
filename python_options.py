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


        self.logml = FILLIN
		""" log ml for loops """
		##int
		('LogML', 'log_m_l')
        
    
        self.joinconc = FILLIN
		""" conc for V calcs """
		##double
        
        
        self.temperature = FILLIN
		None
		##double
		('Temperature', 'temperature')

        self.intramolecularscaling = FILLIN
		""" scale factors, rename to:  uni-molecular-scaling """
		##double
        

        self.intermolecularscaling = FILLIN
		""" rename to: bi-molecular-scaling """
		##double
        
        
        self.ratemethod = FILLIN
		""" kawasaki / metropolis """
		##int
		('RateMethod', 'rate_method')        
        
        self.dangles = FILLIN
		""" dangle options """
		##int
		('Dangles', 'dangles')        
        
        self.energymodel_type = FILLIN
		None
		##int
        
        
        self.energymodel = FILLIN
		None
		##int
		('EnergyModel', 'energy_model')        
        
        self.energymodelfilename = FILLIN
		None
		##char
        

        ##############################
        #
        # functions used
        #
        ##############################
		('BoltzmannStructure', 'boltzmann_structure')
		('CurrentTrajectoryHaltFlag', 'current_trajectory_halt_flag')
		('CurrentTrajectorySuspendFlag', 'current_trajectory_suspend_flag')
		('EnergyMode', 'energy_mode')
		('InitialSeed', 'initial_seed')

		('NumSimulations', 'num_simulations')
		('OutputInterval', 'output_interval')
		('OutputState', 'output_state')
		('OutputTime', 'output_time')
		('ParameterFile', 'parameter_file')
		('ParameterType', 'parameter_type')

		('Sequence', 'sequence')
		('SimulationMode', 'simulation_mode')
		('SimulationTime', 'simulation_time')
		('StopComplexList', 'stop_complex_list')
		('StopCount', 'stop_count')
		('StopOptions', 'stop_options')
		('Structure', 'structure')

		('TrajectoryType', 'trajectory_type')
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






