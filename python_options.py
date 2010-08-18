################################################################################
#                                                                              #
# Python implementation of the options object.                                 #
# Copyright blah blah blah 2010 Caltech                                        #
# Written by:  Joseph Schaeffer.                                               #
# Some stuff written by:  Chris Berlind                                        #
#                                                                              #
# (others add your names as you modify files)                                  #
#                                                                              #
################################################################################

class _OptionsConstants( object ):
    def __init__(self):
        self.ZERO_C_IN_K = 273.15
        """ 0 degrees Celsius, expressed in Kelvin."""
        
        pass

    @property
    def RATEMETHOD(self):
        return {"Invalid" :0, "Metropolis"    :1, \
                "Kawasaki":2,"EntropyEnthalpy":3}

    @property
    def DANGLES(self):
        return {"None":  0, "Some" : 1, \
                "All" :  2, "NupackDefault": 1}

    @property
    def ENERGYMODEL_TYPE(self):
        return {"Vienna":0, "Nupack":1, \
                "Others?":2}

    @property
    def SUBSTRATE_TYPE(self):
        return {"Invalid":0, "RNA":1, \
                "DNA":2}
    
    def __setattr__(self):
        pass

_OC = _OptionsConstants()
FILLIN = None

class MultistrandOptions( object ):
    def __init__(self):
        """ Creates an options object with default (presumably useful) values """

        ##################################################
        #                                                #
        # Data Members                                   #
        # ->Members new to the python implementation     #
        #                                                #
        #                                                #
        ##################################################

        self.errorlog = []
        """ Keeps lines relating to possible errors or warnings that
        should be reported to the user. Usually issues relating to the
        input file or parameters with odd values.

        TODO: implement some functions to report the errors found here.
        """
        

        #############################################
        #                                           #
        # Data Members: Energy Model                #
        # ->Members needed by the energy model      #
        #   NOTE: simple members only - custom      #
        #         members with different accessors  #
        #         appear after __init__.            #
        #                                           #
        #############################################
        
        self.gt_enable = False
        """ Allow GT base pairs? If not, penalize by 10000 kcal/mol.
        
        Type         Default
        int          False (0) : Do not allow GT base pairs.
        """

        self.log_ml = False
        """ Use logarithm to compute multiloop energy contributions?
        
        Type         Default
        int          False (0): Do not use the logarithm.

        If True, uses log to compute one component of the multiloop
        energy, for loops of length > 6. Otherwise, uses the usual
        linear approximation. Using the log formula is slightly more
        expensive as it makes computation time for a multiloop scale
        with the number of adjoining helices.
        """
        
        self.join_concentration = 1.0
        """ conc for V calcs
        
        Type         Default
        double       1.0 (M)

        Units are in M (molar), and represent the concentration of a
        single unique strand in the system. The volume simulated is
        then chosen using this parameter.
        """

        ###
        ### See the temperature property way below (after __init__)
        ### for more info on accessors for these data members.
        ###
        self._temperature_celsius = 37.0
        self._temperature_kelvin  = 310.15

        
        self.unimolecular_scaling = 1.6e6
        """ Rate scaling factor for unimolecular reactions.

        Type         Default
        double       1.6e6:
                     Unitless. Details on default in thesis."""
        
        self.bimolecular_scaling = 0.5e6
        """ Rate scaling factor for bimolecular reactions.
        
        Type         Default
        double       0.5e6:
                     Unitless. Details on default in thesis."""


        self.rate_method = _OC.RATEMETHOD['Kawasaki']
        """ Choice of methods for determining forward/reverse rates.
        
        Type         Default
        int          2: Kawasaki

        Should use the values in the _OC.RATEMETHOD dictionary rather
        than the numbers directly, as those should be consistent with
        that used in the energymodel.h headers.
        """

        self.dangles = _OC.DANGLES['NupackDefault']
        """ Dangles options for the energy model.
        
        Type         Default
        int          1: NupackDefault

        None [0]: Do not include any dangles terms in the energy model.
        Some [1]: Some dangles terms.  (Nupack Default)
        All  [2]: Include all dangles terms, including odd overlapping ones.

        See notes elsewhere re: what each term does.
        Should use the values in the _OC.DANGLES dictionary rather
        than the numbers directly, as those should be consistent with
        that used in the energymodel.h headers.
        """

        self.parameter_type = _OC.ENERGYMODEL_TYPE['Nupack']
        """ Which type of energy model parameter file to use.
        
        Type         Default
        int          1: Nupack

        Vienna [0]: No longer well tested. Recommend not using.
        Nupack [1]: Includes some multi-complex parameters, otherwise
                    nearly the same as mfold style files.

        Should use the values in the _OC.ENERGYMODEL_TYPE dictionary rather
        than the numbers directly, as those should be consistent with
        the ones defined in the options_python.h headers.
        """

        self.substrate_type = _OC.SUBSTRATE_TYPE['DNA']
        """ What substrate's parameter files to use. Note that some
        combinations of this and energymodel_type will be invalid and
        Multistrand will complain.
        
        Type         Default
        int          2: DNA

        Invalid [0]: Indicates we should not auto-search for a param file.
        RNA     [1]: RNA parameters are available for Vienna and Nupack, see also
                     the comment in parameter_file_version.
        DNA     [2]: DNA parameters are publicly available via the Nupack distribution,
                     and possibly by the Vienna group upon request.

        Should use the values in the _OC.SUBSTRATE_TYPE dictionary rather
        than the numbers directly, as those should be consistent with
        the ones defined in the options_python.h headers.
        """
        
        self.parameter_file = None
         """ Shortcut for using a very specific parameter file. Usually shouldn't be used.
        
        Type         Default
        str          None

        Should only be set to a string if it's ok to actually search
        for that parameter file. None will pass an error back to
        Multistrand if it gets used.
        """

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


        self.trajectory_type = 0
        """ The type of trajectory data we are concerned with.
        
        Type         Default
        int          0: Normal
        
        Normal          [0]: First passage time only.
        Transition time [1]: Transition time data kept for some set of states.
        """

        self.simulation_mode = 0
        """ The simulation mode.
        
        Type         Default
        int          0: Normal
        
        Normal            [0]:
        First step        [1]:
        Python normal     [2]:
        Python first step [3]:
        """
        
        self.simulation_time = 10000.0
        """ Maximum time (in seconds) allowed for each trajectory.
        
        Type         Default
        double       10000.0
        """
        
        self.num_simulations = 1
        """ Total number of trajectories to run. 
        
        Type         Default
        int          1
        """
        
        self.initial_seed = None
        """ Initial random number seed to use.
        
        Type         Default
        long         None
        
        If None when simulation starts, a seed will be chosen by the RNG however
        it feels like.
        """
        
        self.unique_id = 0
        """ Unique identifier for strands.
        
        Type         Default
        int          0
        
        Incremented for each strand added.
        """
        
        self.complexcount = FILLIN
        """ count of complexes in initial state """
        ##int
        
        self.strandcount = FILLIN
        """ total # of strands """
        ##int

        self.energymode = FILLIN
        """ deprecated, energy info """
        ##int
        ('EnergyMode', 'energy_mode')
        
        self.flagarray = FILLIN
        """ deprecated, bit position info """
        ##int
        
        ####################
        #
        # BEGIN startstop
        #
        ####################
        
        self.start_complexes = []
        """ The start states, i.e. a list of Complex objects.
        
        Type         Default
        list         []
        
        Start states should be added to this list (e.g. by the parser) so
        trajectories know how to start.
        """         
        
        self.stop_conditions = []
        """ The stop states, i.e. a list of StopCondition objects.
        
        Type         Default
        list         []
        
        Stop states should be added to this list (e.g. by the parser) so
        trajectories know when to end.
        """ 
        
        self.use_stop_states = None
        """ Indicates whether trajectories should end when stop states
        are reached.
        
        Type            Default
        boolean         None
        
        Set to True by default if a stop state is defined. Can be set to False
        manually to avoid stopping at defined stop states. Can only be manually
        changed back to True from False.
        """
        
        self.stopcount = 0
        """ The number of stop states. Equivalent to 'len(self.stop_conditions)'.
        
        Type         Default
        int          0
        
        Incremented automatically when a stop state is added. Should not be 
        modified externally.
        """
        
        self.output_time = None
        """ The amount of time (in seconds) to wait between outputs of 
        trajectory information.
        
        Type         Default
        float        None
        
        A value of None corresponds to not basing outputs on output_time
        (but perhaps outputting based on some other condition). A value of 0 
        means output as often as possible.
        """
        
        self.output_interval = None
        """ The number of states between outputs of trajectory information.
        
        Type         Default
        int          None
        
        A value of None corresponds to not basing outputs on output_interval
        (but perhaps outputting based on some other condition). A value of 0 
        means output every state, 1 means every other state, and so on.
        """
        
        self.current_interval = 0
        """ Current value of output state counter.
        
        Type         Default
        int          0
        
        When current_interval is equal to output_interval, the output state is 
        True, and otherwise the output state is False. This is modified by 
        increment_output_state, and probably shouldn't be used externally."""
        
        self.output_state = False
        """ Indicates whether output should be reported.
        
        Type         Default
        boolean      False
        
        Value should be True if self.current_interval == self.output_interval 
        and False otherwise.        
        """
        

        ####################
        #
        # BEGIN pythondata
        #
        ####################
        
        # TODO: change some of these to use properties and possibly rename
        
        self.python_trajectory_time = None
        """ The total elapsed time of the most recently completed trajectory.
        
        Type         Default
        double       None
        
        Set by ssystem.cc when a trajectory completes.
        """
        
        self.python_trajectory_tag = None
        """ The tag of the stop state of the most recently completed trajectory.
        
        Type         Default
        string       None
        
        Set by ssystem.cc when a trajectory completes.
        """
        
        self.python_collision_rate = None
        """ The collision rate of the most recently completed trajectory.
        
        Type         Default
        double       None
        
        Set by ssystem.cc when a trajectory completes.
        """
        
        self.python_current_seed = None
        """ The seed used by the random number generator in the most recently
        completed trajectory.
        
        Type         Default
        long         None
        
        Set by MultistrandOptions at the end of each trajectory.
        """
        
        self.python_current_time = 0.0
        """ The current elapsed time of the currently running trajectory.
        
        Type         Default
        double       0.0
        
        Set by ssystem.cc at certain points during the simulation.
        """
        
        self.python_trajectory_completion_flag = False
        """ Indicates whether the trajectory has completed.
        
        Type         Default
        boolean      False
        
        Should be set by MultistrandOptions when a trajectory completes.
        """
        
        self.python_halt_trajectory_flag = False
        """ Indicates whether Multistrand was told to halt its trajectory by an
        external Python program.
        
        Type         Default
        boolean      False
        
        Read by ssystem.cc. There should probably be a function that sets it.
        """
        
        self.python_suspend_trajectory_flag = False
        """ Indicates whether Multistrand was told to suspend its trajectory by 
        an external Python program.
        
        Type         Default
        boolean      False
        
        Read by ssystem.cc. There should probably be a function that sets it.
        """
        
        
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
    
    
    def increment_output_state(self):
            """ Modifies self.current_interval and self.output_state as 
            necessary based on self.output_interval.
            """
            if self.output_interval == None or self.output_interval < 0:
                raise ValueError("output_interval has invalid value: %s" % self.output_interval)
            
            elif self.current_interval > self.output_interval:
                raise ValueError("current_interval has invalid value: %s" % self.current_interval)
            
            elif self.current_interval == self.output_interval:
                self.current_interval == 0
            
            else:
                self.current_interval += 1
            
            self.output_state = (self.current_interval == self.output_interval)
    
    
    @property
    def temperature(self):
        return self._temperature_kelvin
    """ Temperature, in degrees Kelvin.
    
    Type         Default
    double          310.15 (K)
    
    Standard units of degrees K, and default is the usual 37(C)
    that's the base value in all the parameter files. Multistrand
    uses Kelvin internally so we use it here as well.  used for
    scaling parameters and those are always relative to degrees K.
    """

    @temperature.setter
    def temperature(self,val):
        """ performs some sanity checking and provides a log error if
        perhaps the user was confused.

        Assumptions: Input should be in Kelvin, if it's not in a
        'reasonable' range for Kelvin, convert to Celsius if it's in a
        reasonable range for C [sending an output to the error log
        warning that it did so], otherwise error.

        Reasonable ranges:
            [0,100]  : Celsius
            [273,373]: Kelvin
            Others:    If you want a Fahrenheit reasonable range, I think you
                       might be unreasonable. Also, it overlaps with Celsius a bit much.

<<<<<<< local
            Yes, these ranges are quite generous.
        """
        if 273.0 < val < 373.0:
            self._temperature_kelvin = val
            self._temperature_celsius = val - _OC.ZERO_C_IN_K
        elif 0.0 < val < 100.0:
            self._temperature_celsius = val
            self._temperature_kelvin = val + _OC.ZERO_C_IN_K
            self.errorlog.append("Warning: Temperature was set at the value [{0}]. We expected a value in Kelvin, or with appropriate units.\n         Temperature was automatically converted to [{1}] degrees Kelvin.\n".format(val, self._temperature_kelvin))


