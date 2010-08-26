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

from options_interface import Interface
from python_objects import Strand, Complex, RestingState, StopCondition
import copy

class _OptionsConstants( object ):
    def __init__(self):
        self.ZERO_C_IN_K = 273.15
        pass
    
    @property
    def RATEMETHOD(self):
        return {"Invalid" :0, "Metropolis"    :1, \
                "Kawasaki":2, "EntropyEnthalpy":3}

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

    @property
    def SIMULATION_MODE(self):
        return {"Normal":                   0x0010,
                "First Step":               0x0030,
                "Python Module":            0x0040,
                "Python Module:First Step": 0x0060,
                "Energy Only":              0x0100}

    @property
    def SIMULATION_MODE_FLAG(self):
        return {"Normal":                   0x0010,
                "First Bimolecular":        0x0020,
                "Python Module":            0x0040}

    @property
    def STOPRESULT(self):
        return {"Normal":                   0x0011,
                "Time":                     0x0012,
                "Forward":                  0x0021,
                "Time (First Step)":        0x0022,
                "Reverse":                  0x0024,
                "Error":                    0x0081,
                "NaN":                      0x0082,
                "No Moves":                 0x0084}
    
    def __setattr__(self, name, value):
        if hasattr(self, name):
            pass
        else:
            object.__setattr__(self, name, value)
    
    def __delattr__(self, *args, **kargs):
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

        self.simulation_mode = _OC.SIMULATION_MODE['Normal']
        """ The simulation mode: how we want the simulation system to
        perform the main loop.

        'Normal'          : Normal markov chain process.
        'First Step'      : Markov chain where the first move chosen is
                            always a bimolecular join step. See
                            thesis/other docs for more info,
                            especially on the statistics gathered by
                            this mode.
                            
        'Python Module'   : The simulator is compiled as a python
                            module.  This means the main loop is
                            outside the simulator, and it should
                            provide feedback to the controlling
                            process by way of this options object.
                            
        'Python Module:First Step'   : Python module, running in first
                                       step mode. (see above)

        'Energy Only'     : Compute the energy of start structure only,
                            printing the result and finishing.
                            
        Should use the values in the _OC.SIMULATION_MODE dictionary rather
        than the numbers directly, as those should be consistent with
        the ones defined in the options_python.h headers.
        """
        
        self.simulation_time = 600.0
        """ Maximum time (in seconds) allowed for each trajectory.
        
        Type         Default
        double       600.0
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
        
        self.name_dict = {}
        """ Dictionary from strand name to a list of unique strand objects
        having that name.
        
        Type         Default
        dict         {}
        
        Modified when start state is added. Used as a lookup when stop states 
        are added.
        """
        
        ####################
        #
        # BEGIN startstop
        #
        ####################
        
        # See accessors below
        self._start_state = []
        
        self.use_resting_states = False
        """ Indicates whether the start state will be determined by sampling
        from resting states or by using single structures.
        
        Type         Default
        boolean      False
        
        Automatically set to True if RestingState objects are given as the 
        start state.
        """

        # See accessors below.
        self._boltzmann_sample = None
                
        # See accessors below
        self._stop_conditions = []
        
        self.use_stop_conditions = True
        """ Indicates whether trajectories should end when stop states
        are reached.
        
        Type            Default
        boolean         True: End trajectory upon reaching stop state
        
        Defaults to ending trajectories upon reaching stop states, but can be
        manually changed to False to avoid stopping at defined stop states.
        """
        
        self.stop_count = 0
        """ The number of stop states. Equivalent to 'len(self.stop_conditions)'.
        
        Type         Default
        int          0
        
        Incremented automatically when a stop state is added. Should not be 
        modified externally.
        """
        
        self.output_time = -1.0
        """ The amount of time (in seconds) to wait between outputs of 
        trajectory information.
        
        Type         Default
        float        -1.0
        
        A value of -1.0 corresponds to not basing outputs on output_time
        (but perhaps outputting based on some other condition). A value of 0 
        means output as often as possible.
        """
        
        self.output_interval = -1
        """ The number of states between outputs of trajectory information.
        
        Type         Default
        int          -1
        
        A value of -1 corresponds to not basing outputs on output_interval
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
        
        self.interface = Interface()
        
        
        ##############################
        #
        # functions used
        #
        ##############################
        
        ('resetCompleted', 'reset_completed')
        ('resetCompleted_Python', 'reset_completed__python')
        ('setCollisionRate_Python', 'set_collision_rate__python')
        ('setCurSimTime', 'set_cur_sim_time')


    @property
    def boltzmann_sample(self):
        """ Indicates whether the start state will be determined by Boltzmann
        sampling or by using exact structures.
        
        Type         Default
        boolean      False
        
        Must be set to True when RestingState objects are given as the
        start state in order to activate the boltzmann sampling.
        """
        if self._boltzmann_sample is None:
            return False
        else:
            return self._boltzmann_sample

    @boltzmann_sample.setter
    def boltzmann_sample(self, val):
        # Set all resting states to use this boltzmann flag.
        if not isinstance( val, bool ):
            raise ValueError("the boltzmann_sample property can only be a boolean, sorry. When set, it applies globally to all resting state used in the start complexes, unless they have already been set.")
        for c,s in self._start_state:
            if s is not None:
                s.set_boltzmann( val )
        
    
    @property
    def start_state(self):
        """ Get the start state, i.e. a list of Complex objects.
        
        Type         Default
        list         []
        
        This should be used by ssystem.cc to get the (potentially sampled) 
        start state.
        """
        def process_state(x):
            cmplx, rest_state = x
            if rest_state is None:
                return cmplx
            else:
                return rest_state.sample()
            
        return [process_state(s) for s in self._start_state]
    
    @start_state.setter
    def start_state(self, *args):
        """ Set the start state, i.e. a list of Complex or RestingState objects.
        
        Type         Default
        list         []
        
        The start state should be set (e.g. by the parser) so trajectories know 
        how to start.
        """
        # Error checking first
        if self._start_state !=  []:
            raise Exception("Start state should only be set once.")
        if len(args) == 0 or len(args[0]) == 0;
            raise ValueError("No start state given.")
        
        # deduce our input from the type of args[0].
        # Copy the input list because it's easy to do and it's safer
        
        if isinstance(args[0],Complex) or isinstance(args[0],RestingState):
            #args is a list of complexes or resting states.
            vals = copy.deepcopy(args) 
        elif len(args) == 1 and hasattr(args[0],"__iter__"):
            vals = copy.deepcopy(args[0])
        else:
            raise ValueError("Could not comprehend the start state you gave me.")

        # vals is now an iterable over our starting configuration, be
        # it complexes or resting states.
        
        for i in vals:
            if not isinstance(i,Complex) and not isinstance(i,RestingState):
                raise ValueError("Start states must be Complexes or RestingStates. Received something of type {0}.".format( type(i)) )
        
            self._add_start_complex( i )

    def _add_start_complex( self, item ):
        l = len(item)
        item.make_unique( self.name_dict, range(self.unique_id,self.unique_id+l) )

        self.unique_id += l
        
        item.set_boltzmann( self.boltzmann_sample )
        if isinstance(item, Complex):
            self._start_state.append( (item,None))
        else:
            self._start_state.append( (item[0], item) )

    @property
    def initial_seed_flag(self):
        return self.initial_seed != None

    
    @property
    def stop_conditions(self):
        """ The stop states, i.e. a list of StopCondition objects.
        
        Type         Default
        list         []
        
        Stop states should be added to this list (e.g. by the parser) so
        trajectories know when to end.
        """
        return self._stop_conditions
    
    
    @stop_conditions.setter
    def stop_conditions(self, stop_list):
        """ The stop states, i.e. a list of StopCondition objects.
        
        Type         Default
        list         []
        
        Stop states should be added to this list (e.g. by the parser) so
        trajectories know when to end.
        """
        # Error checking
        if self._stop_conditions != []:
            raise Exception("Stop conditions should be set only once.")
        if self._start_state == []:
            raise Exception("Start state must be set before stop conditions.")
        
        # Type checking
        for item in stop_list:
            if not isinstance(item, StopCondition):
                raise TypeError("All items must be 'StopCondition', not '{0}'.".format(type(item)))
        
        # Copy the input list because it's easy to do and it's safer
        stop_list = copy.deepcopy(stop_list)
        
        # Assign the unique ids
        for sc in stop_list:
            counts = {}
            for cmplx, st, cn in sc.complex_items:
                for i, s in enumerate(cmplx.strand_list):
                    if s.name not in self.name_dict:
                        raise ValueError("Stop state contains a strand not present in start state.")
                    try:
                        cmplx.strand_list[i] = self.name_dict[s.name][counts[s.name]]
                        counts[s.name] += 1
                    except KeyError:
                        cmplx.strand_list[i] = self.name_dict[s.name][0]
                        counts[s.name] = 1
        
        # Set the internal data members
        self.stop_count = len(stop_list)
        self._stop_conditions = stop_list
    
    @property
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

        return None
    
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

        Yes, these ranges are quite generous.
        """
        if 273.0 < val < 373.0:
            self._temperature_kelvin = val
            self._temperature_celsius = val - _OC.ZERO_C_IN_K
        elif 0.0 < val < 100.0:
            self._temperature_celsius = val
            self._temperature_kelvin = val + _OC.ZERO_C_IN_K
            self.errorlog.append("Warning: Temperature was set at the value [{0}]. We expected a value in Kelvin, or with appropriate units.\n         Temperature was automatically converted to [{1}] degrees Kelvin.\n".format(val, self._temperature_kelvin))
    
    
    def make_unique(self, strand):
        """Returns a new Strand object with a unique identifier replacing the 
        old id. Also adds the new strand to self.name_dict[strand.name].
        """
        new_strand = Strand(self.unique_id, strand.name, strand.sequence, strand.domain_list)
        self.unique_id += 1
        
        try:
            self.name_dict[strand.name].append(new_strand)
        except KeyError:
            self.name_dict[strand.name] = [new_strand]

        return new_strand
    
    @property
    def interface_reset_completion_flag(self):
        tmp = self.interface.trajectory_completion_flag
        self.interface.trajectory_completion_flag = False
        self.interface.collision_rate = -1.0
        return tmp

    @property
    def interface_collision_rate(self):
        return self.interface.collision_rate

    @interface_collision_rate.setter
    def interface_collision_rate(self, rate ):
        self.interface.collision_rate = rate

    @property
    def interface_current_time(self):
        return self.interface_current_time

    @interface_current_time.setter
    def interface_current_time(self, time):
        self.interface.current_time = time
    
    @property
    def interface_halt_flag(self):
        return self.interface.trajectory_halt_flag

    @property
    def interface_suspend_flag(self):
        return self.interface.trajectory_suspend_flag

    @property
    def increment_trajectory_count(self):
        self.interface.increment_trajectory_count()

    @property
    def interface_current_seed(self):
        return self.interface.current_seed

    @interface_current_seed.setter
    def interface_current_seed(self, val):
        self.interface.current_seed = val

    @property
    def interface_current_tag(self):
        return self.interface.current_tag

    @interface_current_tag.setter
    def interface_current_tag(self, val):
        self.interface.current_tag = val

    @property
    def interface_current_completion_type(self):
        return self.interface.current_completion_type


    @interface_current_completion_type.setter
    def interface_current_completion_type(self, val):
        self.interface.current_completion_type = val
