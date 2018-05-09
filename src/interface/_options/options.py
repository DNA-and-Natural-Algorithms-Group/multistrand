# Options object.                                 
# Copyright 2010-2017 Caltech                                                  
# Joseph Schaeffer                                                             
# Chris Berlind                                                                
# Frits Dannenberg                                                             

from interface import Interface
from ..objects import Strand, Complex, RestingState, StopCondition
# from ..utils import JSKawasaki25, JSKawasaki37, JSMetropolis25, JSMetropolis37, JSDefault     # this is for 2.0 support only

import copy


class Options(object):
    """ The main wrapper for controlling a Multistrand simulation. Has an interface for returning results. """
       
    '''Constants:'''
    # FD: some could be enums, using constants instead.
    
    ZERO_C_IN_K = 273.15
    
    STR_ALL = "ALL"
    STR_FAILURE = "FAILURE"
    STR_SUCCESS = "SUCCESS"
    STR_ALT_SUCCESS = "ALT_SUCCESS"
    
    """
         protected results may occupy the [x.tag in a x in options.interface.results]
    """
    STR_TIMEOUT = "timeout"
    STR_NOINITIAL = "noinitial"
    STR_NAN = "nan"
    STR_ERROR = "error"
        
    # rate_method
    metropolis = 1
    kawasaki = 2
    RateMethodToString = [ "None", "Metropolis", "Kawasaki"]

    # Nupack dangle options
    none = 0
    some = 1
    all = 2
    dangleToString = [ "None", "Some", "All"]

    # Parameter type. Vienna is depreciated.
    viennaModel = 0 
    nupackModel = 1 
    parameterTypeToString = [ "Vienna", "Nupack" ]

    # Substrate type.
    substrateRNA = 1
    substrateDNA = 2
    substrateToString = [ "Invalid", "RNA", "DNA"]    
    
    # Simulation mode
    firstPassageTime = 16  # 0x0010
    firstStep = 48  # 0x0030
    transition = 256  # 0x0100
    trajectory = 128  # 0x0080
    
    # translation
    simulationMode = {  "Normal"    :               firstPassageTime,
                        "First Step":               firstStep,
                        "Transition":               transition,
                        "Trajectory":               trajectory,
                        "First Passage Time":       firstPassageTime}


    """
        FD, May 8th, 2018:
        
        StopCondition and Macrostate definitions:
            Exact - match a secondary structure exactly (i.e. any system state that has a complex with this exact structure)
            Bound - match any system state in which the given strand is bound to another strand
            Dissoc - match any system state in which there exists a complex with exactly the given strands, in that order
            Loose - match a secondary structure, allowing a certain number of disagreements, allowing domain state to be unspecified (* notation)
            Count - match a secondary structure, allowing a certain number of disagreements
            
            Compuational expense, low to high :           dissoc -- exact -- count -- loose     
    """

    exactMacrostate = 0  # 
    boundMacrostate = 1  # 
    dissocMacrostate = 2  # 
    looseMacrostate = 3  # 
    countMacrostate = 4  # 


    
    cotranscriptional_rate_default = 0.001  # 1 nt added every 1 ms
    
    activestatespace = False;
    
    def __init__(self, *args, **kargs):
        """
        Initialization of an Options object:

        Keyword Arguments:
        dangles -- Specifies the dangle terms used in the energy model.
                   Can be the strings 'None', 'Some', or 'All'.
        start_state  [type=list]     -- A list of Complex or RestingStates to
                                        use as the initial state of the system.
        simulation_time [type=float] -- Cap on the maximum simulation time.
        num_simulations [type=int]   -- Number of trajectories to run
        biscale         [type=float] -- Bimolecular scaling constant
        uniscale        [type=float] -- Unimolecular scaling constant
        parameter_type               -- Which set of energy parameters is
                                        used. Available options: 'Nupack',
                                        'Vienna'
        substrate_type               -- Whether we want 'DNA' or 'RNA' energy
                                        parameters.
        rate_method                  -- Whether we want 'Kawasaki' or 'Metropolis'
                                        rate method for unimolecular steps.
        useArrRates     [type=bool]  -- if TRUE, use Arrhenius rate model. If using, please set lnAEnd, lnALoop, lnAStack, lnAStackStack, lnALoopEnd, lnAStackEnd, lnAStackLoop
EEnd, ELoop, EStack, EStackStack, ELoopEnd, EStackEnd, EStackLoop (double value).
        """

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
        self.full_trajectory = []
        self.full_trajectory_times = []
        self.full_trajectory_arrType = []
        self.trajectory_complexes = []
        self.trajectory_state_count = 0
        self._current_end_state = []
        self._current_transition_list = []
        self.special_count = 0
        self.trajectory_current_time = 0.0
        self.current_graph = None

        self.verbosity = 1
        """ Indicates how much output will be generated for each trajectory run.
        Value = 0:  No end state reported, no warnings for timeout and nonitial steps
        Value = 1:  No end states reports, warnings active   
        Value = 2:  warnings and end states reports to stdout
        """
        
        self.print_initial_first_step = False
        """
        If True, this value will print the initial state in First Step Mode to the trajectory with a timestamp of -1.0
        """

        self.cotranscriptional = False
        """
        If True, enables the cotranscriptional simulation mode. The mode works only when a single strand is supplied.
        Starting with the initial 8 nucleotides, the simulation adds a nucleotide on the 3' end every 1 millisecond. 
        """
        
        self.cotranscriptional_rate = self.cotranscriptional_rate_default
        """
        By default, the cotranscriptional mode adds one nucleotide every 1 millisecond.
        """
        
        #############################################
        #                                           #
        # Data Members: Energy Model                #
        #                                           #
        #############################################
        
        self.gt_enable = False
        """ Allow GT base pairs? If not, penalize by 10000 kcal/mol.
            False (0) : Do not allow GT base pairs.
        """

        self.log_ml = False
        """ Use logarithm to compute multiloop energy contributions?
        int          False (0): Do not use the logarithm.

        If True, uses log to compute one component of the multiloop energy, for loops of length > 6. Otherwise, uses the usual
        linear approximation. Using the log formula is slightly more expensive as it makes computation time for a multiloop scale
        with the number of adjoining helices.
        """
        
        self.join_concentration = 1.0
        """ concentration for V calcs
        Units are in M (molar), and represent the concentration of a single unique strand in the system. The volume simulated is
        then chosen using this parameter.
        """

        # ##
        # ## See the temperature property way below (after __init__)
        # ## for more info on accessors for these data members.
        # ##
        self._temperature_celsius = 37.0
        self._temperature_kelvin = 310.15

        self.rate_scaling = None
        """FD: This is a legacy option that sets unimolecular and bimolecular scaling automatically if set"""

        self.unimolecular_scaling = -1.0 
        """ Rate scaling factor for unimolecular reactions."""
        
        self.bimolecular_scaling = -1.0 
        """ Rate scaling factor for bimolecular reactions."""

        self.rate_method = self.kawasaki
        """ Choice of methods for determining forward/reverse rates. """

        self.dangles = self.some
        """ Dangles options for the energy model.
        
        None [0]: Do not include any dangles terms in the energy model.
        Some [1]: Some dangles terms.  (Nupack Default)
        All  [2]: Include all dangles terms, including odd overlapping ones.
        """

        self.parameter_type = self.nupackModel
        """ Which type of energy model parameter file to use.

        Vienna [0]: No longer well tested. Recommend not using.
        Nupack [1]: Includes some multi-complex parameters, otherwise
                    nearly the same as mfold style files.
        """

        self.substrate_type = self.substrateDNA
        """ What substrate's parameter files to use. 

        Invalid [0]: Indicates we should not auto-search for a param file.
        RNA     [1]: RNA parameters are available for Vienna and Nupack, see also
                     the comment in parameter_file_version.
        DNA     [2]: DNA parameters are publicly available via the Nupack distribution,
                     and possibly by the Vienna group upon request.
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
        
        self.simulation_mode = self.firstPassageTime
        """ The simulation mode: how we want the simulation system to
        perform the main loop.
        """
        
        self.simulation_time = 600.0
        """ Maximum time (in seconds) allowed for each trajectory.
        
        Type         Default
        double       600.0
        """
        
        self.num_simulations = 1
        """ Total number of trajectories to run. 
        """
        
        self.initial_seed = None
        """ Initial random number seed to use.
        If None when simulation starts, a random seed will be chosen
        """
        
        self.name_dict = {}
        """ Dictionary from strand name to a list of unique strand objects
        having that name.
        
        Type         Default
        dict         {}
        
        Modified when start state is added. Used as a lookup when stop states 
        are added.
        """
        
        self.useArrRates = False;
        
        self.lnAEnd = -0.1;
        self.lnALoop = -0.1;
        self.lnAStack = -0.1;
        self.lnAStackStack = -0.1;
        self.lnALoopEnd = -0.1;
        self.lnAStackEnd = -0.1;
        self.lnAStackLoop = -0.1;        
        self.EEnd = -0.1;
        self.ELoop = -0.1;
        self.EStack = -0.1;
        self.EStackStack = -0.1;
        self.ELoopEnd = -0.1;
        self.EStackEnd = -0.1;
        self.EStackLoop = -0.1; 
        
        self.dSA = -0.0;
        self.dHA = -0.0;

        # Buffer conditions
        self.sodium = 1.0;
        self.magnesium = 0.0;
        
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
        self._use_stop_conditions = False
        
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
        # End of __init__: call the keyword hook fn. 
        #
        ##############################

        self.__init_keyword_args(self, *args, **kargs)

    def legacyRates(self):
                    
        warningmsg = "Warning! rate_scaling is set, enabling support for legacy code. Now setting rate defaults for "
                
        if (self.temperature == 298.15) & (self.rate_method == self.kawasaki) :
            warningmsg += "Kawasaki 25 C"
            self.JSKawasaki25()
        elif (self.temperature == 310.15) & (self.rate_method == self.kawasaki) :
            warningmsg += "Kawasaki 37 C"
            self.JSKawasaki37()
        elif (self.temperature == 298.15) & (self.rate_method == self.metropolis) :
            warningmsg += "Metropolis 25 C"
            self.JSMetropolis25()
        elif (self.temperature == 310.15) & (self.rate_method == self.metropolis) :
            warningmsg += "Metropolis 37 C"
            self.JSMetropolis37()
        else:
            warningmsg += "JS-Default"
            self.JSDefault()
            
        print warningmsg 
        self.rate_scaling = None       
        
    # FD, May 5th 2017
    # Supplying rate options for Metropolis and Kawasaki methods,
    # all using the dangles = some option. Also:  one general default,
    # and one setting for Metropolis rates derived for DNA23. 
    
    def JSDefault(self): 
        """ Default rates from Joseph Schaeffer's thesis  """
        
        self.unimolecular_scaling = 1.50e+08;
        self.bimolecular_scaling = 1.38e+06;
    
    def JSMetropolis25(self): 
        """ Default rates for Metropolis at 25 degree Celcius, from Joseph Schaeffer's thesis
        
        """
        
        self.unimolecular_scaling = 4.4e8;
        self.bimolecular_scaling = 1.26e6;
    
    def JSKawasaki25(self): 
        """ Default rates for Kawasaki at 25 degree Celcius, from Joseph Schaeffer's thesis
        
        """
        
        self.unimolecular_scaling = 6.1e7;
        self.bimolecular_scaling = 1.29e6;
    
    def JSKawasaki37(self):
        """ Default rates for Kawasaki at 37 degree Celcius, from Joseph Schaeffer's thesis
        """
        
        self.unimolecular_scaling = 1.5e8;
        self.bimolecular_scaling = 1.38e6;
    
    def JSMetropolis37(self): 
        """ Default rates for Metropolis at 37 degree Celcius, from Joseph Schaeffer's thesis
        """
        
        self.unimolecular_scaling = 7.3e8;
        self.bimolecular_scaling = 1.40e6;
    
    def DNA23Metropolis(self):
        """ A default rate for Metropolis at 25 degree Celcius, from the DNA23 conference
        """

        # FD march 2018: Setting these to be the 55th walker from the dna23 inference
        # previous values:
#         self.unimolecular_scaling = 5.0e6;
#         self.bimolecular_scaling = 1.4e6;
        
        self.unimolecular_scaling = 2.41686715e+06;
        self.bimolecular_scaling = 8.01171383e+05 ;
     
    def uniformRates(self):
        """ uniform rates without a source
        """
            
        self.unimolecular_scaling = 1.0e6;
        self.bimolecular_scaling = 1.0e6;

# FD: We are implementing an observer pattern. 
# FD: After temperature, substrate (RNA/DNA) or danlges is updated, we attempt to update boltzmann samples.
    def updateBoltzmannSamples(self):

        if len(self._start_state) > 0:
            for c, s in self._start_state:
                c.set_boltzmann_parameters(self.dangleToString[self.dangles], self.substrateToString[self.substrate_type], self._temperature_celsius)
    
#     FD: We are using shadow variables for biscale, uniscale only because we want to 
#     flag a warning when rate_scaling is used
    @property
    def bimolecular_scaling(self):

        if self.rate_scaling != None :
            self.legacyRates()
         
        return self._bimolecular_scaling
 
    @bimolecular_scaling.setter
    def bimolecular_scaling(self, val):
                  
        self._bimolecular_scaling = float(val)
 
    @property
    def unimolecular_scaling(self):
         
        if self.rate_scaling != None :
            self.legacyRates()
         
        return self._unimolecular_scaling  
    
    @unimolecular_scaling.setter
    def unimolecular_scaling(self, val):
                  
        self._unimolecular_scaling = float(val)
        
# FD: Shadow variables for danlges only because we need to observe changes (and update boltzmann samples accordingly)
    @property
    def dangles(self):
        return self._dangles
    
    @dangles.setter
    def dangles(self, val):
        self._dangles = int(val)
        self.updateBoltzmannSamples
        
# FD: shadow parameter so that boltzmann samples can be updated when this parameter is set
# FD: In a better control flow, complexes themselves might fetch the right constants just before evaluating their boltzmann samples 
    @property
    def substrate_type(self):
        return self._substrate_type

    @substrate_type.setter
    def substrate_type(self, value):
        self._substrate_type = int(value)
        self.updateBoltzmannSamples

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
        if not isinstance(val, bool):
            raise ValueError("the boltzmann_sample property can only be a boolean, sorry. When set, it applies globally to all resting state used in the start complexes, unless they have already been set.")
        for c, s in self._start_state:
            if s is not None:
                s.boltzmann_sample = val
            else:
                c.boltzmann_sample = val 
    
    @property
    def start_state(self):
        """ Get the start state, i.e. a list of Complex objects.
        
        Type         Default
        list         []
        
        This should be used by ssystem.cc to get the (potentially sampled) 
        start state.
        """

        #        import pdb
        #        pdb.set_trace()
        def process_state(x):
            cmplx, rest_state = x
            if rest_state is None:
                return cmplx
            else:
                return rest_state.get_starting_complex()
            
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
        if self._start_state != []:
            raise Exception("Start state should only be set once.")
        if len(args) == 0 or len(args[0]) == 0:
            raise ValueError("No start state given.")
        
        # deduce our input from the type of args[0].
        # Copy the input list because it's easy to do and it's safer
        
        if isinstance(args[0], Complex) or isinstance(args[0], RestingState):
            # args is a list of complexes or resting states.
            vals = copy.deepcopy(args) 
        elif len(args) == 1 and hasattr(args[0], "__iter__"):
            vals = copy.deepcopy(args[0])
        else:
            raise ValueError("Could not comprehend the start state you gave me.")

        # vals is now an iterable over our starting configuration, be
        # it complexes or resting states.
        
        for i in vals:
            if not isinstance(i, Complex) and not isinstance(i, RestingState):
                raise ValueError("Start states must be Complexes or RestingStates. Received something of type {0}.".format(type(i)))
        
            self._add_start_complex(i)

    def _add_start_complex(self, item):
        if isinstance(item, Complex):
            self._start_state.append((item, None))
            item.set_boltzmann_parameters(self.dangleToString[self.dangles], self.substrateToString[self.substrate_type], self._temperature_celsius)
        else:
            self._start_state.append((item[0], item))
            # JS 8/30/2014 Not entirely sure about this one. If I remember right, resting states may have more than one
            # component, so setting only the first to have proper Boltzmann sampling might be an error.
            # 
            # If so, it should probably be 'for i in item: i.set_boltzmann_parameters'... {etc}
            item[0].set_boltzmann_parameters(self.dangleToString[self.dangles], self.substrateToString[self.substrate_type], self._temperature_celsius)

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

        # Type checking
        for item in stop_list:
            if not isinstance(item, StopCondition):
                raise TypeError("All items must be 'StopCondition', not '{0}'.".format(type(item)))
        
        # Copy the input list because it's easy to do and it's safer
        stop_list = copy.deepcopy(stop_list)
        
        # Set the internal data members
        self.stop_count = len(stop_list)
        self._stop_conditions = stop_list
        self._use_stop_conditions = True

    @property
    def use_stop_conditions(self):
        """ Indicates whether trajectories should end when stop states
        are reached.
        
        Type            Default
        boolean         False: End trajectory upon reaching max time only.
        
        Defaults to ending trajectories only on the max simulation
        time, but setting any stop conditions automatically changes
        this to True, and it will stop whenever it reaches a stop
        condition, or at max time [whichever comes first].
        """
        return self._use_stop_conditions

    @use_stop_conditions.setter
    def use_stop_conditions(self, val):
        if val == True and len(self._stop_conditions) == 0:
            import warnings
            warnings.warn("Options.use_stop_conditions was set to True, but no stop conditions have been defined!")

        self._use_stop_conditions = val
    
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
        """
        Temperature, in degrees Kelvin.
        
        Arguments:
        temperature [type=float,default=310.15] -- Standard units of Kelvin.
               Default value corresponds to 37(C).
        
        This is used mostly in scaling energy model terms, and those
        scaling factors always use Kelvin. Note that when set, the
        Options object will try to make sure it's actually a 'sane'
        value, as follows:
        
        Temperatures in the range [0,100] are assumed to be Celsius,
        and are converted to Kelvin.
        
        Temperatures in the range [273,373] are assumed to be in
        Kelvin.
        
        Any temperature outside these ranges is set as the temperature
        in Kelvin, and a warning is raised. If any conversion takes
        place, a message is added to the Options object's errorlog.
        """
        return self._temperature_kelvin

    @temperature.setter
    def temperature(self, val):
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
            self._temperature_celsius = val - self.ZERO_C_IN_K
            self.updateBoltzmannSamples

        elif 0.0 < val < 100.0:
            self._temperature_celsius = val
            self._temperature_kelvin = val + self.ZERO_C_IN_K
            self.updateBoltzmannSamples
            self.errorlog.append("Warning: Temperature was set at the value [{0}]. We expected a value in Kelvin, or with appropriate units.\n         Temperature was automatically converted to [{1}] degrees Kelvin.\n".format(val, self._temperature_kelvin))

        else:
            self._temperature_kelvin = val
            self._temperature_celsius = val - self.ZERO_C_IN_K
            self.updateBoltzmannSamples
            self.errorlog.append("Warning: Temperature was set at the value [{0}]. This is outside the normal range of temperatures we expect, so it was assumed to be in Kelvin.\n".format(val))
            raise Warning("Temperature did not fall in the usual expected ranges. Temperatures should be in units Kelvin, though the range [0,100] is assumed to mean units of Celsius.")
    
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
    def add_result_status_line(self):
        return None

    @add_result_status_line.setter
    def add_result_status_line(self, val):
        """ Takes a 4-tuple as the only value type, it should be:
            (random number seed, stop result flag, completion time, stop result tag)

            Prints thingies. Also sets useful values."""
        if not isinstance(val, tuple) or len(val) != 4:
            raise ValueError("Print status line needs a 4-tuple of values.")
        self.interface.add_result(val, res_type='status_line')
        if len(self._current_end_state) > 0:
            self.interface.end_states.append(self._current_end_state)
            self._current_end_state = []
    
    @property
    def add_result_status_line_firststep(self, val):
        return None

    @add_result_status_line_firststep.setter
    def add_result_status_line_firststep(self, val):
        """ Takes a 5-tuple as the only value type, it should be:
            (random number seed, stop result flag, completion time, collision rate, stop result tag)

            Prints thingies. Also sets useful values."""
        if not isinstance(val, tuple) or len(val) != 5:
            raise ValueError("Print status line needs a 5-tuple of values.")
        self.interface.add_result(val, res_type='firststep')
        if len(self._current_end_state) > 0:
            self.interface.end_states.append(self._current_end_state)
            self._current_end_state = []
            
    @property
    def add_complex_state_line(self):
        return None

    @add_complex_state_line.setter
    def add_complex_state_line(self, val):
        """ Takes a 6-tuple as only value, it should be:
            (random number seed, unique complex id, strand names, sequence, structure, energy )
            Adds this data to the interface's results object."""

        self._current_end_state.append(val)
        if self.verbosity > 1:
            print("{0[0]}: [{0[1]}] '{0[2]}': {0[5]} \n{0[3]}\n{0[4]}\n".format(val))

    @property
    def add_transition_info(self):
        return None

    @add_transition_info.setter
    def add_transition_info(self, val):
        """ Takes a 2-tuple, 1st value is current time, 2nd value is a
            list of boolean values indicating which stop conditions we
            currently meet."""
        # print( "Time: {0[0]} Membership: {0[1]}".format( val ))
        self._current_transition_list.append(val)

    @property
    def add_trajectory_complex(self):
        return None

    @add_trajectory_complex.setter
    def add_trajectory_complex(self, val):
        """ Takes a 6-tuple as only value, it should be:
            (random number seed, unique complex id, strand names, sequence, structure, energy )
            Adds this data to the interface's results object."""

        self.trajectory_complexes.append(val)
        # print "Number of complexes in current traj output:" + str(len(self.trajectory_complexes))

    @property
    def add_trajectory_current_time(self):
        return None

    @add_trajectory_current_time.setter
    def add_trajectory_current_time(self, val):
        self.trajectory_current_time = val
        self.trajectory_state_count += 1
        self.full_trajectory.append(self.trajectory_complexes)
        self.full_trajectory_times.append(self.trajectory_current_time)
        self.trajectory_complexes = []
        
    @property
    def add_trajectory_arrType(self):
        return None

    @add_trajectory_arrType.setter
    def add_trajectory_arrType(self, val):
        self.full_trajectory_arrType.append(val)
        
    @property
    def interface_current_seed(self):
        """ This is the current random number seed for the trajectory currently being
        simulated by multistrand.

        This property and its setter exist mainly for use by the
        multistrand module to provide some feedback on what seed it's
        working on, in case that trajectory crashes before completion,
        and other similar cases."""
        return self.interface.current_seed

    @interface_current_seed.setter
    def interface_current_seed(self, val):
        
        self.interface.current_seed = val
        
        def get_structure(s):
            if s.boltzmann_sample:
                return s._last_boltzmann_structure
            else:
                return s._fixed_structure
        
        def process_state(x):
            cmplx, rest_state = x
            if rest_state is None:
                return get_structure(cmplx)
            else:
                return get_structure(rest_state.get_starting_complex())
        
        structures = [process_state(s) for s in self._start_state]
        
        self.interface.start_structures[val] = structures

    @property
    def increment_trajectory_count(self):
        self.interface.increment_trajectory_count()
        if len(self._current_transition_list) > 0:
            self.interface.transition_lists.append(self._current_transition_list)
            self._current_transition_list = []

    def __init_keyword_args(self, *args, **kargs):
        """ Helper subfunction. """
        
        """ Create an options object [with default (presumably useful) values]

        Now with new and improved argument lists!

        Any default Options attribute can be set just by passing a
        keyword argument with the desired value, e.g.
        Options(simulation_time=100.0)

        Listed below are some shortcuts, for attributes that have a
        range of options, or shortened names for long attributes.
        

        Keyword Argument  |  Options
        dangles           |  'None', 'Some', 'All'
        parameter_type    |  'Nupack', 'Vienna'
        substrate_type    |  'DNA','RNA'
        
        sim_time          | [simulation_time] Max time to simulate
        num_sims          | [num_simulations] Number of trajectories to run
        biscale           | [bimolecular_scaling] Bimolecular scaling constant
        uniscale          | [unimolecular_scaling] Unimolecular scaling constant

        start_state       |  List of Complex or RestingState

        ...
        More to come!"""
        arg_lookup_table = {
            'biscale': lambda x: self.__setattr__('bimolecular_scaling', x),
            'uniscale': lambda x: self.__setattr__('unimolecular_scaling', x),
            'num_sims': lambda x: self.__setattr__('num_simulations', x),
            'sim_time': lambda x: self.__setattr__('simulation_time', x),
            'concentration': lambda x: self.__setattr__('join_concentration', x)
            }
        
        # FD: Start throwing errors if not in the right format
        # FD: This does not prevent the user to set them to ints after options 
        # FD: initialization (could use overloading via @property to prevent this).
        for key, value in kargs.items():
            
            if key == "simulation_time":
                if not isinstance(value, (float)):
                    raise Warning("Please provide simulation_time as float")
                
            if key == "bimolecular_scaling":
                if not isinstance(value, (float)):
                    raise Warning("Please provide bimolecular_scaling as float")
                
            if key == "unimolecular_scaling":
                if not isinstance(value, (float)):
                    raise Warning("Please provide unimolecular_scaling as float")
        
        for k in kargs.keys():
            
            if k in arg_lookup_table:
                arg_lookup_table[k](kargs[k])          
                
            # FD: Do some additional parsing for legacy support            
            # FD: This code simply translates the string calls to the numerical constants 
            elif k == 'rate_method':
                if isinstance(kargs[k], basestring):
                    self.rate_method = self.RateMethodToString.index(kargs[k])
                    
            elif k == 'dangles':
                if isinstance(kargs[k], basestring):
                    self.dangles = self.dangleToString.index(kargs[k])

            elif k == 'parameter_type':
                if isinstance(kargs[k], basestring):
                    self.parameter_type = self.parameterTypeToString.index(kargs[k])

            elif k == 'substrate_type':
                if isinstance(kargs[k], basestring):
                    self.substrate_type = self.substrateToString.index(kargs[k])
                    
            elif (k == 'simulation_mode') & (isinstance(kargs[k], basestring)):
                    self.simulation_mode = self.simulationMode[kargs[k]]

            else:
                self.__setattr__(k, kargs[k])
