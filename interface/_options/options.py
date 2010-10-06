################################################################################
#                                                                              #
# Python implementation of the options object.                                 #
# Copyright 2010 Caltech                                                       #
# Written by:  Joseph Schaeffer.                                               #
# Some stuff written by:  Chris Berlind                                        #
#                                                                              #
# (others add your names as you modify files)                                  #
#                                                                              #
################################################################################

from interface import Interface
from ..objects import Strand, Complex, RestingState, StopCondition
from constants import _OptionsConstants

_OC = _OptionsConstants()

import copy

FILLIN = None

class Options( object ):
    """ The main wrapper for controlling a Multistrand simulation. Has information about the energy model, simulation options, and an interface for returning results. """

    defaults = {"num_simulations":1,
                "log_ml":False,
                "substrate_type":2,
                "join_concentration":1.0,
                "simulation_mode":16,
                "parameter_type":1,
                "simulation_time":600.0,
                "temperature":310.15}
    
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

        
        self._unimolecular_scaling = 1.6e6
        """ Rate scaling factor for unimolecular reactions.

        Type         Default
        double       1.6e6:
                     Unitless. Details on default in thesis."""
        
        self._bimolecular_scaling = 0.5e6
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

        self.__init_keyword_args(self, *args, **kargs )
        

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
                s.boltzmann_sample = val
            else:
                c.boltzmann_sample = val 
        
    @property
    def unimolecular_scaling( self ):
        """ Rate scaling factor for unimolecular reactions.

        Type         Default
        double       1.6e6:
            Unitless. Details on default in thesis."""
        return self._unimolecular_scaling

    @unimolecular_scaling.setter
    def unimolecular_scaling( self, val ):
        self._unimolecular_scaling = float( val )
        
    @property
    def bimolecular_scaling( self ):
        """ Rate scaling factor for bimolecular reactions.
        double       0.5e6:
                     Unitless. Details on default in thesis."""
        return self._bimolecular_scaling

    @bimolecular_scaling.setter
    def bimolecular_scaling( self, val ):
        self._bimolecular_scaling = float( val )
    
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
        if self._start_state !=  []:
            raise Exception("Start state should only be set once.")
        if len(args) == 0 or len(args[0]) == 0:
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
        import warnings
        warnings.warn("Unique ID assignment likely incorrect here.")
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
        self._use_stop_conditions = True


    @property
    def use_stop_conditions( self ):
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
    def use_stop_conditions( self, val ):
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
        else:
            self._temperature_kelvin = val
            self._temperature_celsius = val - _OC.ZERO_C_IN_K
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
    def add_result_status_line(self,val):
        """ Takes a 4-tuple as the only value type, it should be:
            (random number seed, stop result flag, completion time, stop result tag)

            Prints thingies. Also sets useful values."""
        if not isinstance(val, tuple) or len(val) != 4:
            raise ValueError("Print status line needs a 4-tuple of values.")
        self.interface.add_result( val, res_type = 'status_line' )
        
        # seed,com_type, time, tag = val
        # self.interface_current_seed = seed  #uses property to get the right sub-object.
        # self.interface_current_completion_type = com_type
        # self.interface_current_time = time
        # self.interface_current_tag = tag

    @property
    def add_result_status_line_firststep(self, val):
        return None

    @add_result_status_line_firststep.setter
    def add_result_status_line_firststep( self, val ):
        pass

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

    @property
    def increment_trajectory_count( self ):
        self.interface.increment_trajectory_count()

    def __repr__( self ):
        items_to_save = {}
        for k in Options.defaults.keys():
            if self.__getattribute__(k) != Options.defaults[k]:
                items_to_save[k] = self.__getattribute__(k)
        def prep(k,v):
            return "{key}={value},\n".format(key=k,value=v)
        if len(items_to_save) == 0:
            res = "Options()\n"
        else:
            res = "Options( " + "".join( [prep(*items_to_save.items()[0])] +
                                     [pad + prep(*item) for item,pad in zip( items_to_save.items()[1:], ["         "] * (len(items_to_save) - 1))])
            res = res[:-2]
            res = res + ")\n"
        return res

    def __init_keyword_args( self, *args, **kargs ):
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
            'dangles': lambda x: self.__setattr__('dangles',_OC.DANGLES[x]),
            'parameter_type': lambda x: self.__setattr__('parameter_type', _OC.ENERGYMODEL_TYPE[x]),
            'substrate_type': lambda x: self.__setattr__('substrate_type', _OC.SUBSTRATE_TYPE[x]),
            'biscale': lambda x: self.__setattr__('bimolecular_scaling', x),
            'uniscale': lambda x: self.__setattr__('unimolecular_scaling', x),
            'num_sims': lambda x: self.__setattr__('num_simulations', x),
            'sim_time': lambda x: self.__setattr__('simulation_time',x),
            'concentration': lambda x: self.__setattr__('join_concentration',x)
            }
        
        for k in kargs.keys():
            if k in arg_lookup_table:
                arg_lookup_table[k]( kargs[k] )
            else:
                self.__setattr__(k, kargs[k])
