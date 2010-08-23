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
from python_objects import Strand, Complex, RestingState
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
        
        self.boltzmann_sample = True
        """ Indicates whether the start state will be determined by Boltzmann
        sampling or by using exact structures.
        
        Type         Default
        boolean      True
        
        Automatically set to True if RestingState objects are given as the 
        start state.
        """
        
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
        
        ('finalizeInput', 'finalize_input')
        ('printStatusLine', 'print_status_line')
        ('printStatusLine_Final_First_Bimolecular', 'print_status_line__final__first__bimolecular')
        ('printStatusLine_First_Bimolecular', 'print_status_line__first__bimolecular')
        ('printStatusLine_Warning', 'print_status_line__warning')
        ('printTrajLine', 'print_traj_line')
        ('resetCompleted', 'reset_completed')
        ('resetCompleted_Python', 'reset_completed__python')
        ('setCollisionRate_Python', 'set_collision_rate__python')
        ('setCurSimTime', 'set_cur_sim_time')
    
    
    @property
    def start_state(self):
        """ Get the start state, i.e. a list of Complex objects.
        
        Type         Default
        list         []
        
        This should be used by ssystem.cc to get the (potentially sampled) 
        start state.
        """
        if self.use_resting_states:
            complex_list = [rs[0] for rs in self._start_state]
        else:
            complex_list = self._start_state
        
        if self.boltzmann_sample:
            complex_list = [boltzmann_sample(c) for c in complex_list]
        
        return complex_list
    
    
    @start_state.setter
    def start_state(self, start_list):
        """ Set the start state, i.e. a list of Complex or RestingState objects.
        
        Type         Default
        list         []
        
        The start state should be set (e.g. by the parser) so trajectories know 
        how to start.
        """
        # Error checking first
        if self._start_state !=  []:
            raise Exception("Start state should only be set once.")
        if start_list == []:
            raise ValueError("No start state given.")
        
        # Copy the input list because it's easy to do and it's safer
        start_list = copy.deepcopy(start_list)
        
        # Make sure all types match
        t = type(start_list[0])
        if any([type(item) is not t for item in start_list]):
            raise TypeError("All items in list must be the same type.")
        elif t is not Complex and t is not RestingState:
            raise TypeError("List items should be 'Complex' or 'RestingState', not '%s'." % t)
        
        # Set the appropriate flags, generate unique strand ids, and store
        if t is RestingState:
            self.use_resting_states = True
            self.boltzmann_sample = True
            
            for resting_state in start_list:
                id_dict = {}
                for strand in resting_state[0].strand_list:
                    id_dict[strand.id] = self.make_unique(strand)
                for cmplx in resting_state:
                    cmplx.strand_list = [id_dict[s.id] for s in cmplx.strand_list]
        
        else:
            self.use_resting_states = False
            
            for cmplx in start_list:
                cmplx.strand_list = [self.make_unique(s) for s in cmplx.strand_list]

    @property
    def initial_seed_flag(self):
        return initial_seed != None

    
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
                raise TypeError("All items must be 'StopCondition', not '%s'." % type(item))
        
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
        
        # Set the internal data member
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
    
    


    


def boltzmann_sample(cmplx):
    """Returns a new Complex object with a structure boltzmann sampled based on
    NUPACK's sample function.
    """
    import os, subprocess
    cwd = os.path.abspath(os.curdir)
    prefix = "temp_boltzmann___"
    
    f = open("%s/%s.in" % (cwd, prefix), "w")
    f.write("%d\n" % len(cmplx.strand_list))
    for strand in cmplx.strand_list:
        f.write(strand.sequence + "\n")
    for i in range(len(cmplx.strand_list)):
        f.write("%d " % (i+1))
    f.write("\n")
    f.close()
    
    subprocess.check_call(["/research/src/sample_dist/bin/sample", "-multi", "-material", "dna", "-count", "1", "%s/%s" % (cwd, prefix)], stdout=subprocess.PIPE)
    
    f = open("%s/%s.sample" % (cwd, prefix), "r")
    for i in range(11):
        line = f.readline()
    f.close()
    
    return Complex(cmplx.id, cmplx.name, cmplx.strand_list, line.strip())







