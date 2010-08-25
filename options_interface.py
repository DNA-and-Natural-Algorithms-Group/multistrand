class Interface(object):
    def __init__(self):
        ####################
        #
        # BEGIN pythondata
        #
        ####################
        
        # TODO: change some of these to use properties and possibly rename
        
        # self.trajectory_time = None
        # """ The total elapsed time of the most recently completed trajectory.
        
        # Type         Default
        # double       None
        
        # Set by ssystem.cc when a trajectory completes.
        # """
        
        self.current_tag = None
        """ The tag of the stop state of the most recently completed trajectory.
        
        Type         Default
        string       None
        
        Set by ssystem.cc when a trajectory completes.
        """

        self.current_completion_type = None
        """ The completion type of the most recent stop state of the trajectory.
        
        Type         Default
        long         None
        
        Set by ssystem.cc when a trajectory completes.
        """
        
        self.collision_rate = None
        """ The collision rate of the most recently completed trajectory.
        
        Type         Default
        double       None
        
        Set by ssystem.cc when a trajectory completes.
        """
        
        self.current_seed = None
        """ The seed used by the random number generator in the most recently
        completed trajectory.
        
        Type         Default
        long         None
        
        Set by MultistrandOptions at the end of each trajectory.
        """
        
        self.current_time = 0.0
        """ The current elapsed time of the currently running trajectory.
        
        Type         Default
        double       0.0
        
        Set by ssystem.cc at certain points during the simulation.
        """
        
        self.trajectory_completion_flag = False
        """ Indicates whether the trajectory has completed.
        
        Type         Default
        boolean      False
        
        Should be set by MultistrandOptions when a trajectory completes.
        """
        
        self.trajectory_halt_flag = False
        """ Indicates whether Multistrand was told to halt its trajectory by an
        external Python program.
        
        Type         Default
        boolean      False
        
        Read by ssystem.cc. There should probably be a function that sets it.
        """
        
        self.trajectory_suspend_flag = False
        """ Indicates whether Multistrand was told to suspend its trajectory by 
        an external Python program.
        
        Type         Default
        boolean      False
        
        Read by ssystem.cc. There should probably be a function that sets it.
        """

        self._trajectory_count = 0

    def __str__(self):
        res = ""
        res = res + "Most recent trajectory time: {0}\n\
        Most recent trajectory tag: {1}\n\
        Most recent trajectory completion type: {2}\n\
        Most recent collision rate: {3}\n\
        Most recent random seed: {4}\n\
        Flags: Completed [{5}] Halt [{6}] Suspend [{7}]\n\
        # of trajectories completed: {8}".format( self.current_time, self.current_tag,
                                                  self.current_completion_type,
                                                  self.collision_rate, self.current_seed,
                                                  self.trajectory_completion_flag, self.trajectory_halt_flag,
                                                  self.trajectory_suspend_flag, self.trajectory_count )
        return res

    @property
    def trajectory_count( self ):
        return self._trajectory_count
    @trajectory_count.setter
    def trajectory_count( self, val ):
        if val == 0:
            self._trajectory_count = 0
        else:
            raise ValueError("Current trajectory count can only be reset to 0, not arbitrarily set.")
        
    def increment_trajectory_count(self):
        self._trajectory_count = self._trajectory_count + 1
