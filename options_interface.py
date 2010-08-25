class Interface(object):
    def __init__(self):



        ####################
        #
        # BEGIN pythondata
        #
        ####################
        
        # TODO: change some of these to use properties and possibly rename
        
        self.trajectory_time = None
        """ The total elapsed time of the most recently completed trajectory.
        
        Type         Default
        double       None
        
        Set by ssystem.cc when a trajectory completes.
        """
        
        self.trajectory_tag = None
        """ The tag of the stop state of the most recently completed trajectory.
        
        Type         Default
        string       None
        
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
   
