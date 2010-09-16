from .options import Constants

class Interface(object):
    def __init__(self):
        """ Sets some default values for the Interface, in addition to
        initializing the current results list, etc."""
               
        self.current_seed = None
        """ The seed used by the random number generator in the most recently
        completed trajectory.
        
        Type         Default
        long         None
        
        Set by MultistrandOptions at when it begins each trajectory
        [is repeated in the results it prints at the end].
        """
        
        self._trajectory_count = 0
        # Current number of trajectories completed, is an internal that gets incremented
        # by the simsystem as it completes trajectories.

        self._results = []
        # hidden member that has a list of Result objects.
        
    @property
    def results( self ):
        return self._results

    def add_result( self, val, type_hint = None ):
        self._results.append( Result( *val, type_hint = type_hint ) )

    def __str__(self):
        res = "# of trajectories completed: {0}\n\
        Most recent trajectory information:\n{1}".format( self.trajectory_count, str( self._results[-1] ))
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
        self._trajectory_count += 1
    

class Result( object ):
    """ Holds the results of a single trajectory run by multistrand.

    How to look at the data of a Result object 'r':

    str(r):  printable representation of the data (e.g. print(r) looks nice)
    repr(r): raw data - the python tuple representing this Result object."""

    def __init__(self, *args, type_hint=None ):
        """ Takes the input list and parses it into something intelligible. """
        if type_hint == 'status_line':
            self.seed, self.com_type, self.time, self.tag = args
        else:
            self.seed = -1
            self.com_type = Constants.STOPRESULT['Error']
            self.time = -1.0
            self.tag = "Result that was not initialized properly."

    def __str__( self ):
        type_name = [k for k,v in Constants.STOPRESULT.iteritems() if v == self.com_type]
        if len(type_name) == 0:
            type_name = "Invalid Type: {0}".format( self.com_type )
        else:
            type_name = type_name[0]
        res = "Trajectory Seed [{0}]\n\
        Result: {1}\n\
        Completion Time: {2}\n\
        Completion Tag: {3}".format( self.seed, type_name, self.time, self.tag )
        
        return res

    def __repr__( self ):
        return "({0}, {1}, {2}, {3}, type_hint='status_line' )".format(self.seed, self.com_type, self.time, self.tag )


