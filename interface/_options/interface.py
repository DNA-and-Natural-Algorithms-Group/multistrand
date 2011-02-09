from constants import _OptionsConstants

Constants = _OptionsConstants()

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
        
        self.start_structures = {}
        """ A dictionary storing the start structures for each trajectory. The
        keys are trajectory seeds and each value is a list of the starting 
        structures for each complex of the trajectory associated with that seed.
        
        Type         Default
        dict         {}
        
        Modified by MultistrandOptions whenever the current seed is set (which
        should occur immediately after the initial structures are chosen). Items
        are removed whenever a trajectory is completed.
        """
        
        self._trajectory_count = 0
        # Current number of trajectories completed, is an internal that gets incremented
        # by the simsystem as it completes trajectories.

        self._results = ResultList([])
        # hidden member that has a list of Result objects.
        
    @property
    def results( self ):
        return self._results

    def add_result( self, val, res_type= None ):
        if res_type == "status_line":
            seed, com_type, time, tag = val
            start = self.start_structures[seed]
            del self.start_structures[seed]
            new_result = Result( value_list=val, result_type=res_type, start_state=start )
        else:
            seed, com_type, time, rate, tag = val
            start = self.start_structures[seed]
            del self.start_structures[seed]
            new_result = FirstStepResult( value_list = val, start_state = start)
        self._results.append( new_result )

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

    def __init__(self, value_list=None, result_type=None, start_state=None):
        """ Takes the input list and parses it into something intelligible. """
        if result_type == 'status_line':
            self.seed, self.com_type, self.time, self.tag = value_list
            self.result_type = result_type
            self.start_state = start_state
            
        else:
            self.seed = -1
            self.com_type = Constants.STOPRESULT['Error']
            self.time = -1.0
            self.tag = "Result that was not initialized properly."
            self.result_type = result_type

        try:
            self.type_name = Constants.STOPRESULT_inv[ self.com_type ]
        except KeyError:
            self.type_name = "Invalid Type: {0}".format( self.com_type )

    def __str__( self ):
        res = "Trajectory Seed [{0}]\n\
        Result: {1}\n\
        Completion Time: {2}\n\
        Completion Tag: {3}".format( self.seed, self.type_name, self.time, self.tag )
        
        return res

    def format( self, columns ):
        """ Returns a list of formatted strings according to columns. """
        res = []
        for k,fmtstring in columns:
            res.append("{0:{1}}".format( (hasattr( self,k) and getattr( self, k ) ) or "N/A",
                                     fmtstring ))
        return res

    def __repr__( self ):
        return "({0}, {1}, {2}, '{3}', result_type='status_line' )".format(self.seed, self.com_type, self.time, self.tag )


class FirstStepResult( object ):
    """ Holds the results of a single trajectory run by multistrand.

    How to look at the data of a Result object 'r':

    str(r):  printable representation of the data (e.g. print(r) looks nice)
    repr(r): raw data - the python tuple representing this Result object."""

    def __init__(self, value_list=None, start_state=None):
        """ Takes the input list and parses it into something intelligible. """
        self.seed, self.com_type, self.time, self.collision_rate, self.tag = value_list
        self.start_state = start_state
        
        try:
            self.type_name = Constants.STOPRESULT_inv[ self.com_type ]
        except KeyError:
            self.type_name = "Invalid Type: {0}".format( self.com_type )

    def __str__( self ):
        res = "Trajectory Seed [{0.seed}]\n\
        Result: {0.type_name}\n\
        Completion Time: {0.time}\n\
        Collision Rate: {0.collision_rate}\n\
        Completion Tag: {0.tag}".format( self )
        
        return res

    def format( self, columns ):
        """ Returns a list of formatted strings according to columns. """
        res = []
        for k,fmtstring in columns:
            res.append("{0:{1}}".format( (hasattr( self,k) and getattr( self, k ) ) or "N/A",
                                     fmtstring ))
        return res

    def __repr__( self ):
        return "({0}, {1}, {2}, {4}, '{3}', result_type='firststep' )".format(self.seed, self.com_type, self.time, self.tag, self.collision_rate)


class ResultList( list ):
    """ Wrapper class to print a list of results nicely. """
    def __init__( self, *args, **kargs ):
        super(ResultList,self).__init__( *args, **kargs )

    def __str__( self ):
        res = "Current List of Results:\n"
        res += ("-" * 70) + "\n"
        res += "{0:<10} | {1:>12} | {2:<25} [{3:<8}]\n".format('RNG Seed', 'Stop Time', 'Stopping Condition Name', 'Condition Type')
        for item in self:
            data = item.format([
                ('seed',  '<#08x'),
                ('time',  '>12e'),
                ('tag',   '<25'),
                ('type_name', '<8')
                ])
            res += "{0[0]:<10} | {0[1]} | {0[2]} [{0[3]}]\n".format( data )
        return res
