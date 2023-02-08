

class StopCondition(object):
    """Represents a trajectory stopping condition.  Plug-and-play equivalent to a Macrostate."""
    
   
                 
    def __init__(self, tag, complex_items):
        """ __init__(self, tag, complex_items)
            tag is the reported name for this stop condition (string)
            complex_items is a list of 3-tuples, which are
                (complex, stoptype, count)
                where complex is a Complex object,
                stoptype is one of:
                   0 - ("exact") exact structure stopping condition 
                   1 - ("bound") the strand given in the Complex is bound to some other strand
                   2 - ("dissoc") the strands given in the Complex are connected 
                       [structure does not matter as long as there is a complex in the system with exactly those strands]
                   3 - ("loose") loose structure stopping condition 
                       [*'s in the stopping structure are considered wildcards and anything in those positions is allowable, 
                        all other items must be correctly paired, with at most 'count' bases deviating]
                   4 - ("count") count-based structure stopping condition
                       [the count parameter is the integer number of base-pairs allowed to differ and still match this structure]
                count is ignored for stoptype 0, 1, and 2.
            The StopCondition is considered to be satisfied if _all_ complex_items tuples are satified.
            For Transition Mode, to indicate that the simulation should actually stop, the tag string should begin with "stop:".
        """
        self.tag = str(tag)
        self.complex_items = complex_items  # List of (complex, stoptype, count) tuples

        """ Need to prevent circular import """
        from ..options import Literals    
        protected = [ Literals.no_initial_moves, Literals.sim_error, Literals.time_out]

        if self.tag in protected:
            raise ValueError('Please do not use a protected simulation result tag as a name for a stopping condition \n Protected are: ' + str(protected))

    def __str__(self):
        str_ret = "Stop Condition, tag:" + self.tag
        for i, val in enumerate(self.complex_items):
            str_ret += "\nSequence  %d: %s" % (i, val[0].sequence) 
            str_ret += "\nStructure %d: %s\n" % (i, val[0].structure) 

        return str_ret


class Macrostate(StopCondition):
    pass
