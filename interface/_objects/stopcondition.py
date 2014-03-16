import warnings

class StopCondition(object):
  """Represents a trajectory stopping condition."""

  def __init__(self, tag, complex_items):
    """ __init__(self, tag, complex_items)
        tag is the reported name for this stop condition (string)
        complex_items is a list of 3-tuples, which are
            (complex, stoptype, count)
            where complex is a Complex object,
            stoptype is one of:
               0 - exact structure stopping condition
               1 - the strand given in the Complex is bound to some other strand
               2 - the strands given in the Complex are connected [structure does not matter as long as there is a complex in the system with exactly those strands]
               3 - loose structure stopping condition (not currently supported in our Complex objects - *'s in the stopping structure are considered wildcards and anything in those positions is allowable, all other items must be correctly paired)
               4 - count-based structure stopping condition: the count parameter is the integer number of base-pairs allowed to differ and still match this structure

    """
    self.tag = str(tag)
    self.complex_items = complex_items  # List of (complex, stoptype, count) tuples


class Macrostate( StopCondition ):
    pass
