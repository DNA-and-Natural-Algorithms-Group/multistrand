import warnings

class Domain(object):
  """Represents a Multistrand Domain object."""

  def __init__(self, id, name, length, is_complement=False):
    self.id = id
    self.name = name
    self.length = length
    self.is_complement = is_complement


class Strand(object):
  """Represents a Multistrand Strand object."""
  unique_id = 0
  
  def __init__(self, *args, **kargs ):
    """ Initializes a new Strand object. """
    if len(args) == 4:
      warnings.warn( DeprecationWarning("Passing the strand ID is deprecated, it's a private matter internal to Multistrand.") )
      self.name,self.sequence,self.domain_list = args[1:4]
    elif len(args) == 3:
      self.name,self.sequence,self.domain_list = args[0:3]
    elif 'sequence' in kargs:
      self.name = "Automatic_" + str(Strand.unique_id)
      self.sequence = kargs['sequence']
      self.domain_list = []
    else:
      raise ValueError("Strands must be initialized with 3 arguments: name, sequence and domain list.")

    # removed id from the interface possibilities.
    self.id = Strand.unique_id
    Strand.unique_id += 1
  def __str__(self):
    if len(self.sequence) > 30:
      return "Strand:         Name:'{0.name}'\n      : Sequence [{2}]:{0.sequence}\n             : Domains:{1}".format( self, [i.name for i in self.domain_list], len(self.sequence) )
    else:
      return "Strand: Name:'{0.name}' Sequence [{2}]:{0.sequence} Domains:{1}".format( self, [i.name for i in self.domain_list], len(self.sequence) )

class RestingState(tuple):
  """Represents a resting state, i.e. a named set of complexes that exists as a
  strongly connected component with no outward fast transitions in the reaction
  graph."""

  def __new__(cls, *args):
    return tuple.__new__(cls, args[-1])
  
  def __init__(self, name, complex_set, boltzmann_sample=False):
    self.name = name
    self._boltzmann_sample = boltzmann_sample

  def __len__(self):
    return super(RestingState,self).__len__()

  @property
  def boltzmann_sample( self ):
    """ This property indicates whether or not this resting state should be Boltzmann sampled.

    Setting it to true implies that every time any contained complex
    is queried for the 'structure' property, it will return a sampled
    structure.
    """
    return self._boltzmann_sample

  @boltzmann_sample.setter
  def boltzmann_sample(self, value ):
    self._boltzmann_sample = value
    # we are a tuple of complexes, iterate over those to set boltzmann flags.
    for i in self:
      i.boltzmann_sample = value

  def get_starting_complex( self ):
    """ Retreives a suitable complex to use as a start condition -> it
    must have sequence and structure defined, whether it generates
    it from a boltzmann distribution or have a fixed sequence
    structure.

    The current form always returns the 1st complex in the resting
    state and lets that complex do any sub-sampling (boltzmann or
    fixed) as previously set. In the future we probably want to define
    how to sample out of a resting state that contains more than one
    possible form.
    """
    return self[0]
  
class StopCondition(object):
  """Represents a trajectory stopping condition."""

  def __init__(self, tag, complex_items):
    self.tag = tag
    self.complex_items = complex_items  # List of (complex, stoptype, count) tuples

