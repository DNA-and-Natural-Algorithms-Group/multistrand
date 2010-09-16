class Domain(object):
  """Represents a Multistrand Domain object."""

  def __init__(self, id, name, length, is_complement=False):
    self.id = id
    self.name = name
    self.length = length
    self.is_complement = is_complement


class Strand(object):
  """Represents a Multistrand Strand object."""

  def __init__(self, id, name, sequence, domain_list):
    self.id = id
    self.name = name
    self.sequence = sequence
    self.domain_list = domain_list


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

