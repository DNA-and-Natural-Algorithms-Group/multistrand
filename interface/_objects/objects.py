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


class StopCondition(object):
  """Represents a trajectory stopping condition."""
  
  def __init__(self, tag, complex_items):
    self.tag = tag
    self.complex_items = complex_items  # List of (complex, stoptype, count) tuples














