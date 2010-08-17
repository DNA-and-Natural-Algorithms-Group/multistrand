

class Domain(object):
  """Represents a Multistrand Domain object. The multi_format attribute
  contains the C++ object."""
  
  def __init__(self, id, name, length, is_complement):
    self.id = id
    self.name = name
    self.length = length
    self.is_complement = is_complement


class Strand(object):
  """Represents a Multistrand Strand object. The multi_format attribute
  contains the C++ object."""
  
  def __init__(self, number, name, sequence, domain_ids):
    self.id = number
    self.name = name
    self.sequence = sequence
    self.domain_ids = domain_ids


class Complex(object):
  """Represents a Multistrand Complex object. The multi_format attribute
  contains the C++ object."""
  
  def __init__(self, id, name, strand_ids, bonds):
    self.id = id
    self.name = name
    self.strand_ids = strand_ids
    self.bonds = bonds


class RestingState(tuple):
  """Represents a resting state, i.e. a named set of complexes that exists as a
  strongly connected component with no outward fast transitions in the reaction
  graph. The multi_format attribute contains the C++ object."""
  
  def __new__(typ, name, complex_set):
    return super(RestingState, typ).__new__(typ, complex_set)
  
  def __init__(self, name, complex_set):
    self.name = name


class StopCondition(object):
  """Represents a trajectory stopping condition."""
  
  def __init__(self, tag, complex_items):
    self.tag = tag
    self.complex_items = complex_items  # List of (complex, stoptype) pairs
















