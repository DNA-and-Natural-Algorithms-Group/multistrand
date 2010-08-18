

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


class Complex(object):
  """Represents a Multistrand Complex object."""
  
  def __init__(self, id, name, strand_list, structure):
    self.id = id
    self.name = name
    self.strand_list = strand_list
    self.structure = structure


class RestingState(tuple):
  """Represents a resting state, i.e. a named set of complexes that exists as a
  strongly connected component with no outward fast transitions in the reaction
  graph."""
  
  def __new__(typ, name, complex_set):
    return super(RestingState, typ).__new__(typ, complex_set)
  
  def __init__(self, name, complex_set):
    self.name = name


class StopCondition(object):
  """Represents a trajectory stopping condition."""
  
  def __init__(self, tag, complex_items):
    self.tag = tag
    self.complex_items = complex_items  # List of (complex, stoptype, count) tuples














