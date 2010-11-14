import warnings

class StopCondition(object):
  """Represents a trajectory stopping condition."""

  def __init__(self, tag, complex_items):
    self.tag = tag
    self.complex_items = complex_items  # List of (complex, stoptype, count) tuples

