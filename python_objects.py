

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
    # what is this?

    self.name = name
    # This too.

    self.strand_list = strand_list
    # Ok, list of strand names (possibly replicated)

    self.structure   = structure
    # single flat sequence? or list?

    self.unique_ids  = []
    # list of strand's unique ids that map to strand_list's names.
  def __len__(self):
    return len(self.strand_list)
  
  @property
  def sequence(self):
    return "+".join([strand.sequence for strand in self.strand_list])

  @property
  def boltzmann_sequence(self):
    return self.sequence  # For now...


class RestingState(tuple):
  """Represents a resting state, i.e. a named set of complexes that exists as a
  strongly connected component with no outward fast transitions in the reaction
  graph."""

  def __new__(cls, *args):
    return tuple.__new__(cls, args[-1])

  def __init__(self, name, complex_set):
    self.name = name
    self.boltzmann_sample = False

  def __len__(self):
    return super(RestingState,self).__len__()

  def set_boltzmann(self, val):
    self.boltzmann_sample = val
    # possibly update other components.

  def sample(self):
    if self.boltzmann_sample:
      return self.sample_boltzmann()
    else:
      return self[0]

  def sample_boltzmann(self):
    """Returns a new Complex object with a structure boltzmann sampled based on
    NUPACK's sample function.
    """
    import os, subprocess
    cwd = os.path.abspath(os.curdir)
    prefix = "temp_boltzmann___"

    f = open("%s/%s.in" % (cwd, prefix), "w")
    cmplx = self[0]
    f.write("%d\n" % len(cmplx.strand_list))
    for strand in cmplx.strand_list:
      f.write(strand.sequence + "\n")
    for i in range(len(cmplx.strand_list)):
      f.write("%d " % (i+1))
    f.write("\n")
    f.close()

    subprocess.check_call(["/research/src/sample_dist/bin/sample", "-multi", "-material", "dna", "-count", "1", "%s/%s" % (cwd, prefix)], stdout=subprocess.PIPE)

    f = open("%s/%s.sample" % (cwd, prefix), "r")
    for i in range(11):
      line = f.readline()
    f.close()

    return Complex(cmplx.id, cmplx.name, cmplx.strand_list, line.strip())



class StopCondition(object):
  """Represents a trajectory stopping condition."""

  def __init__(self, tag, complex_items):
    self.tag = tag
    self.complex_items = complex_items  # List of (complex, stoptype, count) tuples














