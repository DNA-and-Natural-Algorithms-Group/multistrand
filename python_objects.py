

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
  
  @property
  def sequence(self):
    return "+".join([strand.sequence for strand in self.strand_list])
  
  @property
  def boltzmann_structure(self):
    import os, subprocess, random
    cwd = os.path.abspath(os.curdir)
    prefix = "temp_boltzmann___" + str(random.random())
    
    f = open("%s/%s.in" % (cwd, prefix), "w")
    f.write("%d\n" % len(self.strand_list))
    for strand in self.strand_list:
        f.write(strand.sequence + "\n")
    for i in range(len(self.strand_list)):
        f.write("%d " % (i+1))
    f.write("\n")
    f.close()
    
    subprocess.check_call(["/research/src/sample_dist/bin/sample", "-multi", "-material", "dna", "-count", "1", "%s/%s" % (cwd, prefix)], stdout=subprocess.PIPE)
    
    f = open("%s/%s.sample" % (cwd, prefix), "r")
    for i in range(11):
        line = f.readline()
    f.close()
    os.remove("%s/%s.in" % (cwd, prefix))
    os.remove("%s/%s.sample" % (cwd, prefix))
    
    sampled_structure = line.strip()    
    return sampled_structure


class RestingState(tuple):
  """Represents a resting state, i.e. a named set of complexes that exists as a
  strongly connected component with no outward fast transitions in the reaction
  graph."""
  
  def __new__(cls, *args):
    return tuple.__new__(cls, args[-1])
  
  def __init__(self, name, complex_set):
    self.name = name


class StopCondition(object):
  """Represents a trajectory stopping condition."""
  
  def __init__(self, tag, complex_items):
    self.tag = tag
    self.complex_items = complex_items  # List of (complex, stoptype, count) tuples














