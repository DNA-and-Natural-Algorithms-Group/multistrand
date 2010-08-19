from python_objects import *

class Options(object):
  
  def __init__(self):
    self.name = "Multistrand Options object"
    self.decimal = 9.1
    self._integer = 6
    self.neg_integer = -4
    self.list_of_strings = ["ACT", "CAT", "GAC", "TAG"]
    
    strand1 = Strand("S1", "S1", "GATTACA", [Domain("d1", "d1", 7)])
    strand2 = Strand("S2", "S2", "GATTACC", [Domain("d2", "d2", 7)])
    complex1 = Complex("C1", "C1", [strand1, strand2], "(())")
    complex2 = Complex("C2", "C2", [strand1], "(())")
    complex3 = Complex("C3", "C3", [strand2], "(())")
    
    self.stop_conditions = [StopCondition("TAG1", [(complex1, 4, 5), (complex2, 4, 5)]),
                            StopCondition("TAG2", [(complex2, 6, 0), (complex3, 6, 1)])]
    self.start_complexes = [complex1, complex2, complex3]
    
  
  @property
  def integer(self):
    return self._integer

  def no_args_no_return(self):
    print "This returns nothing."
  
  def one_arg_no_return(self, number):
    print "You gave me %.1f" % number

