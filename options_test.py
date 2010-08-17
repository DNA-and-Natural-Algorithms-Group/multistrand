

class Options(object):
  
  def __init__(self):
    self.name = "Multistrand Options object"
    self.decimal = 9.1
    self._integer = 6
    self.list_of_strings = ["ACT", "CAT", "GAC", "TAG"]
    
  @property
  def integer(self):
    return self._integer

  def no_args_no_return(self):
    print "This returns nothing."
  
  def one_arg_no_return(self, number):
    print "You gave me %.1f" % number


