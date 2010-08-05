

class Options(object):
  
  def __init__(self):
    self.name = "Multistrand Options object"
    self.integer = 39999
    self.decimal = 9.1
    self.list_of_strings = ["ACT", "CAT", "GAC", "TAG"]
  
  def no_args_no_return(self):
    print "This returns nothing."
  
  def one_arg_no_return(self, number):
    print "You gave me %.1f" % number

'''  

def multiply(a,b):
  print "Will compute", a, "times", b
  c = 0
  for i in range(0, a):
    c = c + b
  return c

'''
