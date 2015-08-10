####################################################################
#                                                                  #
#  Copyright (c) 2010-2015 California Institute of Technology.     #
#  Distributed under the MIT License.                              #
#  (See accompanying file LICENSE or copy at                       #
#  http://opensource.org/licenses/MIT)                             #
#                                                                  #
####################################################################
#                                                                  #
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)         #
#                                                                  #
#                                                                  #
# Options module to grab the options object into this namespace.   #
#                                                                  #
#                                                                  #
####################################################################

from _options.options   import Options
from _options.constants import Constants
from _options.interface import Result

Options.__module__ = 'multistrand.options'
Result.__module__ = 'multistrand.options'

