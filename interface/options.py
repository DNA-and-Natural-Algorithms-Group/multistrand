#
#
# Copyright (c) 2010 Caltech. All rights reserved.
# Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
#
#
# Options module to grab the options object into this namespace.
#

from _options.options   import Options
from _options.interface import Result

Options.__module__ = 'multistrand.options'
Result.__module__ = 'multistrand.options'

