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
# Objects module to get the useful names into the correct          #
# namespace.                                                       #
#                                                                  #
#                                                                  #
####################################################################


__module__ = 'multistrand.objects'

from _objects.domain  import Domain
from _objects.complex import Complex
from _objects.strand import Strand
from _objects.restingstate import RestingState
from _objects.stopcondition import StopCondition, Macrostate

__all__ = ['Strand','RestingState','StopCondition','Macrostate','Complex','Domain']

# The following appears to be necessary [at the moment] as otherwise
# we can't generate appropriate documentation for these sub objects as
# they look like they are in a different package. Perhaps we should
# change the directory hierarchy - rename _objects to objects, and
# then moving this file into the new directory as __init__.py and
# importing from local space.

Strand.__module__ = 'multistrand.objects'
Complex.__module__ = 'multistrand.objects'
RestingState.__module__ = 'multistrand.objects'
Macrostate.__module__ = 'multistrand.objects'
StopCondition.__module__ = 'multistrand.objects'
Domain.__module__ = 'multistrand.objects'
