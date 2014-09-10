# multistrand_setup.py
#
# This is really unnecessary on most Multistrand installations -- you can just import the multistrand components and be done with it.
# But using this module instead ("from multistrand_setup import *") might catch a few installation glitches.  So it's useful.

import os, sys, subprocess

multihome = None
if 'MULTISTRANDHOME' in os.environ:
    if not os.path.isfile( os.path.join( os.environ['MULTISTRANDHOME'], 'setup.py') ):
        warnings.warn( ImportWarning("Could not find the file 'setup.py' in your MULTISTRANDHOME [{0}]; this environment variable is possibly out of date or not referring to the new Multistrand distribution."))
    else:
        if os.environ['MULTISTRANDHOME'] not in sys.path:
            multihome= os.environ['MULTISTRANDHOME']     

if multihome != None:
    sys.path.append(multihome)

try:
    from multistrand.objects import *
    from multistrand.options import Options
    from multistrand.system import SimSystem

except ImportError:
    print("Could not import Multistrand.")
    raise

nupackbin = None
if 'NUPACKHOME' in os.environ:
    if not os.path.isfile( os.path.join( os.environ['NUPACKHOME'] + '/bin', 'sample') ):
        warnings.warn( ImportWarning("Could not find the file 'sample' in your NUPACKHOME [{0}] 'bin' directory; this environment variable is possibly out of date or not referring to the latest NUPACK distribution."))
    else:
        if (os.environ['NUPACKHOME']+'/bin') not in sys.path:
            nupackbin= os.environ['NUPACKHOME']+'/bin'

if nupackbin != None:
    sys.path.append(nupackbin)

p = subprocess.Popen(['which','sample'],stdout=subprocess.PIPE)
result,error=p.communicate()
if os.environ['NUPACKHOME'] not in result:
    warnings.warn( ImportWarning("The NUPACK executable 'sample' is not located in system calls; Boltzmann sampling may not work."))

