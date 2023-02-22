from multistrand.objects import *
from multistrand.options import Options
from multistrand.system import *

TEMPERATURE=25
o = Options(temperature=TEMPERATURE,dangles="Some", substrate_type="RNA")
initialize_energy_model(o)
