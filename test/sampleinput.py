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
####################################################################

# #Strands
# S,AGTACGGACACTAGCTGGATCTGAGGATTAGT
# Q,ACTAATCCTCAGATCCAGCTAGTGTC
# T6,ACTAATCCTCAGATCCAGCTAGTGTCCGTACT
# #StartStructure
# S,Q
# ......((((((((((((((((((((((((((_))))))))))))))))))))))))))
# T6
# ................................
# #SimTime=2000000
# #NumSims=1
# #OutputInterval=0
# #StopOptions=2
# #Energymodel=NUPACK_DNA_2_3
# #Temperature = 20.0
# #Concentration = .001
# #StopStructures
# Q
# ..........................
# TAG: finalstop
# ##
# ## end old file

start = [Complex("Starting Gate",
                   [['S', 'AGTACGGACACTAGCTGGATCTGAGGATTAGT', '......(((((((((((((((((((((((((('],
                    ['Q', 'ACTAATCCTCAGATCCAGCTAGTGTC', '))))))))))))))))))))))))))']]),
         Complex("Input","T6","ACTAATCCTCAGATCCAGCTAGTGTCCGTACT")]
stop = StopCondition("Output","Q")
my_options = Options(sim_time=2e6,
                     num_sims=1,
                     parameter_type='Nupack',
                     substrate_type='DNA',
                     temperature=20.0,
                     concentration=1e-6,  # 1 microM, equivalent to the
                                          # old which was in units of mM
                     start_structure = start,
                     stop_conditions = [stop])  #or we can use `stop_condition=stop`

run_system( my_options ) # run_system by default uses the
                        #  standard logging/interfacing settings.

forward = StopCondition("Forward Stop State","Q")
reverse = StopCondition(start, tolerance=2)

firststep_options = my_options.copy()
firststep_options.stop_conditions = [forward,reverse]
firststep_options.simulation_mode = 'First Step'
run_system( firststep_options )

start = [Complex("Starting Gate", ["S","Q"], ["AGTACGGACACTAGCTGGATCTGAGGATTAGT",
                                                        "ACTAATCCTCAGATCCAGCTAGTGTC"],
                                                       ["......((((((((((((((((((((((((((",
                                                        "))))))))))))))))))))))))))"]),        
class my_foo(object):
    """ A very simple markov chain initialization object.

    Usage:
    >>> b = my_foo()
    >>> b.temperature = 320.15
    >>> b.markov_start = [0,15]
    >>> print(b)
    Temperature: 320.15
    Starting States: [0, 15]
    >>> c = my_foo()
    >>> c.markov_start = 1
    >>> print(c)
    Temperature: 310.15
    Starting States: [1]
    >>> c.temperature = 'foo'
    Traceback (most recent call last):
        ...
    ValueError: Temperature should be a floating point value, in degrees K.
    >>> c.temperature = -10
    Traceback (most recent call last):
        ...
    Warning: Your temperature was negative, which is likely an error as we expect degrees K.
    """
    
    def __init__(self):
        self._temperature = 273.15 + 37.0

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        try:
            self._temperature = float(value)
        except ValueError:
            raise ValueError("Temperature should be a floating point value, in degrees K.")
        if self._temperature < 0.0:
            raise Warning("Your temperature was negative, which is likely an error as we expect degrees K.")
        
    @property
    def markov_start( self ):
        """ The start state(s) for this simulation. Should be a list of valid states. """
        return self._start_states

    @markov_start.setter
    def markov_start(self, *args):
        if len(args) == 1:
            try:
                states = iter(args[0])
            except TypeError:
                states = [args[0]]
        else:
            states = args
        self._start_states = [i for i in states]

    def __str__(self):
        return "Temperature: {0._temperature}\nStarting States: {0._start_states}".format(self)
