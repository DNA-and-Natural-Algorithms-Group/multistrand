################################################################################
#                                                                              #
# Important Constants for Multistrand.                                         #
# Copyright 2010-217 Caltech                                                   #
# Written by:  Joseph Schaeffer.  Frits Dannenberg                             #
#                                                                              #
#                                                                              #
#                                                                              #
################################################################################



class OptionsConstants( object ):
    def __init__(self):
        self.ZERO_C_IN_K = 273.15
        pass
    
    @property
    def SUBSTRATE_TYPE(self):
        return {"Invalid":0, "RNA":1, \
                "DNA":2}

    @property
    def SUBSTRATE_TYPE_inv(self):
        return {0:"Invalid", 1:"RNA", \
                2:"DNA"}

    @property
    def SIMULATION_MODE(self):
        return {"Normal":                   0x0010,
                "First Step":               0x0030,
                "Transition":               0x0100,
                "Trajectory":               0x0080,
                "First Passage Time":       0x0010}


    @property
    def STOPRESULT(self):
        return {"Normal":                   0x0011,  # 17
                "Time":                     0x0012,  # 18
                "Forward":                  0x0021,
                "Time (First Step)":        0x0022,
                "Reverse":                  0x0024,
                "Error":                    0x0081,
                "NaN":                      0x0082,
                "No Moves":                 0x0084}

    @property
    def STOPRESULT_inv(self):
        return {0x0011:                   "Normal",
                0x0012:                     "Time",
                0x0021:                  "Forward",
                0x0022:        "Time (First Step)",
                0x0024:                  "Reverse",
                0x0081:                    "Error",
                0x0082:                      "NaN",
                0x0084:                 "No Moves"}


    def __setattr__(self, name, value):
        if hasattr(self, name):
            pass
        else:
            object.__setattr__(self, name, value)
    
    def __delattr__(self, *args, **kargs):
        pass




Constants = OptionsConstants()