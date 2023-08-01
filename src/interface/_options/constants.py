# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Important Constants for Multistrand.
"""

class OptionsConstants:

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
