################################################################################
#                                                                              #
# Important Constants for Multistrand.                                         #
# Copyright 2010 Caltech                                                       #
# Written by:  Joseph Schaeffer.                                               #
#                                                                              #
#                                                                              #
#                                                                              #
################################################################################


class _OptionsConstants( object ):
    def __init__(self):
        self.ZERO_C_IN_K = 273.15
        pass
    
    @property
    def RATEMETHOD(self):
        return {"Invalid" :0, "Metropolis"    :1, \
                "Kawasaki":2, "EntropyEnthalpy":3}
    @property
    def RATEMETHOD_inv(self):
        return {0:"Invalid",1:"Metropolis", \
                2:"Kawasaki",3:"EntropyEnthalpy"}
        
    @property
    def DANGLES(self):
        return {"None":  0, "Some" : 1, \
                "All" :  2, "NupackDefault": 1}

    @property
    def DANGLES_inv(self):
        return {0:"None", 1:"Some", \
                2:"All" }

    @property
    def ENERGYMODEL_TYPE(self):
        return {"Vienna":0, "Nupack":1, \
                "Others?":2}
    
    @property
    def ENERGYMODEL_TYPE_inv(self):
        return {0:"Vienna", 1:"Nupack", \
                2:"Others?"}

    @property
    def energy_type( self ):
        pass

    @property
    def macrostate_type( self ):
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
                "Python Module":            0x0040,
                "Python Module:First Step": 0x0060,
                "Energy Only":              0x0200,
                "Transition":               0x0100,
                "Trajectory":               0x0080,
                "First Passage Time":       0x0010}

#     @property
#     def SIMULATION_MODE_FLAG(self):
#         return {"Normal":                   0x0010,
#                 "First Bimolecular":        0x0020,
#                 "Python Module":            0x0040,
#                 "Trajectory":               0x0080,
#                 "Transition":               0x0100,
#                 "First Passage Time":       0x0010}

    @property
    def STOPRESULT(self):
        return {"Normal":                   0x0011,
                "Time":                     0x0012,
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

#     @property
#     def rate_scaling_sets(self):
#         return {"Bimolecular_Calibrate": ## WETMUR unimolecular parameters
#                 {"description":   "Calibration data set for determining bimolecular parameters. Uses the Wetmur unimolecular calibration.\n Unimolecular: {0}  Bimolecular: {1}. Concentration: {3:e}M\n Model: [{2[1]}] Substrate: [{2[0]}] Dangles: [{2[2]}] Rate Method: [{2[3]}]\n",
#                  "DNA:Nupack:None:Kawasaki"  :    {"uni": 2.1e8,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:Some:Kawasaki"  :    {"uni": 1.5e8,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:All:Kawasaki"   :    {"uni": 1.5e8,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:None:Metropolis":    {"uni": 6.8e8,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:Some:Metropolis":    {"uni": 7.3e8,
#                                                    "bi":  1.0e3},
#                  "DNA:Nupack:All:Metropolis" :    {"uni": 7.3e8,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:None:Kawasaki:298.15"  :    {"uni": 9.5e7,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:Some:Kawasaki:298.15"  :    {"uni": 6.1e7,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:All:Kawasaki:298.15"   :    {"uni": 6.1e7,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:None:Metropolis:298.15":    {"uni": 4.2e8,
#                                                    "bi": 1.0e3},
#                  "DNA:Nupack:Some:Metropolis:298.15":    {"uni": 4.4e8,
#                                                    "bi":  1.0e3},
#                  "DNA:Nupack:All:Metropolis:298.15" :    {"uni": 4.4e8,
#                                                    "bi": 1.0e3},
# 
#                  "default": {"uni": 7.3e8,
#                              "bi": 1.0e3}
#                  },
#         "Default": ## Calibrated via Wetmur/Morrison  # should be same as "Calibrated"
#                 {"description":   "Calibrated parameters. Uses the Wetmur unimolecular calibration.\n Unimolecular: {0}  Bimolecular: {1}. Concentration: {3:e}M\n Model: [{2[1]}] Substrate: [{2[0]}] Dangles: [{2[2]}] Rate Method: [{2[3]}] Temperature: [{2[4]:.2f}K]\n",
#                  "DNA:Nupack:None:Kawasaki"  :    {"uni": 2.1e8,
#                                                    "bi":  8.07e5},
#                  "DNA:Nupack:Some:Kawasaki"  :    {"uni": 1.5e8,
#                                                    "bi": 1.38e6},
#                  "DNA:Nupack:All:Kawasaki"   :    {"uni": 1.5e8,
#                                                    "bi": 1.38e6},
#                  "DNA:Nupack:None:Metropolis":    {"uni": 6.8e8,
#                                                    "bi": 8.18e5},
#                  "DNA:Nupack:Some:Metropolis":    {"uni": 7.3e8,
#                                                    "bi":  1.40e6},
#                  "DNA:Nupack:All:Metropolis" :    {"uni": 7.3e8,
#                                                    "bi":  1.41e6},
#                  "DNA:Nupack:None:Kawasaki:298.15"  :    {"uni": 9.5e7,
#                                                    "bi": 6.07e5},
#                  "DNA:Nupack:Some:Kawasaki:298.15"  :    {"uni": 6.1e7,
#                                                    "bi": 1.29e6},
#                  "DNA:Nupack:All:Kawasaki:298.15"   :    {"uni": 6.1e7,
#                                                    "bi": 1.28e6},
#                  "DNA:Nupack:None:Metropolis:298.15":    {"uni": 4.2e8,
#                                                    "bi": 6.02e5},
#                  "DNA:Nupack:Some:Metropolis:298.15":    {"uni": 4.4e8,
#                                                    "bi": 1.26e6},
#                  "DNA:Nupack:All:Metropolis:298.15" :    {"uni": 4.4e8,
#                                                    "bi": 1.29e6},
# 
#                  "default": {"uni": 1.5e8,
#                              "bi": 1.38e6}
#                  },
#         "Calibrated": ## Calibrated via Wetmur/Morrison  
#                 {"description":   "Calibrated parameters. Uses the Wetmur unimolecular calibration.\n Unimolecular: {0}  Bimolecular: {1}. Concentration: {3:e}M\n Model: [{2[1]}] Substrate: [{2[0]}] Dangles: [{2[2]}] Rate Method: [{2[3]}] Temperature: [{2[4]:.2f}K]\n",
#                  "DNA:Nupack:None:Kawasaki"  :    {"uni": 2.1e8,
#                                                    "bi":  8.07e5},
#                  "DNA:Nupack:Some:Kawasaki"  :    {"uni": 1.5e8,
#                                                    "bi": 1.38e6},
#                  "DNA:Nupack:All:Kawasaki"   :    {"uni": 1.5e8,
#                                                    "bi": 1.38e6},
#                  "DNA:Nupack:None:Metropolis":    {"uni": 6.8e8,
#                                                    "bi": 8.18e5},
#                  "DNA:Nupack:Some:Metropolis":    {"uni": 7.3e8,
#                                                    "bi":  1.40e6},
#                  "DNA:Nupack:All:Metropolis" :    {"uni": 7.3e8,
#                                                    "bi":  1.41e6},
#                  "DNA:Nupack:None:Kawasaki:298.15"  :    {"uni": 9.5e7,
#                                                    "bi": 6.07e5},
#                  "DNA:Nupack:Some:Kawasaki:298.15"  :    {"uni": 6.1e7,
#                                                    "bi": 1.29e6},
#                  "DNA:Nupack:All:Kawasaki:298.15"   :    {"uni": 6.1e7,
#                                                    "bi": 1.28e6},
#                  "DNA:Nupack:None:Metropolis:298.15":    {"uni": 4.2e8,
#                                                    "bi": 6.02e5},
#                  "DNA:Nupack:Some:Metropolis:298.15":    {"uni": 4.4e8,
#                                                    "bi": 1.26e6},
#                  "DNA:Nupack:All:Metropolis:298.15" :    {"uni": 4.4e8,
#                                                    "bi": 1.29e6},
# 
#                  "default": {"uni": 1.5e8,
#                              "bi": 1.38e6}
#                  },
# 
# 
# 
#             "Unitary":
#                 {"description": "Default 1.0/1.0.",
#                  "default": {"uni":1.0,
#                              "bi":1.0}
#                 }
#             }

    def __setattr__(self, name, value):
        if hasattr(self, name):
            pass
        else:
            object.__setattr__(self, name, value)
    
    def __delattr__(self, *args, **kargs):
        pass

Constants = _OptionsConstants()
