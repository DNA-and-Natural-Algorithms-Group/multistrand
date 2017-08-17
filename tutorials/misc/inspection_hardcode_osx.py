
from multistrand.objects import Strand, Complex, Domain, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem

def createOptions(start_complex, stop_complex, simMode):

    full_sc = StopCondition("CLOSED", [(stop_complex, Options.dissocMacrostate, 2)])    
    
    o1 = Options(simulation_mode=simMode,  # "First Passage Time", 
                 parameter_type="Nupack", substrate_type="DNA", temperature=273.15 + 25.0,
                num_simulations=10,
                simulation_time=0.00001,
                rate_scaling='Calibrated',
                verbosity=0,
                join_concentration=1e-9,
                rate_method="Metropolis",
                start_state=[start_complex],
                stop_conditions=[full_sc])
   

    
    return o1

    

def makeComplex(seq, dotparen):
    
    strandList = []
    
    for seq in seq:
    
        onedomain = Domain(name="domain" + str(makeComplex.counter), sequence=seq)
        makeComplex.counter += 1;

        onestrand = Strand(domains=[onedomain])
        strandList.append(onestrand)
    
    return Complex(strands=strandList, structure=dotparen)

makeComplex.counter = 0

def create_test10():


    seq0 = "GTCACTGCTTTT"
    seq1 = "GCAGTGAC"
    dotparen1 = "..((((((....+)))))).."
    
    
    seq2 = "GTCACTGC"
    dotparen2 = "........"
    
    complex1 = makeComplex([seq0, seq1], dotparen1)
    complex2 = makeComplex([seq2], dotparen2 )
        
    return createOptions(complex1, complex2, "First Passage Time")



def main():
    

    o1 = create_test10()      # half open duplex

    s = SimSystem(o1)
    s.initialInfo()


main()








