from multistrand.objects import Complex, Domain, Strand
from multistrand.options import Options
from multistrand.system import SimSystem

""" Run this test twice times to check if Multistrand returns a deterministic
    trajectory """

def create_setup(seed):

    # build complexes with domain-level information
    toehold_seq = "GTGGGT"
    bm_design_B = "ACCGCACGTCACTCACCTCG"
    toehold_extra = "TTT"

    toehold = Domain(name="toehold", sequence=toehold_seq, length=6)
    branch_migration_B = Domain(name="bm_B", sequence=bm_design_B, seq_length=20) 

    substrate_B = toehold + branch_migration_B

    incumbent_B = Strand(name="incumbent", domains=[branch_migration_B.C])
    incoming_B = substrate_B.C


    start_complex_B1 = Complex(strands=[incoming_B, substrate_B, incumbent_B],
                              structure="..(.((.....)).).....((((((+))))))((((((((((((((((((((+))))))))))))))))))))")
    
    
    o2 = Options()
    o2.simulation_mode = 0x0080  # trajectory mode
    o2.current_seed = 20
    o2.initial_seed = seed
    o2.num_simulations = 1
    o2.simulation_time = 0.0001  
    o2.temperature = 37.0
    o2.dangles = 1
    o2.start_state = [start_complex_B1]
    o2.output_interval = 1      
    
    o2.JSDefault()
    
     
    return o2

def process_trajectory(o, testing, seed):

    seqstring = ''
    
    times = o.full_trajectory_times
    tubestructs = []
    energies = []
    
    for i in range(len(o.full_trajectory)):  # go through each output microstate of the trajectory

        states = o.full_trajectory[i]  # this is a list of the complexes present in this tube microstate

        structs = []
        for state in states: structs += [ state[4] ]  # similarly extract the secondary structures for each complex
        tubestruct = ' '.join(structs)
        tubestructs.append(tubestruct)  # give the dot-paren secondary structure for the whole test tube

        dG = 0
        for state in states: dG += state[5]
                      
        energies.append(dG);
        
        newseqs = []
        for state in states: newseqs += [ state[3] ]  # extract the strand sequences in each complex (joined by "+" for multistranded complexes)
        newseqstring = ' '.join(newseqs)  # make a space-separated string of complexes, to represent the whole tube system sequence
        
        arrType = o.full_trajectory_arrType[i]
        
        if not newseqstring == seqstring : 
#             print newseqstring
            seqstring = newseqstring  # because strand order can change upon association of dissociation, print it when it changes        
        
        
        
        
    nVal = 18;
    
    if(testing & (seed == 19)):
     
        print(times[-1])
        print(times[-2])
                 
        print(tubestructs[-1])
        print(tubestructs[-2])
     
        print(energies[-1])
        print(energies[-2])
                 
#         print(round(times[-1], nVal) == round(2.00311181181e-06, nVal))
#         print(round(times[-2], nVal) == round(1.99980667462e-06, nVal))
#          
#         print(tubestructs[-1] == "((.(.(....).)...))...(((((+))))).((((((((((((((((((((+))))))))))))))))))))")
#         print(tubestructs[-2] == "((.(.(....).)...))..((((((+))))))((((((((((((((((((((+))))))))))))))))))))")
#      
#         print(round(energies[-1], 2) == round(-31.7, 2))
#         print(round(energies[-2], 2) == round(-32.13, 2))
         
        print 
        

# Perform the simulations




o2 = create_setup(seed=19)
s2 = SimSystem(o2)
s2.start()
# see below about the energy model


# process_trajectory(o2, False)
process_trajectory(o2, True, 19)

