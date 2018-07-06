
from multistrand.objects import Strand, Complex, Domain, StopCondition
from multistrand.options import Literals, standardOptions
from multistrand.system import SimSystem
from multistrand.experiment import makeComplex, hybridization, dissociation


def test0():

    o1 = standardOptions();
    dissociation("CCCCCATTAAC");
    o1.simulation_mode = Literals.first_passage_time
    
    return o1


def create_test3():
    
    o1 = standardOptions();
    hybridization("CCC");
    o1.simulation_mode = Literals.first_passage_time
    
    return o1


def main():
        
    b = Builder(test0, None)
    print b
        
#     o1 = create_test0()  # just a fully hybridized strand.
#     o1 = create_test3()  # this is the  bi-molecular test

# #     setArrheniusConstantsNM1(o1)
#     o1.activestatespace = True
#     o1.output_interval = 1
#     s = SimSystem(o1)
#     # s.start()
#     s.localTransitions()
#     # s.calculate_rate
# #     print "Testing loop internals"
#     
#     # s.start()
#     # s.InitializeSystem()


main()

