"""
This test exercises the Python interface.
"""

import unittest

from multistrand.objects import Domain, Strand, Complex, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem


class Base_Objects_TestCase(unittest.TestCase):
    """
    This test case handles all the basic objects testing for python_objects.py.
    Domain, Strand, Complex, StopCondition and RestingState should all be tested
    here.
    """
    def checkDuplicates(self, seq ):
        """ Helper function to see if we duplicated the 'id' attribute or 'tag' attribute in an interable."""
        if hasattr(seq[0], 'id' ):
            return len(seq) == len( set( [i.id for i in seq] ))
        else:
            return len(seq) == len( set( [i.tag for i in seq] ))

    def test_domains(self):
        """ Test [Objects]: Create several Domain objects

        Currently, it checks to make sure domains are uniquely named after addition, but it may be useful to test other properties as well.
        """
        if hasattr(self,"domains"):
            del self.domains

        self.domains = []
        self.domains.append( Domain("d1", "d", 5, False) )
        self.domains.append( Domain("d1p", "d", 5, True) )

        self.assertTrue( self.checkDuplicates(self.domains), "Duplicate Domains were added.")

    def test_strands(self):
        """ Test [Objects]: Create some Strand objects

        This test creates some strands using the default domains [created by test_domains], makes sure the strands are not duplicate names."""

        if not hasattr(self,"domains"):
            self.test_domains()

        if hasattr(self,"strands"):
            del self.strands

        self.strands = []
            
        self.strands.append( Strand("s1", "s1",  "ACTTG", [self.domains[0]]))
        self.strands.append( Strand("s2", "s2",  "CAAGT", [self.domains[1]]))

        self.assertTrue( self.checkDuplicates(self.strands), "Duplicate strands were added.")

    def test_complexes(self):
        """ Test [Objects]: Create several Complexes

        This test creates some complexes, using the default strands - it makes three for later use as start structures and stop structures.
        Current the only 'test' implemented is making sure they have no duplicate names.
        """

        if not hasattr(self,"strands"):
            self.test_strands()

        if hasattr(self,"complexes"):
            del self.complexes

        self.complexes = []

        self.complexes.append( Complex("c1", "c1", [self.strands[0]], "....."))
        self.complexes.append( Complex("c2", "c2", [self.strands[1]], "....."))
        self.complexes.append( Complex("c3", "c3", [self.strands[0], self.strands[1]], "(((((+)))))"))
                               
        self.assertTrue( self.checkDuplicates(self.complexes), "Duplicate complexes were added.")

    def test_conditions(self):
        """ Test [Objects]: Create two StopConditions [forward + reverse] 

        This test creates some conditions, makes sure they have no duplicate names."""

        if not hasattr(self,"complexes"):
            self.test_complexes()

        if hasattr(self,"conditions"):
            del self.conditions
        self.conditions = []
        
        self.conditions.append( StopCondition("REVERSE", [(self.complexes[0], 2, 0), (self.complexes[1], 2, 0)]))
        self.conditions.append( StopCondition("END", [(self.complexes[2], 4, 1)]))
        
        self.assertTrue (self.checkDuplicates(self.conditions), "Duplicate conditions were added.")


class Options_Object_TestCase(unittest.TestCase):
    """ This test case handles creation of an options object, and should test all functionality in python_options.py (eventually).

    """
    def setUp(self):
        self.domains = []
        self.strands = []
        self.complexes = []
        self.conditions = []
        
        self.domains.append( Domain("d1", "d", 5, False) )
        self.domains.append( Domain("d1p", "d", 5, True) )
        self.strands.append( Strand("s1", "s1",  "ACTTG", [self.domains[0]]))
        self.strands.append( Strand("s2", "s2",  "CAAGT", [self.domains[1]]))
        self.complexes.append( Complex("c1", "c1", [self.strands[0]], "....."))
        self.complexes.append( Complex("c2", "c2", [self.strands[1]], "....."))
        self.complexes.append( Complex("c3", "c3", [self.strands[0], self.strands[1]], "(((((+)))))"))
        self.conditions.append( StopCondition("REVERSE", [(self.complexes[0], 2, 0), (self.complexes[1], 2, 0)]))
        self.conditions.append( StopCondition("END", [(self.complexes[2], 4, 1)]))

    def tearDown(self):
        self.domains[:] = []
        self.strands[:] = []
        self.complexes[:] = []
        self.conditions[:] = []

    def test_options(self):
        """ Test [Options]: Create an options object using some default values and start/stop states"""
        o = Options()
        o.simulation_mode = 3
        o.use_stop_states = True
        o.parameter_type  = 1
        o.substrate_type  = 2
        o.num_simulations = 3
        o.simulation_time = 0.5
        o.start_state = [self.complexes[0], self.complexes[1]]
        o.stop_conditions = [self.conditions[0], self.conditions[1]]
        #Do some tests?
        # For the moment, if creation doesn't fail, we're happy. :)


class System_Object_TestCase(unittest.TestCase):
    """ This test case handles creation of an multistrand SimulationSystem object (via Boost), and then tries a few test cases to see if they work.

    """
    def setUp(self):
        self.domains = []
        self.strands = []
        self.complexes = []
        self.conditions = []
        
        self.domains.append( Domain("d1", "d", 5, False) )
        self.domains.append( Domain("d1p", "d", 5, True) )
        self.strands.append( Strand("s1", "s1",  "ACTTG", [self.domains[0]]))
        self.strands.append( Strand("s2", "s2",  "CAAGT", [self.domains[1]]))
        self.complexes.append( Complex("c1", "c1", [self.strands[0]], "....."))
        self.complexes.append( Complex("c2", "c2", [self.strands[1]], "....."))
        self.complexes.append( Complex("c3", "c3", [self.strands[0], self.strands[1]], "(((((+)))))"))
        self.conditions.append( StopCondition("REVERSE", [(self.complexes[0], 2, 0), (self.complexes[1], 2, 0)]))
        self.conditions.append( StopCondition("END", [(self.complexes[2], 4, 1)]))

        self.options = Options()
        self.options.simulation_mode = 3
        self.options.use_stop_states = True
        self.options.parameter_type  = 1
        self.options.substrate_type  = 2
        self.options.num_simulations = 3
        self.options.simulation_time = 0.5
        self.options.start_state = [self.complexes[0], self.complexes[1]]
        self.options.stop_conditions = [self.conditions[0], self.conditions[1]]

    def tearDown(self):
        self.domains[:] = []
        self.strands[:] = []
        self.complexes[:] = []
        self.conditions[:] = []
        self.options = None

    def test_create_system(self):
        """ Test [System]: Create a simulation system object

        Yep, that's it."""
        system = SimSystem(self.options)

    def test_create_system_repeated(self):
        """ Test [System]: Create three distinct system objects

        Note that this uses the /same/ options object for each, which
        is theoretically a really bad idea - an options object should
        be tied to a system in future versions."""
        
        system1 = SimSystem(self.options)
        system2 = SimSystem(self.options)
        system3 = SimSystem(self.options)

    def test_run_system(self):
        """ Test [System]: Create a system object and then run the system

        This test creates a SimulationSystem and then runs it, printing the options object after completion."""
        system = SimSystem(self.options)
        system.start()

        System_Object_TestCase.str_run_system = str(self.options.interface)

    def test_run_system_several_times(self):
        """ Test [System]: Create three system objects, then run each in sequence

        We run a system several times (using the same options object, but in sequence)
        Needs changing in the future."""
        system = SimSystem(self.options)
        system.start()

        System_Object_TestCase.str_run_system_several_times = "First run results:\n{0}\n".format(str(self.options.interface))
        
        system2 = SimSystem(self.options)
        system2.start()

        System_Object_TestCase.str_run_system_several_times += "Second run results [different system]:\n{0}\n".format(str(self.options.interface))

        system3 = SimSystem(self.options)
        system3.start()
        
        System_Object_TestCase.str_run_system_several_times += "Third run results [yet another system]:\n{0}\n".format(str(self.options.interface))


class SetupSuite( object ):
    """ Container for default set of tests and standard method for running them."""

    def __init__(self):
        self._suite = unittest.TestSuite()
        self._suite.addTests(
            unittest.TestLoader().loadTestsFromTestCase(
                Base_Objects_TestCase ))
        # the basic object test cases.
        self._suite.addTests(
            unittest.TestLoader().loadTestsFromTestCase(
                Options_Object_TestCase ))

    def runTests(self):
        if hasattr(self, "_suite") and self._suite is not None:
            unittest.TextTestRunner(verbosity=2).run( self._suite )


if __name__ == '__main__':
    suite = SetupSuite()
    suite.runTests()
