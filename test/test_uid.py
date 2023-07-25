
from multistrand.utils import pairType


class Test_UID:
    """
    This demonstrates the pair-type functionality.
    For a given complex, Pairtype returns a unique representation.
    This is important for hashing functions.
    """
    def test_eq(self):
        names1 = "strand1,strand2"
        names2 = "strand2,strand1"
    
        struct1 = "..((....((..+))....))"
        struct2 = "((....((+..))....)).."
    
        uID1 = pairType(names1, struct1)
        uID2 = pairType(names2, struct2)
        
        print(str(uID1))
        print(str(uID2))
        
        assert uID1 == uID2
        
    def test_neq(self):
        names1 = '11:invader,5:invader*,10:top'
        names2 = '11:invader,5:invader*,10:top'

        struct1 = ".(.+.)..(((+)))."
        struct2 = ".(.+.)..((.+.))."
        
        uID1 = pairType(names1, struct1)
        uID2 = pairType(names2, struct2)
        
        print(str(uID1))
        print(str(uID2))
        
        assert uID1 != uID2
