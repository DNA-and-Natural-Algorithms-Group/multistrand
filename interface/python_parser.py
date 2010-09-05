#
#!/usr/bin/python2.6
#
#

import os
from util_regexp import NamedRE as NamedRE



class MultistrandInputFile( file ):
    def __init__(self, *args, **kargs):
        super(MultistrandInputFile,self).__init__(*args, **kargs)

    def __str__( self ):
        if self.closed:
            raise ValueError("ERROR: MultistrandInputFile Parser: Cannot print from a closed file.")
        self.seek(0, os.SEEK_SET)
        # seek to start of file

        res = ""
        for line in self:
            res = res + str(line)
        return res


class MultistrandInputFile_BaseRegexps( object ):
    """ This class should only ever be instantiated once.
    """

    unique_count = 0

    def __init__(self):
        if unique_count > 0:
            raise UserWarning("There only be one MultistrandInputFile_BaseRegexps object.")
        super(MultistrandInputFile_BaseRegexps,self).__init__()
        
        unique_count = 1

        self.regexp_names = ["_regexp_comment","_regexp_strand","_regexp_startst","_regexp_stopst","_regexp_onelines"]
        self.regexp_list = {}

        for i in self.regexp_names:
            if i in dir(self):
                self.__getattribute__( i )()

        # done setting the object up, if it can't find the attribute.

    def _regexp_comment( self ):
        
            
            

