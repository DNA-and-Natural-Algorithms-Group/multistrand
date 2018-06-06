from multistrand.builder import hybridizationString

import sys


# # The actual main method
if __name__ == '__main__':

    print sys.argv

    for stateState in hybridizationString("ACTG"):
        for complex in stateState:
            print complex