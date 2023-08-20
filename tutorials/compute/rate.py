# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

import sys
from pathlib import Path
sys.path.append(Path(__file__).resolve().as_posix())

import multiprocess


import sys, time


def errorPrint():
    print("Specify the experiment in the first argument: <dissociation> or <hybridization>  ")
    print("Please provide a DNA sequence as the second argument")
    print("Add -bootstrap to do a boostrap ")
    print("Example: rate.py dissociation ATGCAGT -bootstrap")
    exit()


def main():
    if(len(sys.argv) < 2):
        errorPrint()

    type = sys.argv[1]
    mySequence = sys.argv[2]
    start_time = time.time()

    doBootstrap = False
    if(len(sys.argv) > 3):
        if(str(sys.argv[3])=="-bootstrap"):
            doBootstrap = True

    if type == "dissociation":
        from dissociation import computeAndWriteToCL as compute
    elif type == "hybridization":
        from anneal import computeAndWriteToCL as compute
    else:
        errorPrint()
    result = compute(mySequence, doBootstrap)
    print(f"Computing took {time.time() - start_time:.4f} s")


if __name__ == "__main__":
    main()
