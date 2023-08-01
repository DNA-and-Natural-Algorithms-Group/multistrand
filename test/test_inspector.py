# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from functools import partial
from types import FunctionType

import pytest

from multistrand.options import Literals
from multistrand.experiment import hybridization, dissociation, standardOptions
from multistrand.builder import Builder


class Test_Inspector:

    def create_config(self, reaction: FunctionType, seq: str, _):
        o = standardOptions()
        reaction(o, seq)
        o.simulation_mode = Literals.first_passage_time
        return o

    @pytest.mark.parametrize(
        "reaction, seq",
        [(hybridization, "CCC"), (dissociation, "CCCCCATTAAC")])
    def test_inspector(self, reaction: FunctionType, seq: str):
        b = Builder(partial(self.create_config, reaction, seq), [])
        b.genAndSavePathsFile(inspecting=True)
        print(b)
