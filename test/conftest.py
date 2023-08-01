# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

import sys
from pathlib import Path
from importlib import import_module
from types import ModuleType

import pytest

import multistrand


@pytest.fixture
def tutorials() -> ModuleType:
    """
    Access modules in the `tutorials/` folder.
    """
    sys.path.append(Path(multistrand.__file__).parents[2].resolve().as_posix())
    return import_module("tutorials")
