# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from types import ModuleType, FunctionType
from typing import Callable, Dict
from functools import partial
from importlib import import_module

import pytest

from multistrand.system import SimSystem
from multistrand.options import Options


class Test_Tutorials:
    """
    Check that the small tutorials finish successfully.
    """

    @staticmethod
    def check_output(main: FunctionType, capfd: pytest.CaptureFixture) -> None:
        options, simsystem = main()
        assert isinstance(options, Options)
        assert isinstance(simsystem, SimSystem)

        capture = capfd.readouterr()
        print(capture.out)
        assert capture.out != ""
        assert capture.err == ""

    @staticmethod
    def build_tutorial_test(group: str, name: str):
        def test(self, tutorials: ModuleType, capfd: pytest.CaptureFixture):
            self.check_output(import_module(f"tutorials.{group}.{name}").main,
                              capfd)
        return test

    @classmethod
    def build_test_subclass(cls, name: str,
                            tests: Dict[str, Callable[[str], FunctionType]]):
        klass_name = f"{cls.__name__}_{name}"
        globals()[klass_name] = type(
            klass_name, (cls,),
            {f"test_{tut}": test(tut) for tut, test in tests.items()})

    @classmethod
    def build_all_tests(cls):
        """
        Create a test for each target tutorial, grouped into test classes by
        import directory.
        """
        cls.build_test_subclass("Misc", {
            tut: partial(Test_Tutorials.build_tutorial_test, "misc")
            for tut in ["sample_trace", "sample_trace_no_function",
                        "sample_trace_arrhenius", "debug_mac"]})


Test_Tutorials.build_all_tests()
