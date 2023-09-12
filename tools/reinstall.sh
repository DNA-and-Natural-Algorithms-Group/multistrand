#!/usr/bin/env bash

# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

## Reinstall the `multistrand` package with `pip3`, forcing full recompilation.

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd $root
echo $root
echo
pip3 --version
pip3 uninstall multistrand
echo
rm -rfv ./build
find . -type f -regex ".*\.so" -print | xargs rm -v
echo
pip3 install -vvv -e .
