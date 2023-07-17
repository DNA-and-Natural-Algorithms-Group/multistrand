#!/usr/bin/env bash
#
# Reinstall the `multistrand` package with `pip3`, forcing full recompilation

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd $root
echo $root
echo
pip3 --version
pip3 uninstall multistrand
echo
rm -rfv ./build
echo
pip3 install -vvv -e .
