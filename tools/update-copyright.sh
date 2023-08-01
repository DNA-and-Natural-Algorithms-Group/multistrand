#! /usr/bin/env bash

# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

# default configuration
first_year=2008
old_last_year=""
new_last_year=""

# evaluate arguments
if !(( $# == 1 || $# == 2 )); then
    echo "Expecting 1 or 2 arguments:"
    echo "   old_last_year [new_last_year]"
    echo "   (int)         (int)"
    exit 1
fi;
old_last_year=$(($1))
if (( $# == 2 )); then
    new_last_year=$(($2)); else
    new_last_year=$(($old_last_year + 1))
fi;

# determine absolute path of repository
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd $repo

# update copyright notices
grep -rlIE --exclude-dir=".git" --exclude-dir=".tox" --exclude-dir="build" \
    "$first_year-+$old_last_year" \
    | xargs sed -i"" -E -e "s:($first_year-+)$old_last_year:\1$new_last_year:g"
