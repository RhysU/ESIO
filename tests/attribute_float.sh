#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x attribute_float ]; then
    echo "attribute_float binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping attribute_float"
    exit 0
fi

set -e # Fail on first error
for cmd in "mpiexec -np 1 ./attribute_float -n 8" \
           "mpiexec -np 2 ./attribute_float -n 4"
do
    echo $cmd
    $cmd
done
