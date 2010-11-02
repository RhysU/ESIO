#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x attribute_double ]; then
    echo "attribute_double binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping attribute_double"
    exit 0
fi

set -e # Fail on first error
for cmd in "mpiexec -np 1 ./attribute_double -n 11" \
           "mpiexec -np 2 ./attribute_double -n  5" \
           "mpiexec -np 3 ./attribute_double -n  2"
do
    echo $cmd
    $cmd
done
