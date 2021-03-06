#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x line_float ]; then
    echo "line_float binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping line_float"
    exit 0
fi

set -e # Fail on first error
for auxstride in "--auxstride-a=7"
do
    for cmd in "mpiexec -np 3 ./line_float -p  11"
    do
        echo $cmd $auxstride
        $cmd $auxstride
    done
done
