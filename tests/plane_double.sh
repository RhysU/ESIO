#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x plane_double ]; then
    echo "plane_double binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping plane_double"
    exit 0
fi

set -e # Fail on first error
for auxstride in ""                \
                 "--auxstride-b=3" \
                 "--auxstride-a=5"
do
    for cmd in "mpiexec -np 1 ./plane_double -p 11 -u 13" \
               "mpiexec -np 2 ./plane_double -p  5 -u  7" \
               "mpiexec -np 3 ./plane_double -p  7 -u  5"
    do
        for dir in "B" "A"
        do
            echo -n "Distribute $dir: "
            echo $cmd -d $dir $auxstride
            $cmd -d $dir $auxstride
        done
    done
done
