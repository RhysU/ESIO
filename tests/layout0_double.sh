#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x layout0_double ]; then
    echo "layout0_double binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping layout0_double"
    exit 0
fi

set -e # Fail on first error
for auxstride in ""                \
                 "--auxstride-c=3" \
                 "--auxstride-b=5" \
                 "--auxstride-a=7"
do
    for cmd in "mpiexec -np 1 ./layout0_double -p 11 -u 13" \
               "mpiexec -np 2 ./layout0_double -p  5 -u  7" \
               "mpiexec -np 3 ./layout0_double -p  7 -u  5"
    do
        for dir in "C" "B" "A"
        do
            echo -n "Distribute $dir: "
            echo $cmd -d $dir $auxstride
            $cmd -d $dir $auxstride
        done
    done
done
