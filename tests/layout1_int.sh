#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x layout1_int ]; then
    echo "layout1_int binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping layout1_int"
    exit 0
fi

set -e # Fail on first error
for auxstride in "--auxstride-c=7 --auxstride-b=5 --auxstride-a=3"
do
    for cmd in "mpiexec -np 3 ./layout1_int -p  11 -u  13"
    do
        for dir in "C" "B" "A"
        do
                echo -n "Distribute $dir:"
                echo $cmd -d $dir $auxstride
                $cmd -d $dir $auxstride
        done
    done
done
