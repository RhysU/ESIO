#!/bin/bash
# Must use mpiexec to run serial-but-MPI-enabled tests on some MPI stacks.  In
# particular, mvapich seems to exhibit this problem.  Moreover, we cannot
# always use mpiexec on some login nodes.  Better to warn the user that a test
# was skipped then worry them when make check fails as a result.

if ! [ -x plane_int ]; then
    echo "plane_int binary not found or not executable"
    exit 1
fi

if ! [ -x plane_int_f ]; then
    echo "plane_int_f binary not found or not executable"
    exit 1
fi

if ! which mpiexec > /dev/null ; then
    echo "WARNING: Unable to find mpiexec; skipping plane_int"
    exit 0
fi

set -e # Fail on first error
for auxstride in "--auxstride-b=5 --auxstride-a=3"
do
    for cmd in "mpiexec -np 3 ./plane_int -p  11 -u  13"
    do
        for dir in "B" "A"
        do
                echo -n "Distribute $dir:"
                echo $cmd -d $dir $auxstride
                $cmd -d $dir $auxstride
        done
    done
done
for cmd in "mpiexec -np 1 ./plane_int_f" \
           "mpiexec -np 2 ./plane_int_f"
do
    echo $cmd
    $cmd
done
