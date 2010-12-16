#!/bin/bash
# Unit tests for the restart_rename routine in restart_rename.c
# Much, much easier to test these from the shell than from C
set -e # Fail on first error

CMD=`readlink -f restart_rename`
if ! [ -x "$CMD" ]; then
    echo "restart_rename binary not found or not executable"
    exit 1
fi

# Create and change to a temporary directory
# Automatically clean up this temporary directory on exit
TMPTMPL=test.`basename $0`.XXXXXXXXXX
TMPDIR=`mktemp -d --tmpdir $TMPTMPL`
trap 'rm -rf "$TMPDIR"' 0
trap ' exit ' 1 2 15
cd "$TMPDIR"


# Note that many of the filename variations are tested # within the
# restart_helpers test case.  Here we focus on basic operation only.

echo $0 running in `pwd`

echo -n Simple test with a single hash sign
echo b > test0
echo c > test1
echo e > test4
echo junk1 > junk
$CMD junk 'test#' 3
test `cat test0` == 'junk1'
test `cat test1` == 'b'
test `cat test2` == 'c'
test `cat test4` == 'e'
echo junk2 > junk
$CMD junk 'test#' 3
test `cat test0` == 'junk2'
test `cat test1` == 'junk1'
test `cat test2` == 'b'
test `cat test4` == 'e'
echo junk3 > junk
$CMD junk 'test#' 3
test `cat test0` == 'junk3'
test `cat test1` == 'junk2'
test `cat test2` == 'junk1'
test `cat test4` == 'e'
rm test*
echo " OK"

echo -n Simple test with multiple hash signs
echo b > test0
echo c > test1
echo e > test4
echo junk1 > junk
$CMD junk 'test##' 3
test `cat test00` == 'junk1'
test `cat test01` == 'b'
test `cat test02` == 'c'
test `cat test4` == 'e'
echo junk2 > junk
$CMD junk 'test##' 3
test `cat test00` == 'junk2'
test `cat test01` == 'junk1'
test `cat test02` == 'b'
test `cat test4` == 'e'
echo junk3 > junk
$CMD junk 'test##' 3
test `cat test00` == 'junk3'
test `cat test01` == 'junk2'
test `cat test02` == 'junk1'
test `cat test4` == 'e'
rm test*
echo " OK"

echo -n Simple test where keep_howmany sets number of digits retained
echo b > test0
echo c > test1
echo junk1 > junk
$CMD junk 'test#' 100
test `cat test00` == 'junk1'
test `cat test01` == 'b'
test `cat test02` == 'c'
echo junk2 > junk
$CMD junk 'test#' 1000
test `cat test000` == 'junk2'
test `cat test001` == 'junk1'
test `cat test002` == 'b'
test `cat test003` == 'c'
echo junk3 > junk
$CMD junk 'test#' 10000
test `cat test0000` == 'junk3'
test `cat test0001` == 'junk2'
test `cat test0002` == 'junk1'
test `cat test0003` == 'b'
test `cat test0004` == 'c'
rm test*
echo " OK"
