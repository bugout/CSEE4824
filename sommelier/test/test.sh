#!/bin/bash

set -e

SOMMELIER=../sommelier
MIDFILE=all.mat

for m in 1 2 5 30 51; do
    for n in 3; do
	for s in 10 55; do
	    $SOMMELIER -m $m -n $n -s $s -d > $MIDFILE
	    ./test.pl $MIDFILE
	done
    done
done

echo "test OK"
