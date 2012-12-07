#!/bin/sh

PROGRAM=$1
CONF=$2
SCALE=$3

echo $PROGRAM, $CONF, $SCALE

cd sommelier
make clean
make ALG=$PROGRAM

cd ../workspace
./wattchify ../confs/$CONF ../confs/tmp.conf
./cactify ../confs/tmp.conf ../confs/test.conf

rm sim$SCALE.result

./sesc.smp -c../confs/test.conf -dsim$SCALE -fresult  ../sommelier/sommelier.sesc -t $SCALE

../scripts/report.pl sim$SCALE.result