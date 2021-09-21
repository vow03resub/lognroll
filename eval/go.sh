#!/bin/bash

LOGFILE=$1
COUNT=$2

# Method will be allowed to run only 10 minutes.
CUTOFF=600

if [ -z "$LOGFILE" ]; then
    echo "Please designate a log file to use."
    exit 0
fi

if [ -z "$COUNT" ]; then
    echo "Please specify how many times to repeat."
    exit 0
fi

TMPFILE="__GOSH38048392__"

for ((i=1;i<=$COUNT;i++));
do

    echo "Round <"$i"> for " $LOGFILE >&2

    #METHOD_LIST="Drain IPLoM AEL LFA LKE LenMa LogCluster LogMine LogSig SHISO SLCT Spell"
    METHOD_LIST="AEL LogCluster LFA IPLoM Spell Drain SHISO LenMa LogSig"
    for METHOD in $METHOD_LIST; do
        echo "**** Method_name:" $METHOD `date +%H:%M:%S` >&2
        /usr/bin/time -o $TMPFILE timeout $CUTOFF ./run_method.py --method $METHOD --logfile $LOGFILE
        echo DONE $METHOD $LOGFILE $i `cat $TMPFILE`
        rm -rf $TMPFILE
    done

    ## MoLFI must use python3. So, we run it separately.
    #METHOD_LIST="MoLFI"
    #for METHOD in $METHOD_LIST; do
    #    echo "**** Method_name:" $METHOD `date +%H:%M:%S` >&2
    #    /usr/bin/time -o $TMPFILE timeout $CUTOFF ./run_MoLFI.py --method $METHOD --logfile $LOGFILE
    #    echo DONE $METHOD $LOGFILE $i `cat $TMPFILE`
    #    rm -rf $TMPFILE
    #done

done

# Delete all pyc files
cd ..
find . -name *.pyc -exec rm {} \;

