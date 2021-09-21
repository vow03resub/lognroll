#!/bin/bash

LOGNAME=$1
COUNT=$2

if [ -z "$LOGNAME" ]; then
    echo "Please provide log file to use."
    echo "   Choice: hadoop, openstack, cassandra"
    exit 0
fi

if [ -z "$COUNT" ]; then
    echo "Please specify how many times to repeat each log size case."
	echo "One round will run all 20 sizes of all 4 logs once."
    exit 0
fi

TIMESTAMP=`date +%Y%m%d_%H%M%S`

for ((i=1;i<=$COUNT;i++));
do
    for ((j=1;j<=20;j++)); # 20 is the number of log sizes (from 1k to 20k)
    do
        echo Round $i Processing ${LOGNAME}_${j}k ...
        ./go.sh ${LOGNAME}_${j}k.log 1 >> ../output/__logpai.${LOGNAME}.${TIMESTAMP}.txt

#        echo Processing hadoop_${j}k ...
#        ./go.sh hadoop_${j}k.log 1 >> ../output/__logpai.hadoop_${j}k.${TIMESTAMP}.txt
#
#        echo Processing openstack_${j}k ...
#        ./go.sh openstack_${j}k.log 1 >> ../output/__logpai.openstack_${j}k.${TIMESTAMP}.txt
#
#        echo Processing cassandra_${j}k ...
#        ./go.sh cassandra_${j}k.log 1 >> ../output/__logpai.cassandra_${j}k.${TIMESTAMP}.txt
    done
done

