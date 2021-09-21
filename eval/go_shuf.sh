#!/bin/bash

LOGNAME=$1
COUNT=$2
K_START=$3
K_END=$4
METHOD_LIST=$5
#METHOD_LIST="AEL LogCluster LFA IPLoM Spell Drain SHISO LenMa LogSig LogMine SLCT MoLFI LKE"

if [ ! $# -eq 5 ]; then
    echo "Please specify all the followings:"
	echo "   (1) Log to use. One of these - hadoop, openstack, and cassandra."
	echo "   (2) Repeat count"
	echo "   (3) Start size of a log file in 1000 line unit."
	echo "   (4) End size of a log file in 1000 line unit."
	echo "   (5) List of methods in double quotation. ex) \"Drain IPLoM AEL\""
    exit 0
fi

# Build a mapping between the log (type) names and the log file names.
declare -A LogfileMap=( ["hadoop"]="hadoop_clean.log" ["openstack"]="openstack_small_clean.log" ["cassandra"]="cassandra_debug.log" )
LOGFILE=${LogfileMap[$LOGNAME]}

TimeLimit=900 # Each method will be allowed to run only 15 minutes.
TempFile="__GOSH380COWX48392__"
ShufLog="__K83KDFOX3034FO22__"
TimeStamp=`date +%Y%m%d_%H%M%S`


for ((i=1;i<=$COUNT;i++));
do
    echo -e "\033[93;42mCount:${i}\033[0m">&2
    for Method in $METHOD_LIST;
	do

            for ((j=${K_START};j<=${K_END};j++));
            do

                #rm -rf Result.HDFS

                # Create new log file using specified size each time we run it
                shuf -n ${j}000 ../log_sample/$LOGFILE > ../log_sample/${ShufLog}

                echo -e "\033[0;34mUsing "$j"k lines of "$LOGFILE for $Method"\033[0m">&2
        
            if [ "$Method" = "MoLFI" ]; then # MoLFI must use python3. So, we run it separately.
                echo "Method_name:" $Method `date +%H:%M:%S` >&2
                /usr/bin/time -o $TempFile -f "%e" timeout $TimeLimit ./run_MoLFI.py --method $Method --logfile $ShufLog >/dev/null
				echo "Elapsed time:" `cat $TempFile`
                echo DONE $Method ${LOGNAME}_${j}k $i _ _ 0:`cat $TempFile` >> ../output/__logpai.${LOGNAME}.${TimeStamp}.txt
                rm -rf $TempFile
    
            else
                echo "Method_name:" $Method `date +%H:%M:%S` >&2
                /usr/bin/time -o $TempFile -f "%e" timeout $TimeLimit ./run_method.py --method $Method --logfile $ShufLog >/dev/null
				echo "Elapsed time:" `cat $TempFile`
                echo DONE $Method ${LOGNAME}_${j}k $i _ _ 0:`cat $TempFile` >> ../output/__logpai.${LOGNAME}.${TimeStamp}.txt
                rm -rf $TempFile
            fi

            rm -rf ../log_sample/${ShufLog}

        done
    done
done

rm -rf $TempFile
# Delete all pyc files
cd ..
find . -name *.pyc -exec rm {} \;
