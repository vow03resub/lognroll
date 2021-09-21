#!/bin/bash

LOGNAME_LIST=$1
METHOD_LIST=$2

if [ -z "$LOGNAME_LIST" ]; then
    echo "Please specify all the followings:"
    echo "   (1) List of logs to use. One or more of these - hadoop, openstack, and cassandra."
    echo "   (2) List of methods in double quotation. ex) \"Drain IPLoM AEL\""
    exit 0
fi

if [ -z "$METHOD_LIST" ]; then
    echo "Please specify all the followings:"
    echo "   (1) List of logs to use. One or more of these - hadoop, openstack, and cassandra."
    echo "   (2) List of methods in double quotation. ex) \"Drain IPLoM AEL\""
    exit 0
fi

#LOGFILE_NAME="hadoop_clean.log"
#LOGFILE_NAME="cassandra_debug.log"
#LOGFILE_NAME="openstack_small_clean.log"
#METHOD_LIST="AEL LogCluster LFA IPLoM Spell Drain SHISO LenMa LogSig"

TimeStamp=`date +%Y%m%d_%H%M%S`
TmpPrefix="__VKED732KDJFW__"
TMPFILE=`echo ${TmpPrefix}${RANDOM}${TimeStamp}`

# Build a mapping between the log (type) names and the log file names.
declare -A LogfileMap=( ["hadoop"]="hadoop_clean.log" ["openstack"]="openstack_small_clean.log" ["cassandra"]="cassandra_debug.log" )

for Logname in $LOGNAME_LIST;
do

	LOGFILE=${LogfileMap[$Logname]}

    for Methodname in $METHOD_LIST;
    do
        echo -e "\033[2;31m<>"  $Logname $LOGFILE $Methodname "\033[0m"

		if [ "$Methodname" = "MoLFI" ]; then
			/usr/bin/time -o $TMPFILE -f "%e" ./run_MoLFI.py --method $Methodname --logfile $LOGFILE >/dev/null
        	echo -e "    Elapsed time: \033[0;103m"`cat $TMPFILE`"\033[0m seconds"
		else
        	/usr/bin/time -o $TMPFILE -f "%e" ./run_method.py --method $Methodname --logfile $LOGFILE >/dev/null
        	echo -e "    Elapsed time: \033[0;103m"`cat $TMPFILE`"\033[0m seconds"
		fi

    done
done

rm -rf $TMPFILE
