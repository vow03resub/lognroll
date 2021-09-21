#!/bin/bash

#Method_list="AEL LogCluster LFA IPLoM Spell Drain SHISO LenMa"
#Method_list="LogSig"

#LOGFILE_NAME="hadoop_clean.log"
#LOGFILE_NAME="cassandra_debug.log"
#LOGFILE_NAME="openstack_small_clean.log"

LOGNAME_LIST=$1
METHOD_LIST=$2

if [ -z "$1" ]; then
    echo "Please specify all the followings:"
    echo "   (1) List of logs to use. One or more of these - hadoop, openstack, and cassandra."
    echo "   (2) List of methods in double quotation. ex) \"Drain IPLoM AEL\""
	exit 0
fi

if [ -z "$2" ]; then
    echo "Please specify all the followings:"
    echo "   (1) List of logs to use. One or more of these - hadoop, openstack, and cassandra."
    echo "   (2) List of methods in double quotation. ex) \"Drain IPLoM AEL\""
	exit 0
fi

# Build a mapping between the log (type) names and the log file names.
declare -A LogfileMap=( ["hadoop"]="hadoop_clean.log" ["openstack"]="openstack_small_clean.log" ["cassandra"]="cassandra_debug.log" )

for Logname in $LOGNAME_LIST;
do

    Logfile_name=${LogfileMap[$Logname]}
	Logname_prefix=`echo $Logfile_name | sed 's/\.log//g'`

    for Methodname in $METHOD_LIST;
    do

        if [ "$Methodname" = "Lognroll" ]; then
            cd ..
            ./lognroll_SLCL.py --logname $Logname
            cd eval
        else
            echo -e "\033[0;46m${method}\033[0m"
            python ../analyze_others.py --logfile ../log_sample/$Logfile_name --template Result.final/${Logname_prefix}.${Methodname}.templates
    
        fi

    done
done

