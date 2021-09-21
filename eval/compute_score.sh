#!/bin/bash

#LOGFILE_NAME="hadoop_clean.log"
#LOGFILE_NAME="cassandra_debug.log"
LOGFILE_NAME="openstack_small_clean.log"

#METHOD_LIST="AEL LogCluster LFA IPLoM Spell Drain SHISO LenMa LogSig"
#METHOD_LIST="LenMa SHISO"
METHOD_LIST="LogCluster LFA IPLoM Spell SHISO LenMa LogSig"

TMPFILE="__VKED732KDJFWHF873__"

for logfile_name in $LOGFILE_NAME;
do
    for methodname in $METHOD_LIST;
    do
        echo -e "\033[1;95m<>"  $logfile_name $methodname "\033[0m"

        /usr/bin/time -o $TMPFILE -f "%e" ./run_method.py --method $methodname --logfile $logfile_name >/dev/null

        echo -e "    Elapsed time: \033[0;103m"`cat $TMPFILE`"\033[0m seconds"

        #prefix=`echo $logfile_name | sed 's/\.log$//g'`
        #template_file=Result.final/${prefix}.${methodname}.templates
        #../34_analyze_logpai.py --logfile ../log_sample/${logfile_name} --template $template_file

    done
done

rm -rf $TMPFILE
