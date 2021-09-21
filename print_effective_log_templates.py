#!/usr/bin/python
#-*- coding:utf-8 -*-
import re
import os
import csv
import sys
import time
import pickle
import argparse
import logutils

STAR_THRESHOLD=25

seqnum=1

MOD_FACTOR = 32

debug_mode=0

def mark_matched_logs(logs, mask, template):
    marked=0
    for i in range(0, len(logs)):
        if mask[i]>-1:
            continue

        #if i%(len(logs)/30)==0:
        #    sys.stdout.write('\r'+"    \033[1;91m "+"{0:.1f}".format(float(i*100)/float(len(logs)))+"% \033[0m"+template)
        #    sys.stdout.flush()

        log = logs[i]

        matched = re.match('^'+template+'$', log)
        if matched!=None:
            mask[i]=99
            marked+=1
    return marked


# Read all logs from multiple files.
def read_log_files(flist):
    logs = []
    for i in range(0,len(flist)):
        for log in flist[i]:
            if len(log.strip())==0:
                continue
            log = " ".join(log.split())
            logs.append(log)
    return logs


def load_csv_template_file(f):
    template_csv = csv.reader(open(f,'r'))
    template_list = []
    for t in template_csv:
        template_list.append(t)
    template_list.pop(0) 

    #re_pattern=['.', '*', '+', '?', '^', '$', '[', ']', '{', '}', '(', ')', '|']
    re_pattern=['*', '+', '?', '^', '$', '[', ']', '{', '}', '(', ')', '|']
    for t in template_list:
        if 'SHISO' in args.template:
            t[1] = t[1].replace('*', '~STAR~')
        #elif 'LFA' in args.template:
        #    t[1] = t[1].replace('*', '~STAR~')
        elif 'LogCluster' in args.template:
            #t[1] = t[1].replace('*', '~STAR~')
            sss = re.sub("\*{\d+,\d+}","~STAR~",t[1])
            t[1] = sss
        elif 'MoLFI' in args.template:
            t[1] = t[1].replace('#spec#', '~STAR~')
            t[1] = t[1].replace(' :', ':') 
            t[1] = t[1].replace(' ;', ';') 
            t[1] = t[1].replace(' ,', ',') 
            t[1] = t[1].replace('< ', '<')
            t[1] = t[1].replace(' >', '>')
            t[1] = t[1].replace('( ', '(')
            t[1] = t[1].replace(' )', ')')
            t[1] = t[1].replace('*', '~STAR~')
        else:
            t[1] = t[1].replace('<*>', '~STAR~')
        
        for c in re_pattern:
            t[1]=t[1].replace(c, '\\'+c)

        #t[1] = t[1].replace('~STAR~', '.*')
        t[1] = t[1].replace('~STAR~', '\S*')
        t[1] = ' '.join(t[1].split())  
    return template_list

if __name__ == '__main__':
    
    try:
        parser = argparse.ArgumentParser(description="")
        parser.add_argument('--logfile', type=argparse.FileType('r'), nargs='+', required=True, help='Name of the log files')
        parser.add_argument('--template', required=True, help='Name of the template file')

        args=parser.parse_args()
        logfile_list=args.logfile
        template_file=args.template

    except Exception, e:
        print('Error: %s' % str(e))

    print('    Input log file: %s' %(logfile_list[0].name))
    #print('Input template file: %s' %(template_file))
    raw_logs = read_log_files(logfile_list)
    print('    Total number of logs: ' + str(len(raw_logs)))
    loaded_templates = load_csv_template_file(template_file)
    print "    Total number of templates loaded:",len(loaded_templates)

    selected = []

    log_templates=[] 
    log_mask=[-1 for _ in range(len(raw_logs))] 
    for i in range(0,len(loaded_templates)):
        t = loaded_templates[i]

        marked = mark_matched_logs(raw_logs, log_mask, t[1])
#            if t[1].count("*")==STAR_THRESHOLD:
#                print "    Finished matching logs.", marked, "logs matched."
        #print "["+str(i)+"]",marked,t[1]

        log_templates.append((marked, t[1])) 
        log_mask=[-1 for _ in range(len(raw_logs))] 

        progress_mod=max(5,len(loaded_templates)/40)
        if i%progress_mod==0:
            sys.stdout.write('\r'+"    \033[0;103mMatch counting processed "+"{0:.1f}".format(float(i*100)/float(len(loaded_templates)))+"% \033[0m")
            sys.stdout.flush()
    sys.stdout.write("\n")
 
    sum_matched=0
    for t in sorted(log_templates,reverse=True):
        #print t[0], t[1]
        sum_matched += int(t[0])
    print "    Sum of matched logs:",sum_matched

    for t in sorted(log_templates, reverse=True):
        # first, build a list of index to delete
        to_delete = []
        for i in range(0,len(raw_logs)):
            log = raw_logs[i]
            matched = re.match("^"+t[1]+"$",log)

            if matched!=None:
                to_delete.append(i)

        # delete matched logs
        before_removal = len(raw_logs)
        to_delete = sorted(to_delete)
        for i in reversed(sorted(to_delete)):
            del raw_logs[i]
        del_count = before_removal - len(raw_logs)
        #print "\033[0;32mRemoved",del_count,"logs.\033[0m",t[1]
        #if del_count>0:
        #    print "\033[33;31m",format(del_count,'5d'), format(int(t[0]),'5d'), t[1], "\033[0m"

        if del_count>0:
            #selected.append(t[1])
            selected.append({"count":del_count,"template":str(t[1])})
        #if len(selected)==2:
        #    break
 
    print "   ",len(raw_logs),"logs remaining."
    print "    Initial template count:", len(log_templates)
    print "    Selected template count:", len(selected)

    mycnt=1
    for x in selected:
        print mycnt, x['count'], x['template']
        mycnt+=1

    sys.exit(0)

    mylogs = read_log_files(logfile_list)
    SL,CPL = logutils.compute_slcpl(mylogs, selected)

    print "    SL= "+str(SL)
    print "    CPL= "+str(CPL)

    log_name = template_file.replace("Result.final/","").replace(".templates","").split('.')[0].split('_')[0]
    if log_name=="hadoop":
        log_name_color= "91;106m"
    elif log_name=="openstack":
        log_name_color= "91;103m"
    elif log_name=="cassandra":
        log_name_color= "1;106m"
    else:
        log_name_color= "32;31m"
    method_name = template_file.replace("Result.final/","").replace(".templates","").split('.')[1]
    if method_name=="Lognroll":
        mname_color= "1;94m"
    elif method_name=="AEL":
        mname_color= "1;95m"
    elif method_name=="Drain":
        mname_color= "1;91m"
    elif method_name=="SHISO":
        mname_color= "1;32m"
    elif method_name=="IPLoM":
        mname_color= "1;106m"
    elif method_name=="LFA":
        mname_color= "1;103m"
    elif method_name=="LenMa":
        mname_color= "45;33m"
    elif method_name=="LogSig":
        mname_color= "32;33m"
    elif method_name=="LogMine":
        mname_color= "33;31m"
    elif method_name=="SLCT":
        mname_color= "43;31m"
    elif method_name=="MolFI":
        mname_color= "91;106m"
    elif method_name=="LKE":
        mname_color= "1;98m"
    else:
        mname_color= "1;103m"

    print "    \033[1;94mscore= "+str(SL-CPL)+"\033[0m", "\033["+log_name_color+log_name+"\033[0m", "\033["+mname_color+method_name+"\033[0m"
    print sum_matched, len(raw_logs), len(log_templates), len(selected), SL, CPL
