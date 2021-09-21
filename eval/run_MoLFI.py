#!/usr/bin/python3
#-*- coding: utf-8 -*-

import os
import sys
sys.path.append('../')
from logparser import MoLFI
import pandas as pd
import argparse
from logparser.logmatch import regexmatch
#MolFI


benchmark_settings = {
    'HDFS': {
        'input_dir': '../log_sample/',        
        'output_dir': 'Result.HDFS/',         
        'log_file': 'hadoop_clean_micro.log', 
        'log_format': '<Content>',            

        'regex': [],                    # Drain, IPLoM
        'st': 0.5,                      # Drain
        'depth': 4,                     # Drain

        'CT': 0.35,         # IPLoM
        'lowerBound': 0.25, # IPLoM

        'minEventCount' : 2,   # AEL
        'merge_percent' : 0.5, # AEL

        'split_threshold' : 3, # LKE

        'rsupport'   : 10, # LogCluster

        'levels'     : 2,     # LogMine
        'max_dist'   : 0.001, # LogMine

        'k'          : 1, # LogMine

        'group_number' : 14, # LogSig

        'maxChildNum' : 4,              # SHISO
        'mergeThreshold' : 0.1,         # SHISO
        'formatLookupThreshold' : 0.3,  # SHISO
        'superFormatThreshold'  : 0.85, # SHISO

        'support' : 10, #SLCT

        'tau' : 0.5, #Spell

        'n_workers' :1 #logmatch
    },
}


if __name__ == '__main__':

    try:
        parser = argparse.ArgumentParser(description="")
        parser.add_argument('--method',  type=str, required=True, help='Template discovery technique: MoLFI')
        parser.add_argument('--logfile',  type=str, required=False, help='Log file to use. Default is hadoop_clean_micro.log.')
        args = parser.parse_args()
        LTD_name = args.method # LTD: Log Template Discovery
        if args.logfile!=None:
            benchmark_settings['HDFS']['log_file'] = args.logfile

    except Exception as e:
        print('Error: %s' % str(e))

    method_list = ["MoLFI"]
    if LTD_name not in method_list:
        print("Unknown technique name. It must be MoLFI.")
        print(method_list)
        sys.exit(0)

    print("Log file:"+benchmark_settings['HDFS']['log_file'])

    for dataset, setting in benchmark_settings.items():

        if LTD_name == "MoLFI":
            parser = MoLFI.LogParser(setting['input_dir'], setting['output_dir'], setting['log_format'],rex=setting['regex'])

        parser.parse(os.path.basename(setting['log_file']))
        
        bs=benchmark_settings['HDFS'] 
        target_file_1 = bs['output_dir']+bs['log_file']+'_structured.csv'
        target_file_2 = bs['output_dir']+bs['log_file']+'_templates.csv'
        new_name_1 = bs['output_dir']+os.path.splitext(bs['log_file'])[0]+'.'+ LTD_name +'.structured'
        new_name_2 = bs['output_dir']+os.path.splitext(bs['log_file'])[0]+'.'+ LTD_name +'.templates'
        os.rename(target_file_1, new_name_1)
        os.rename(target_file_2, new_name_2)
