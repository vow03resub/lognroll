#!/usr/bin/python
#-*- coding: utf-8 -*-

import os
import sys
sys.path.append('../')
from logparser import Drain, IPLoM, AEL, LFA, LKE, LenMa, LogCluster, LogMine, LogSig, SHISO, SLCT, Spell
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
        parser.add_argument('--method',  type=str, required=True, help='Template discovery technique: Drain, IPLoM, AEL, LFA, LKE, LenMa, LogCluster, LogMine, LogSig, SHISO, SLCT, Spell, logmatch')
        parser.add_argument('--logfile',  type=str, required=False, help='Log file to use. Default is hadoop_clean_micro.log.')
        args = parser.parse_args()
        LTD_name = args.method # LTD: Log Template Discovery
        if args.logfile!=None:
            benchmark_settings['HDFS']['log_file'] = args.logfile

    except Exception as e:
        print('Error: %s' % str(e))

    method_list = ["AEL", "Drain", "IPLoM", "LFA", "LKE", "LenMa", "LogCluster", "LogMine", "LogSig", "SHISO", "SLCT", "Spell", "logmatch"]
    if LTD_name not in method_list:
        print("Unknown technique name. It must be one of the followings.")
        print(method_list)
        sys.exit(0)

    print("Log file:"+benchmark_settings['HDFS']['log_file'])

    for dataset, setting in benchmark_settings.items():

        if LTD_name=="Drain":
            parser = Drain.LogParser(log_format=setting['log_format'], indir=setting['input_dir'], outdir=setting['output_dir'], rex=setting['regex'], depth=setting['depth'], st=setting['st'], keep_para=False)
        elif LTD_name=="IPLoM":
            parser = IPLoM.LogParser(log_format=setting['log_format'], indir=setting['input_dir'], outdir=setting['output_dir'], CT=setting['CT'], lowerBound=setting['lowerBound'], rex=setting['regex'], keep_para=False)
        elif LTD_name=="AEL":
            parser = AEL.LogParser(setting['input_dir'], setting['output_dir'], log_format=setting['log_format'], rex=setting['regex'], minEventCount=setting['minEventCount'], merge_percent=setting['merge_percent'], keep_para=False)
        elif LTD_name == "LFA":
            parser = LFA.LogParser(setting['input_dir'], setting['output_dir'], log_format=setting['log_format'], rex=setting['regex'])
        elif LTD_name == "LKE":
            parser = LKE.LogParser(indir=setting['input_dir'], outdir=setting['output_dir'], log_format=setting['log_format'],rex=setting['regex'], split_threshold=setting['split_threshold'])
        elif LTD_name == "LenMa":
            parser = LenMa.LogParser(setting['input_dir'], setting['output_dir'], setting['log_format'],rex=setting['regex'])
        elif LTD_name == "LogCluster":
            parser = LogCluster.LogParser(setting['input_dir'],  setting['log_format'], setting['output_dir'], rsupport=setting['rsupport'])
        elif LTD_name == "LogMine":
            parser = LogMine.LogParser(setting['input_dir'], setting['output_dir'], setting['log_format'],rex=setting['regex'], levels=setting['levels'], max_dist=setting['max_dist'], k=setting['k'])
        elif LTD_name == "LogSig":
            parser = LogSig.LogParser(setting['input_dir'],  setting['output_dir'], setting['group_number'], setting['log_format'], rex=setting['regex'])
        elif LTD_name == "SHISO":
            parser = SHISO.LogParser(setting['log_format'], indir=setting['input_dir'],outdir=setting['output_dir'], rex=setting['regex'], maxChildNum=setting['maxChildNum'], mergeThreshold=setting['mergeThreshold'], formatLookupThreshold=setting['formatLookupThreshold'], superFormatThreshold=setting['superFormatThreshold'])
        elif LTD_name == "SLCT":
            parser = SLCT.LogParser(log_format=setting['log_format'], indir=setting['input_dir'], outdir=setting['output_dir'], support=setting['support'], rex=setting['regex'])
        elif LTD_name == "Spell":
            parser = Spell.LogParser(log_format=setting['log_format'], indir=setting['input_dir'], outdir=setting['output_dir'], tau=setting['tau'], rex=setting['regex'], keep_para=False)
        elif LTD_name == "logmatch": 
            log_filepath = setting['input_dir'] + setting['log_file']
            template_filepath = setting['output_dir'] + setting['log_file'] + '_templates.csv'
            matcher = regexmatch.PatternMatch(outdir=setting['output_dir'], n_workers=setting['n_workers'], logformat=setting['log_format'])
            matcher.match(log_filepath, template_filepath)

        if LTD_name!='logmatch':
            parser.parse(os.path.basename(setting['log_file']))
        
        bs=benchmark_settings['HDFS'] 
        target_file_1 = bs['output_dir']+bs['log_file']+'_structured.csv'
        target_file_2 = bs['output_dir']+bs['log_file']+'_templates.csv'
        new_name_1 = bs['output_dir']+os.path.splitext(bs['log_file'])[0]+'.'+ LTD_name +'.structured'
        new_name_2 = bs['output_dir']+os.path.splitext(bs['log_file'])[0]+'.'+ LTD_name +'.templates'
        os.rename(target_file_1, new_name_1)
        os.rename(target_file_2, new_name_2)
