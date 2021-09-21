#-*- coding:utf-8 -*-

import re
import os
import csv
import sys
import time
import math

STAR_THRESHOLD=25

seqnum=1
MOD_FACTOR = 32


debug_mode=0

delimiter_list = [ ' ', ',', ';', '|', ':', '<', '>', '=', '/','@'] # < and > are added because of '->'. By adding > as delimiter it will be broken up into '-' and the right side string.
keyval_pattern = { "pattern":"[^\s~]+", "label": "~300~" }


def split_by_delimiter(lvl, str_data, delim):
    global seqnum

    #if debug_mode:
    #    print " "*lvl+"\033[0;46m"+"Entering split_by_delimiter"+"\033[0m", "->"+str_data+"<-", "->"+delim+"<-"

    # if input string is just one delimiter character, just return it
    if str_data==delim:
        #if debug_mode:
        #    print " "*lvl+"\033[0;47m"+"Leaving split_by_delimiter"+"\033[0m", "->"+str_data+"<-", "->"+delim+"<-"
        return [delim]

    # TODO: I may need to convert all known patterns first before splitting by delimters.
    #s = detect_all_patterns(lvl+4,s)

    # If I detect key-value pattern, convert the value part into *
    if delim=='=' and '=' in str_data:
        s = ' '+str_data+' '
        # pad front and back with space
        found = True
        while found:
            diff = 0
            found = False
            matched = re.finditer("[\t ](\S+[=])"+keyval_pattern["pattern"]+"[\t ]", s)
            for m in matched:
                label_str = "~KV"+format(seqnum,'09d')+"~"
                seqnum += 1
                seqnum = seqnum % MOD_FACTOR 
                s = s[0:m.start()+1-diff]+m.group(1)+label_str+s[(m.end()-1)-diff:]
                diff = diff + ((m.end()-1)-(m.start()+1)) - len(m.group(1)+label_str)
                found = True
        #print "::::::::::->"+str_data+"<----->"+s+"<-"
        str_data = s[1:len(s)-1] # remove spaces at the front and back

    # split the input string by delimiter
    tokenized = re.split("("+delim+")",str_data) # keep the delimiter within the list

    # remove empty token
    removed = True
    while removed:
        removed = False
        if "" in tokenized:
            tokenized.remove("")
            removed = True
    #if debug_mode:
    #    print " "*lvl+"Removed empty token"+"\033[93;44m"+str(tokenized)+"\033[0m"

    # process each token
    for i in reversed(range(0,len(tokenized))):

        tok = tokenized[i]
        if tok==delim:
            continue

        if debug_mode and len(tok)==0:
            print " "*lvl+"\033[0;31m"+"WARNING 742: zero length token detected -"+"\033[0m", str_data
            print " "*(lvl+4), tokenized
            print log
            sys.exit(0)
            continue

        # if any delimiter in lower priority than current one exists, recursively call the split_by_delimiter
        low_delim_found = False
        cur = delimiter_list.index(delim)
        for j in range(cur+1,len(delimiter_list)):
            if delimiter_list[j] in tok:

                #if debug_mode:
                #    print " "*lvl+"\033[0;42m"+"[split_by_delimiter("+delim+")] CALLING itself("+delimiter_list[j]+")"+"\033[0m ->"+tokenized[i]+"<-"
                #    print " "*lvl, "****",tok, i
                #    print " "*lvl, "****",tokenized
                #    print " "*lvl, "****",tokenized[0:i]
                #    print " "*lvl, "****",tokenized[i+1:]

                tokenized = tokenized[:i] + split_by_delimiter(lvl+4, tok, delimiter_list[j]) + tokenized[i+1:]

                #if debug_mode:
                #    print " "*lvl+"\033[0;42m"+"[split_by_delimiter("+delim+")] Returned"+"\033[0m ->"+str(tokenized)+"<-"
                #    print " "*lvl, "****",tokenized

                low_delim_found = True
                break
        #if not low_delim_found:
        #    if debug_mode:
        #        print " "*lvl+"\033[0;33m"+"[split_by_delimiter("+delim+")] CALLING detect_all_patterns()"+"\033[0m ->"+tokenized[i]+"<-"
        #    tokenized[i] = detect_all_patterns(lvl+4,tokenized[i])

    #if debug_mode:
    #    print " "*lvl+"\033[0;47m"+"Leaving split_by_delimiter"+"\033[0m", tokenized
    return tokenized


def get_bracket_char(log):
    pos = len(log)
    bkt_open = None
    bkt_close = None
    # quote and double-quote are treated like parentheses, but since there is no left and right, it is handled differently in the if statement in Custom_split.
    bkt_pair = {"(":")","{":"}","[":"]","<":">","'":"'","\"":"\""}

    for c in "({[<'\"":
        loc = log.find(c)
        if loc>=0 and loc<pos:
            pos = loc
            bkt_open = c
            bkt_close = bkt_pair[c]
    return bkt_open, bkt_close, pos


def custom_split(log):

    tokenized = []
    # Determine the first occuring parentheses
    bracket_open, bracket_close, pos = get_bracket_char(log)
    if bracket_open==None: # no parentheses found
        return split_by_delimiter(8,log,' ')

    # Extract the content within the brackets pair and call recursively
    pstack = [] # for storing index
    qstack = [] # for storing char in case of ' or "
    for i,c in enumerate(log):
        if (c==bracket_open and c not in ["'","\""]) or (c=="'" and "'" not in qstack) or (c=="\"" and "\"" not in qstack):
            pstack.append(i)
            if c in ["'","\""]:
                qstack.append(c)
        elif (c==bracket_close and c not in ["'","\""]) or (c=="'" and "'" in qstack) or (c=="\"" and "\"" in qstack):
            if len(pstack)>0:
                spos = pstack.pop() # start and end position of parentheses segment
                epos = i
            # There are cases where > exists without opening <.
            # Ex) ...org.apache.hadoop.mapreduce.v2.app.MRAppMaster 1><LOG_DIR>/stdout 2><LOG_DIR>/stderr
            else:
                continue
            if len(qstack)>0 and c in ["'","\""]:
                schar = qstack.pop()
                if c!=schar:
                    print "ERROR 531", schar
                    print log
                    sys.exit(0)
            if len(pstack)==0: # if stack becomes empty
                middle_part = custom_split(log[spos+1:epos])
                ending_part = custom_split(log[epos+1:])

                #if debug_mode:
                #    print " "*4+"Whole:",log
                #    print " "*4+"spos=:",spos
                #    print " "*4+"epos=:",epos
                #    print " "*4+"Begin:",log[0:spos]
                #    print " "*4+"Middle:",middle_part
                #    print " "*4+"Ending:",ending_part


                tokenized = split_by_delimiter(8,log[0:spos],' ') + [bracket_open] + middle_part + [bracket_close] + ending_part

                #tokenized = split_by_delimiter(8,log[0:spos],' ') + [bracket_open]
                #if len(middle_part)>0:
                #    tokenized = tokenized + middle_part
                #tokenized = tokenized + [bracket_close]
                #if len(ending_part)>0:
                #    tokenized = tokenized + ending_part
                break

    # Unclosed parentheses may exist. In such case, there are unprocessed strings.
    # Ex) ucState = COMMITTED, replication# = 0 < minimum = 0
    if len(pstack)>0:
        spos = pstack[0] # take the first parentheses in the pstack and do recursive call
        epos = len(log) # assume that there is imaginary closing parenthesis
        middle_part = custom_split(log[spos+1:epos])
        tokenized = split_by_delimiter(4,log[0:spos],' ') + [bracket_open] + middle_part

    #if debug_mode:
    #    print "\033[0;35mcustom_split() returning: "+str(tokenized)+"\033[0m"
    return tokenized


def do_tokenization(logs):
    tlogs = []
    for i in range(0,len(logs)):
        log = logs[i]
        tlogs.append(custom_split(log))
    return tlogs

common_patterns = [

    {   "pattern": "hdfs://([a-zA-Z0-9_\-\.\*:\+]+/)+([a-zA-Z0-9_\-\.\*:\+]*(\?[a-zA-Z0-9_\-:\+]+=[a-zA-Z0-9_\-:\+]+(&[a-zA-Z0-9_\-:\+]+=[a-zA-Z0-9_\-:\+]+)*)?)",
        "serial": "1",
        "prefix":"hdfs_url" },

    {   "pattern": "hdfs://([a-zA-Z0-9_\-\.]+):([0-9]+)",
        "serial": "1",
        "prefix":"hdfs_url" },

    {   "pattern": "http://([a-zA-Z0-9_\-\.\*:\+]+/)+([a-zA-Z0-9_\-\.\*:\+]*(\?[a-zA-Z0-9_\-:\+]+=[a-zA-Z0-9_\-:\+]+(&[a-zA-Z0-9_\-:\+]+=[a-zA-Z0-9_\-:\+]+)*)?)",
        "serial": "1",
        "prefix":"http_url" },

    {   "pattern": "http://([a-zA-Z0-9_\-\.]+):([0-9]+)",
        "serial": "1",
        "prefix":"http_url" },

    {   "pattern": "https://([a-zA-Z0-9_\-\.\*:\+]+/)+([a-zA-Z0-9_\-\.\*:\+]*(\?[a-zA-Z0-9_\-:\+]+=[a-zA-Z0-9_\-:\+]+(&[a-zA-Z0-9_\-:\+]+=[a-zA-Z0-9_\-:\+]+)*)?)",
        "serial": "1",
        "prefix":"https_url" },

    {   "pattern": "https://([a-zA-Z0-9_\-\.]+):([0-9]+)",
        "serial": "1",
        "prefix":"https_url" },

    {   "pattern": "http://([a-zA-Z0-9_\-\.]+):([0-9]+)(/[a-zA-Z0-9_\-\.\*\+\-]+)+",
        "serial": "1",
        "prefix":"http_url" },

    {   "pattern": "/(var|tmp|home|usr|home|etc|opt|gogo|airwordcount|wikimean|wikimedian|wikistandarddeviation|cluster|jobhistory|node|ws)(/[a-zA-Z0-9_\-\.\*\+\-]+)+",
        "serial": "1",
        "prefix":"file_path" },

    {   "pattern": "\d+\.\d+\.\d+\.\d+(:\d+)?", # IP address and port number
        "serial": "1",
        "prefix":"ipaddr_port" },

#    {   "pattern":"\-?\d+\.\d+ (KB|GB|MB)", # 21.5 MB
#        "label": "~324~" },
#    {   "pattern":"\-?\d+ (KB|GB|MB)", # 5 GB
#        "label": "~325~" },
]



# Longest common subsequence
def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])
    return lcs_set


def recognize_known_patterns(log):
    new_log=log
    for item in common_patterns:
        found = True
        while found:
            diff = 0
            found = False
            matched = re.finditer(item["pattern"], new_log)
            for m in matched:
                label_str = item["prefix"]+"_"+format(int(item["serial"]),'09d')
                item["serial"]=str(int(item["serial"])+1)
                new_log = new_log[0:m.start()-diff]+label_str+new_log[(m.end())-diff:]
                diff = diff + ((m.end())-(m.start())) - len(label_str)
                found = True
    return new_log


def compute_slcpl(thelogs, logtem):

    selected = []
    for t in logtem:

        t['template'] = recognize_known_patterns(t['template'].replace(".*","").replace("\S*","").replace("\\",""))
        t['template'] = " ".join(t['template'].split())
        t['template'] = do_tokenization([t['template']])[0] 
        u = []
        for s in t['template']:
            if s==" ":
                continue
            if s.isdigit():
                continue
            u.append(s)
        t['template'] = u
        selected.append(t)

    SL_sum = 0
    total_log_count = 0 
    for t in selected:
        static_length = len(t['template'])
        #print "   \033[1;94m", static_length, "\033[0m\033[33;36m","".join(t['template']), "\033[0m"

        #SL_sum += (static_length * t['count'])
        #SL_sum += math.sqrt((static_length * t['count'] *t['count']))

        #if t['count']>1:
        #    SL_sum += (static_length * t['count'])
        #else:
        #    print "NOTICE: Match count is only", t['count'], "for", t['template']

        #SL_sum += float(static_length * t['count']) * 5.0 * float(1.0-1.0/float(t['count']))
        SL_sum += float(static_length * t['count']) * float(1.0-1.0/float(t['count']))

        if t['count']>1:
            #print (static_length * t['count']), static_length, t['count'], total_log_count, t['template']
            print "\t",(static_length * t['count']), "\t\t",static_length, "\t\t",t['count'], "\t\t",total_log_count
        total_log_count += t['count']

    if total_log_count==0:
        print "    \033[0;91mError: total_log_count 0", "\033[0m"
        return -1.0,-1.0
    SL = float(SL_sum)/float(total_log_count) 

    total_cpl = 0
    weighted_cpl_sum = 0.0
    for i in range(0,len(selected)):
        set_len_sum = 0
        for j in range(0,len(selected)):
            if i==j: 
                continue

            n=0
            token_d = {}
            basechar="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!@#$%^&*()_+[]\{}|;:,./<>?`~=-ÇüéâäàåçêëèïîìÄÅÉæÆôöòûùÿÖÜ¢£¥₧ƒáíóúñÑªº¿⌐¬½¼¡«»░▒▓│┤╡╢╖╕╣║╗╝╜╛┐└┴┬├─┼╞╟╚╔╩╦╠═╬╧╨╤╥╙╘╒╓╫╪┘┌█▄▌▐▀αßΓπΣσµτΦΘΩδ∞φε∩≡±≥≤⌠⌡÷≈°∙·√ⁿ²■"

            for tok in selected[i]['template']:
                if n<len(basechar):
                    if tok not in token_d:
                        token_d[tok]=basechar[n]
                        n+=1
                else:
                    print "WARNING: basechar ran out!"
            for tok in selected[j]['template']:
                if n<len(basechar):
                    if tok not in token_d:
                        token_d[tok]=basechar[n]
                        n+=1
                else:
                    print "WARNING: basechar ran out!"

            str_i=""
            for tok in selected[i]['template']:
                if tok in token_d.keys():
                    str_i = str_i + token_d[tok]
            str_j=""
            for tok in selected[j]['template']:
                if tok in token_d.keys():
                    str_j = str_j + token_d[tok]

            lcs_set = lcs(str_i,str_j)

            # sum length of all the set elements
            # multiply by the 'count'
            # after the loop, divide by the total log count
            #set_len_sum += sum(len(x) for x in lcs_set)
            for k in lcs_set:
                set_len_sum += len(k)
        set_len_sum = float(set_len_sum)/float(len(selected))

#        set_len_sum = float(set_len_sum)/float(len(logtem)) # divide by the log template count because template i is compared with all the rest

        #print "\n",i,"set_len_sum=",set_len_sum
        #print "    weighted:", float(set_len_sum*selected[i]['count'])/float(total_log_count)
        weighted_cpl_sum += float(set_len_sum*selected[i]['count'])/float(total_log_count) # weighted sum of cpl

        sys.stdout.write('\r'+"    Processed "+"{0:.1f}".format(100.0*float(i)/float(len(selected)))+"%")
        sys.stdout.flush()

    print " "
    if len(selected)>1:
        CPL = float(total_cpl)/float(len(selected))/float(len(selected)-1)
    else:
        CPL=0.0
    #print "\033[1;95mAverage CPL:", CPL, "\033[0m"
    #print "\033[0;32mInverse Average max CPL:", 1.0/(1.0+CPL), "\033[0m"
    #print "\033[0;103mScore:", SL/(1.0+CPL), "\033[0m"

    return SL,weighted_cpl_sum

