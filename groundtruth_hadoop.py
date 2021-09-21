#!/usr/bin/python
#-*- coding: utf-8 -*-

import os
import re
import ast
import sys
import copy
import time
import uuid
import numpy
import argparse
import pickle
from scipy import stats
from random import randint
from collections import defaultdict

global debug_mode
global smask

tm001=0.0
tm002=0.0
tm003=0.0
tm004=0.0
tm005=0.0
tm006=0.0
tm007=0.0
tm008=0.0
tm009=0.0
tm010=0.0
tm011=0.0

UNIFORM_THRESHOLD=0.98

STAR_THRESHOLD=25
LOGLEN_THRESHOLD=int(1024*6) 

seqnum=1

MOD_FACTOR = 32

MIN_NUM_OF_TERMS=1 # minimum required number of terms in the term band, initially I introduced this to filter out one word term band.

RANDOM_SAMPLE_SIZE=4096
# frequent words must be above this cut-off line
CUTOFF_COUNT=20
CC_THRESHOLD=0.1 # select only the word-term pairs that have correlation above this value

WILDCARD_THRESHOLD=0.9 # if almost all filter vectors are filled (above this threshold), then the rest of them will be * regardless of value distribution.
FILTER_THRESHOLD=70.0 # (percentage) above this max runlength proportion, it will be added to the filter

PATTERN_THRESHOLD=3 # there must be more than this number of values to say if there is any pattern among string values

def sanitize_id(id):
    return id.strip().replace(" ", "")

(_ADD, _DELETE, _INSERT) = range(3)
(_ROOT, _DEPTH, _WIDTH) = range(3)

class Node:

    def __init__(self, name, log_count, identifier=None, expanded=True):
        self.__identifier = (str(uuid.uuid1()) if identifier is None else
                sanitize_id(str(identifier)))
        self.name = name
        self.expanded = expanded
        self.__bpointer = None
        self.__fpointer = []
        # custom states
        self.log_templates = []
        self.rep_logs = []
        self.all_vect = [-1]*log_count

    @property
    def identifier(self):
        return self.__identifier

    @property
    def bpointer(self):
        return self.__bpointer

    @bpointer.setter
    def bpointer(self, value):
        if value is not None:
            self.__bpointer = sanitize_id(value)

    @property
    def fpointer(self):
        return self.__fpointer

    def update_fpointer(self, identifier, mode=_ADD):
        if mode is _ADD:
            self.__fpointer.append(sanitize_id(identifier))
        elif mode is _DELETE:
            self.__fpointer.remove(sanitize_id(identifier))
        elif mode is _INSERT:
            self.__fpointer = [sanitize_id(identifier)]

    def is_leaf_node(self):
        if len(self.fpointer)==0:
            return True
        return False

    def print_node(self):
        print "* Node name:", self.name
        print "     log template count:", len(self.log_templates)
        print "     rep_lgos count:", len(self.rep_logs)
        print "     all_vect valid count:", sum(1 for x in self.all_vect if x>0)

class Tree:

    def __init__(self):
        self.nodes = []
        self.serial = 0
        self.coterm_set = set()

    def get_index(self, position):
        for index, node in enumerate(self.nodes):
            if node.identifier == position:
                break
        return index

    def create_node(self, name, log_count, identifier=None, parent=None):

        node = Node(name, log_count, identifier)
        self.nodes.append(node)
        self.__update_fpointer(parent, node.identifier, _ADD)
        node.bpointer = parent
        return node

    def find_inprogress_node(self):
        for node in self.nodes:
            if len(node.fpointer)==0 and -1 in node.all_vect: # if there is no child and -1 is in the all_vect, it is the unfinished leaf node.
                return node
        return None

    def find_leaf_node(self):
        for node in self.nodes:
            if len(node.fpointer)==0:
                return node
        return None

    def find_node(self, identifier):
        for node in self.nodes:
            if node.identifier==identifier:
                return node
        return None


    def show(self, position, level=_ROOT):
        queue = self[position].fpointer
        if level == _ROOT:
            print("{0} [{1}] all_vect(-1)={2} len(log_templates)={3}".format(self[position].name, self[position].identifier, self[position].all_vect.count(-1), len(self[position].log_templates)))
        else:
            print(" "*level*10, "{0} [{1}] all_vect(-1)={2} len(log_templates)={3}".format(self[position].name, self[position].identifier, self[position].all_vect.count(-1), len(self[position].log_templates)))
        if self[position].expanded:
            level += 1
            for element in queue:
                self.show(element, level)  # recursive call

    # Get the list of terms up the parents
    def linage(self, position):
        term_list = []
        term_list.append(self[position].term)

        element = self[position].bpointer
        while element is not None:
            term_list.append(self[element].term)
            element = self[element].bpointer
        return list(reversed(term_list))

    def traverse_leaf(self, position, level=_ROOT):
        queue = self[position].fpointer
        if len(queue)==0:
            data_item = sorted(self.linage(position)[1:])
            if len(data_item)>=MIN_NUM_OF_TERMS:
                self.coterm_set.add(str(data_item))
            return
        level += 1
        for element in queue:
            self.traverse_leaf(element, level)  # recursive call

    def expand_tree(self, position, mode=_DEPTH):
        # Python generator. Loosly based on an algorithm from 'Essential LISP' by
        # John R. Anderson, Albert T. Corbett, and Brian J. Reiser, page 239-241
        yield position
        queue = self[position].fpointer
        while queue:
            yield queue[0]
            expansion = self[queue[0]].fpointer
            if mode is _DEPTH:
                queue = expansion + queue[1:]  # depth-first
            elif mode is _WIDTH:
                queue = queue[1:] + expansion  # width-first

    def is_branch(self, position):
        return self[position].fpointer

    def __update_fpointer(self, position, identifier, mode):
        if position is None:
            return
        else:
            self[position].update_fpointer(identifier, mode)

    def __update_bpointer(self, position, identifier):
        self[position].bpointer = identifier

    def __getitem__(self, key):
        return self.nodes[self.get_index(key)]

    def __setitem__(self, key, item):
        self.nodes[self.get_index(key)] = item

    def __len__(self):
        return len(self.nodes)

    def __contains__(self, identifier):
        return [node.identifier for node in self.nodes if node.identifier is identifier]


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


def is_number(s):
    return s.lstrip('-').replace('.','',1).isdigit()

def are_all_numbers(numlist):
    for n in numlist:
        if not(is_number(n)):
            return False
    return True

def is_hexa(s):
    try:
        int(s, 16)
        return True
    except ValueError:
        return False

def are_all_hexa(numlist):
    for n in numlist:
        if not(is_hexa(n)):
            return False
    return True


def is_all_integer(numlist):
    for n in numlist:
        # sign check of the first character
        if not(ord(n[0]) in range(ord('0'),ord('9')+1) or n[0] in ['+','-']):
            return False

        for c in n[1:]:
            if ord(c)<48: # '0' is 48
                return False
            if ord(c)>57: # '9' is 57
                return False
    return True


def is_all_floatingpoint(numlist):
    for n in numlist:
        # sign check of the first character
        if not(ord(n[0]) in range(ord('0'),ord('9')+1) or n[0] in ['+','-']):
            return False

        point_detected = False
        for c in n[1:]:
            if not point_detected and c=='.':
                point_detected = True
                continue
            if ord(c)<48: # '0' is 48
                return False
            if ord(c)>57: # '9' is 57
                return False
        if n[-1]=='.':
            return False
        if not point_detected:
            return False
    return True



def follows_format(klist):

    global smask

    #if debug_mode:
    #    print "        [follows_format()] Entering. len=",len(klist)

    # if there are less than 3 values, it is too less to say whether there is a pattern or not
    if len(klist)<PATTERN_THRESHOLD:
        return False

    # check for string length
    slen = len(klist[0])
    for w in klist[1:]:
        if slen!=len(w):

    #        if debug_mode:
    #            print "        [follows_format()] Various string length!"

            return False
    
    # check if all entries are alphabet only
    alphabet_only = True
    for w in klist:
        if not w.isalpha():
            alphabet_only = False
    if alphabet_only: # if only alphabet, do not add as new pattern
        return False

    # eliminate any word in the form of ~100~ since these are special token I inserted
    keylist = []
    for w in klist:
        if "~" not in w:
            keylist.append(w)
    if len(keylist)<2:
        return False

    # check if there is any common characters
    seedkey = keylist[0]
    smask = ""
    fixed_count = 0
    for c in range(0,slen):
        pos_fixed = True
        for w in keylist[1:]:
            if w[c]!=seedkey[c]:
                pos_fixed = False
                break
        if pos_fixed:
            fixed_count += 1
            smask = smask + seedkey[c]
        else:
            smask = smask + "."

    # change any numbers in smask to \d
    smask = re.sub("[0-9]","\d",smask)

    if fixed_count==0:
    #    if debug_mode:
    #        print "        [follows_format()] No fixed character position detected!"
        return False

    #if debug_mode:
    #    print "        [follows_format()] detected format:", smask

    pat = { "pattern": smask, "label":"~"+str(400+len(discovered_patterns))+"~"}
    #print "\033[31;46m[follows_format] New pattern detected -\033[0m", smask

    pattern_exists = False
    for p in discovered_patterns:
        if p["pattern"]==smask:
            pattern_exists = True
            pat = p
            break
    if not pattern_exists:
        discovered_patterns.append(pat)
#        print "                 \033[31;46mNew pattern added -\033[0m", pat
#    else:
#        print "                 \033[31;46mPattern already exists.\033[0m"

    return True


def all_terms_exist(s, keywords_list):
    for term in keywords_list:
        #if term not in custom_split(s):
        if term not in s:
            return False
    return True


# kws: keyword set
def multiple_term_inclusion_count(lines, kws):
    cnt = 0
    for s in lines:
        if all_terms_exist(s, kws):
            cnt += 1
    return cnt


def print_correlation(bow,v,w1,w2):
    for i in range(0,len(bow)):

        if bow[i]!=w1:
            continue

        for j in range(0,len(bow)):
            if bow[j]!=w2:
                continue

            val = numpy.corrcoef(v[i], v[j])[0][1]
            return w1+":"+w2,"{0:.3f}".format(val)
 

standalone_patterns = [
#    {   "pattern":"\-?\d+",          # integer
#        "label": "~100~" },
#    {   "pattern":"\-?\d+(ms|msec|millisec|s|sec|second|seconds|us|microsec|KiB|GiB|MB|KB|GB|%)", # millisec, seconds, microsec ... in integer value
#        "label": "~101~" },
#    {   "pattern": "\-?\d+\.\d+",     # FP num
#        "label": "~102~" },
#    {   "pattern":"\-?\d+\.\d+(ms|msec|millisec|s|sec|second|seconds|us|microsec|KiB|GiB|MB|KB|GB|%)", # millisec, seconds, microsec ... in FP value
#        "label": "~103~" },
#    {   "pattern": "\-?\d+\.\d+%",     # FP percent
#        "label": "~104~" },
#    {   "pattern":"\d+\^\d+",         # exponent
#        "label": "~105~" },
#    {   "pattern": "0x[\da-f]+",   # hexa num
#        "label": "~106~" },
#    {   "pattern": "155\.230\.91\.\\d{3}(:\\d)?",   # IP and port
#        "label": "~107~" },

    {   "pattern":"[\da-zA-Z]{8}\-[\da-zA-Z]{4}\-[\da-zA-Z]{4}\-[\da-zA-Z]{4}\-[\da-zA-Z]{12}", # UUID format
        "label": "~108~" },
#    {'pattern': 'container_\\d{13}_\\d{4}_\\d{2}_\\d{6}', 'label': '~109~'},
#    {'pattern': 'blk_\\d{10}_\\d{4}', 'label': '~110~'},
#    {'pattern': 'application_\\d{13}_\\d{4}', 'label': '~111~'},
#    {'pattern': 'DFSClient_NONMAPREDUCE_\-?\\d+_\\d', 'label': '~112~'},
#    {'pattern': 'DFSClient_attempt_\\d+_\\d{4}_._000000_0_\\d+_1', 'label': '~113~'},
#    {'pattern': 'fsimage.ckpt_\\d{19}', 'label': '~114~'},
#    {'pattern': 'BP\-\\d{9}\-\\d+.\\d+.\\d+.\\d+\-\\d{13}', 'label': '~115~'},
#    {'pattern': 'appattempt_\\d{13}_\\d{4}_\\d{6}', 'label': '~116~'},
#    {'pattern': 'job_\\d{13}_\\d{4}', 'label': '~117~'},
#
    {'pattern': 'req\-[a-z0-9]{8}\-[a-z0-9]{4}\-[a-z0-9]{4}\-[a-z0-9]{4}\-[a-z0-9]{12}', 'label': '~115~'},

#    {'pattern': 'DFSClient_NONMAPREDUCE_\-\\d{9}_\\d', 'label': '~111~'},
#    {'pattern': 'DFSClient_attempt_\\d{13}_\\d{4}_r_\\d{6}_\\d_\\d{9}_\\d', 'label': '~112~'},
#    {'pattern': 'edits_tmp_\\d{19}-\\d{19}_\\d{19}', 'label': '~113~'},
#    {'pattern': 'application_\\d{13}_\\d{4}', 'label': '~114~'},
#    {'pattern': 'appattempt_\\d\\d\\d\\d\\d........_\\d\\d\\d\\d_\\d\\d\\d\\d\\d\\d', 'label': '~116~'},
#    {'pattern': 'deimos\\d.', 'label': '~117~'},
#    {'pattern': 'job_\\d\\d\\d\\d........._\\d\\d\\d\\d', 'label': '~119~'},
#    {'pattern': '\\d\\d\\d.\\d\\d\\d.\\d\\d.\\d\\d.', 'label': '~121~'},
#    {'pattern': '#\\d\\d....', 'label': '~122~'},
#    {'pattern': '\\d.\\d.s', 'label': '~123~'},
#    {'pattern': 'DS-........-....-\\d...-....-............', 'label': '~124~'},
#    {'pattern': 'masterappattempt_\\d\\d\\d\\d........._\\d\\d\\d\\d_\\d\\d\\d\\d\\d\\d', 'label': '~125~'},
]


# This list grows as we learn more patterns.
discovered_patterns = [
#    {   "pattern":"\-?\d+ms",
#        "label": "~105~" },
#    {'pattern': 'container_\\d\\d\\d\\d........._\\d\\d\\d\\d_\\d\\d_\\d\\d\\d...', 'label': '~106~'},
#    {'pattern': 'appattempt_\\d\\d\\d\\d\\d........_\\d\\d\\d\\d_\\d\\d\\d\\d\\d\\d', 'label': '~107~'},
#    {'pattern': 'deimos\\d.', 'label': '~108~'},
#    {'pattern': 'application_\\d\\d\\d\\d........._\\d\\d\\d\\d', 'label': '~109~'},
#    {'pattern': 'job_\\d\\d\\d\\d........._\\d\\d\\d\\d', 'label': '~110~'},
#    {'pattern': 'blk_\\d\\d\\d\\d\\d\\d\\d..._....', 'label': '~111~'},
#    {'pattern': '\\d\\d\\d.\\d\\d\\d.\\d\\d.\\d\\d.', 'label': '~112~'},
#    {'pattern': '#\\d\\d....', 'label': '~113~'},
#    {'pattern': '\\d.\\d.s', 'label': '~114~'},
#    {'pattern': 'DS-........-....-\\d...-....-............', 'label': '~115~'},
#    {'pattern': 'masterappattempt_\\d\\d\\d\\d........._\\d\\d\\d\\d_\\d\\d\\d\\d\\d\\d', 'label': '~116~'},
]


common_patterns = [
#    {   "pattern":"<=",
#        "label": "~201~" },
#    {   "pattern":"=>",
#        "label": "~202~" },
#    {   "pattern":"<-",
#        "label": "~203~" },
#    {   "pattern":"->",
#        "label": "~204~" },

#    {   "pattern":"don't",
#        "label": "~205~" },
#    {   "pattern":"won't",
#        "label": "~206~" },
#    {   "pattern":"shouldn't",
#        "label": "~207~" },
#    {   "pattern":"couldn't",
#        "label": "~208~" },
#    {   "pattern":"it's",
#        "label": "~209~" },
#    {   "pattern":"It's",
#        "label": "~210~" },
#    {   "pattern":"Didn't",
#        "label": "~211~" },
#    {   "pattern":"didn't",
#        "label": "~212~" },
#    {   "pattern":"wasn't",
#        "label": "~213~" },
#    {   "pattern":"Wasn't",
#        "label": "~214~" },

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

def preprocess_known_patterns(logs):
#    global seqnum
#    processed = []
#    label_str = "~"+format(seqnum,'09d')+"~"
#    for i in range(0,len(logs)):
#        log1 = logs[i]
#        for item in common_patterns:
#            log2 = re.sub(item["pattern"], label_str, log1)
#            if log1!=log2:
#                seqnum += 1
#                seqnum = seqnum % MOD_FACTOR
#                label_str = "~"+format(seqnum,'09d')+"~"
#            log1 = log2
#        processed.append(log1)

    processed = []
    for i in range(0,len(logs)):
        log = logs[i]
        for item in common_patterns:

            found = True
            while found:
                diff = 0
                found = False
                matched = re.finditer(item["pattern"], log)
                for m in matched:
                    label_str = item["prefix"]+"_"+format(int(item["serial"]),'09d')
                    item["serial"]=str(int(item["serial"])+1)
                    log = log[0:m.start()-diff]+label_str+log[(m.end())-diff:]
                    diff = diff + ((m.end())-(m.start())) - len(label_str)
                    found = True
        processed.append(log)
    return processed


def replace_known_patterns(tlogs):
    global seqnum

    for i in range(0,len(tlogs)):
        tlog = tlogs[i]

        for j in range(0,len(tlog)):

            w = tlog[j]
            if w in [' ',':',',','=','<','>','@']: 
                continue

            # TODO: If I enable this, "assigned" and "released" are not converted into ~10?~ format.
            #if w.isalpha():
            #    continue
            #if not any(i.isdigit() for i in w):
            #    continue

            # preprocess any troublesome characters to special marker
            if '*' in w:
                tlogs[i][j] = re.sub("\*","~200~",w)
                w = tlogs[i][j]

            for pat in standalone_patterns:
                matched = re.match("^"+pat["pattern"]+"$",w)
                if matched!=None:

                    #tlogs[i][j] = pat["label"]
                    tlogs[i][j] = "~AB"+format(seqnum,'09d')+"~"
                    seqnum += 1
                    seqnum = seqnum % MOD_FACTOR

                    #print "\033[0;35m"+matched.group(0)+"\033[0m -->", tlogs[i][j]

            is_date = True
            matched = re.match("^(\\d{4})\-(\\d{2})\-(\\d{2})$",w)
            if matched!=None:
                year = matched.group(1)
                month = matched.group(2)
                day = matched.group(3)

                if int(year)>2020:
                    is_date = False
                if int(month)>12:
                    is_date = False
                if int(day)>31:
                    is_date = False

                if is_date:
                    tlogs[i][j] = "~date_"+format(seqnum,'04d')+"~"
                    seqnum += 1
                    seqnum = seqnum % 1000
                    #print "Date format reassigned as:",tlogs[i][j]


number_patterns = [
#    {   "pattern": "\d+\.\d+\.\d+\.\d+(:\d)?", # IP address and port number
#        "type":"ipaddrport",
#        "increment":"1",
#        "serial":"1"},

    {   "pattern":"\-?\d+",          # integer
        "type":"int",
        "increment":"1",
        "serial": "1" },

    {   "pattern": "\-?\d+\.\d+",     # FP num
        "type":"float",
        "increment":"0.1",
        "serial": "0.1" },

    {   "pattern":"\-?\d+(ms|msec|millisec|s|sec|second|seconds|us|microsec|KiB|GiB|MB|KB|GB|%)", # millisec, seconds, microsec ... in integer value
        "type":"int_time",
        "increment":"1",
        "serial": "1" },

    {   "pattern":"\-?\d+\.\d+(ms|msec|millisec|s|sec|second|seconds|us|microsec|KiB|GiB|MB|KB|GB|%)", # millisec, seconds, microsec ... in FP value
        "type":"float_time",
        "increment":"0.1",
        "serial": "0.1" },

    {   "pattern":"\-?\d+\^\d+", # exponent
        "type":"exponent",
        "increment":"1",
        "serial": "1" },

    {   "pattern": "0x[\da-fA-F]+", # hexa num
        "type":"hexa1", 
        "increment":"1", # not used
        "serial": "1" }, # not used

    {   "pattern": "[\da-fA-F]+",   # hexa num
        "type":"hexa2", 
        "increment":"1", # not used
        "serial": "1" }, # not used

#    {   "pattern": "155\.230\.91\.\\d{3}(:\\d)?",   # IP and port
#        "label": "~107~" },
#    {   "pattern":"[\da-zA-Z]{8}\-[\da-zA-Z]{4}\-[\da-zA-Z]{4}\-[\da-zA-Z]{4}\-[\da-zA-Z]{12}", # UUID format
#        "label": "~108~" },
]


def uniquify_numbers(tlogs):
    for i in range(0,len(tlogs)):
        tlog = tlogs[i]
        for j in range(0,len(tlog)):
            w = tlog[j]
            if w in [' ',':',',','=','<','>']:
                continue
            if '*' in w:
                tlogs[i][j] = re.sub("\*","~200~",w)
                w = tlogs[i][j]

            found = False
            for p in number_patterns: 
                matched = re.match("^"+p["pattern"]+"$", w)
                if matched!=None:
                    found = True
                    break
            if found:
                if p["type"]=="int":
                    tlogs[i][j] = p["serial"] 
                    p["serial"] = str(int(p["serial"])+int(p["increment"])) 
                elif p["type"]=="float":
                    tlogs[i][j] = p["serial"]
                    p["serial"] = str(float(p["serial"])+float(p["increment"]))
                elif p["type"]=="int_time": 
                    tlogs[i][j] = p["serial"]+matched.group(1)
                    p["serial"] = str(int(p["serial"])+int(p["increment"]))
                elif p["type"]=="float_time":
                    tlogs[i][j] = p["serial"]+matched.group(1)
                    p["serial"] = str(float(p["serial"])+float(p["increment"]))
                elif p["type"]=="exponent":
                    tlogs[i][j] = p["serial"]
                    p["serial"] = str(int(p["serial"])+int(p["increment"]))
                elif p["type"]=="hexa1":
                    val = list(matched.group(0)[2:]) 
                    for k in range(0,len(val)):
                        c = val[k]
                        if c in ['0','1','2','3','4','5','6','7','8','9']:
                            val[k]=str(randint(0,9))
                        elif c in ['a','b','c','d','e','f']:
                            val[k]=['a','b','c','d','e','f'][randint(0,5)]
                        elif c in ['A','B','C','D','E','F']:
                            val[k]=['A','B','C','D','E','F'][randint(0,5)]
                        else:
                            print "ERROR 235 c=",c
                            sys.exit(0)
                    tlogs[i][j] = "0x"+"".join(val)
                elif p["type"]=="hexa2": 
                    '''
                    if any(i.isdigit() for i in matched.group(0)):
                        val = list(matched.group(0))
                        for k in range(0,len(val)):
                            c = val[k]
                            if c in ['0','1','2','3','4','5','6','7','8','9']:
                                val[k]=str(randint(0,9))
                            elif c in ['a','b','c','d','e','f']:
                                val[k]=['a','b','c','d','e','f'][randint(0,5)]
                            elif c in ['A','B','C','D','E','F']:
                                val[k]=['A','B','C','D','E','F'][randint(0,5)]
                            else:
                                print "ERROR 835 c=",c
                                sys.exit(0)
                        tlogs[i][j] = "".join(val)
                    '''
                    pass
                else:
                    print "ERROR 283"
                    sys.exit(0)


def apply_all_patterns(tlogs):
    global tm001

    # Convert known patterns in the token to a marker
    #print "    [apply_all_patterns()] Converting tokens of known patterns to markers using newly discovered pattern ..."
    #print "    [apply_all_patterns()] number of tokenized logs:", len(tlogs)
    tm_checkpt = time.time()
    for i in range(0,len(tlogs)):
        tlog = tlogs[i]
        for j in range(0,len(tlog)):
            word = tlog[j]

            # Two performance optimization
            if word in [' ',':',',','=','<','>']:
                continue

            # TODO: If I enable this, "assigned" and "released" are not converted into ~10?~ format.
            #if word.isalpha():
            #    continue

            #if not any(i.isdigit() for i in word):
            #    continue

            # preprocess any troublesome characters to special marker
            if '*' in word:
                tlogs[i][j] = re.sub("\*","~200~",word)
                word = tlogs[i][j]

            for p in standalone_patterns:
                matched = re.match("^"+p["pattern"]+"$",word)
                if matched!=None:
                    tlogs[i][j] = p["label"]
                    #print "\033[0;35m"+matched.group(0)+"\033[0m -->", tlogs[i][j]

    elapsed = time.time() - tm_checkpt
    tm001 += elapsed
    #print "    [apply_all_patterns()] Done converting. It took", elapsed, "seconds"


# Just apply newly added pattern to the tokenized logs instead of going through
# all the patters because regex takes time.
def apply_new_patterns(tlogs):
    global tm002
    # Convert known patterns in the token to a marker
    #print "    [apply_new_patterns()] Converting tokens of the new known pattern to markers ..."
    #print "    [apply_new_patterns()] Pattern:", standalone_patterns[-1]
    tm_checkpt = time.time()
    for i in range(0,len(tlogs)):
        tlog = tlogs[i]
        for j in range(0,len(tlog)):
            word = tlog[j]

            # Two performance optimization
            if word in [' ',':',',','=','<','>']:
                continue

            # TODO: If I enable this, "assigned" and "released" are not converted into ~10?~ format.
            #if word.isalpha():
            #    continue

            #if not any(i.isdigit() for i in word):
            #    continue

            # preprocess any troublesome characters to special marker
#            if '*' in word:
#                tlogs[i][j] = re.sub("\*","~200~",word)
#                word = tlogs[i][j]

            p = standalone_patterns[-1] # just use the newly added pattern, no need to match all
            matched = re.match("^"+p["pattern"]+"$",word)
            if matched!=None:
                tlogs[i][j] = p["label"]
                #print "\033[0;35m"+matched.group(0)+"\033[0m -->", tlogs[i][j]
    elapsed = time.time() - tm_checkpt
    tm002 += elapsed
    #print "    [apply_new_patterns()] Done converting. It took", elapsed, "seconds"


def Determine_runlen_filter_word(token_d, tlogs, ftwords):

    global new_pattern_added

    # if there is only one word, just add it to the filter_words
    if len(token_d)==1:
        return token_d.keys()[0]
    else:
        if is_all_integer(token_d.keys()):
            #print "    All integer keys detected!"
            return "*"
        elif is_all_floatingpoint(token_d.keys()):
            #print "    All floating number keys detected!"
            return "*"
        elif follows_format(token_d.keys()):
#            #print "    Keys follow certain format!", token_d.keys()[0]
#            # Insert newly detected custom pattern to the standalone_patterns dictionary
#            pat = { "pattern": smask, "label":"~"+str(100+len(standalone_patterns))+"~"}
#            pattern_exists = False
#            for p in standalone_patterns:
#                if p["pattern"]==smask:
#                    pattern_exists = True
#                    pat = p
#            if not pattern_exists:
#                standalone_patterns.append(pat)
#                new_pattern_added = True
#                if debug_mode:
#                    print "    [Determine_runlen_filter_word] \033[31;46m New pattern added \033[0m - ", pat
#                    print "    [Determine_runlen_filter_word] \033[37;41m Current filter words:\033[0m - ", ftwords
#                print "\033[31;46m New pattern added \033[0m - ", pat
#
#                # Update all the tokenized logs
#                apply_new_patterns(tlogs)
#                apply_new_patterns([ftwords]) # Update all the tokens in the filter words as well. Input should be list of list.
#
#                if debug_mode:
#                    print "    [Determine_runlen_filter_word] \033[38;42m Updated filter words:\033[0m - ", ftwords
#            else:
#                if debug_mode:
#                    print "    [Determine_runlen_filter_word] Pattern already exists - ", pat, str(ftwords)
            return "*"
        return sorted(token_d, key=lambda k: token_d[k], reverse=True)[0]


def determine_filter_word(token_d, tlen, fillup_ratio):

    pv = compute_uniformity_pvalue(token_d)
    cr = 100.0*float(len(token_d))/float(tlen)
    if debug_mode:
        print "\033[34;46mpv:"+str(pv)+"\033[0m \033[34;42mcr:"+str(cr)+"\033[0m"

    # if there is only one word, just add it to the filter_words
    if len(token_d)==1:
        if debug_mode:
            print "\033[43;5m"+"STATIC STRING because there is only one value in the dictionary."+"\033[0m"
        #print "\033[43;5m"+"STATIC STRING because there is only one value in the dictionary."+"\033[0m"
        return token_d.keys()[0], pv, cr

    #if debug_mode:
    #    print "    Tokens in the dictionary:",token_d.keys()

    if any(w in token_d for w in [" ", "@","<",">","=","(",")"]):
        if debug_mode:
            print "\033[43;5m"+"STATIC STRING because special char (including space) is in the token keys."+"\033[0m"
        #print "\033[43;5m"+"STATIC STRING because special char (including space) is in the token keys."+"\033[0m"
        return sorted(token_d, key=lambda k: token_d[k], reverse=True)[0], pv, cr

    if are_all_numbers(token_d.keys()):
        if debug_mode:
            print "\033[43;5m"+"WILDCARD because they are all numbers."+"\033[0m"
        #print "\033[43;5m"+"WILDCARD because they are all numbers."+"\033[0m"
        return '*', pv, cr
    if follows_format(token_d.keys()):
        if debug_mode:
            print "\033[43;5m"+"WILDCARD because new pattern is detected."+"\033[0m"
        #print "\033[43;5m"+"WILDCARD because new pattern is detected."+"\033[0m"
        return '*', pv, cr
    if pv>UNIFORM_THRESHOLD:
        if debug_mode:
            print "\033[43;5m"+"WILDCARD because it is a uniform distribution."+"\033[0m"
        #print "\033[43;5m"+"WILDCARD because it is a uniform distribution."+"\033[0m"
        return '*', pv, cr
    if are_all_hexa(token_d.keys()):
        if debug_mode:
            print "\033[43;5m"+"WILDCARD because they are all hexadecimal numbers."+"\033[0m"
        #print "\033[43;5m"+"WILDCARD because they are all hexadecimal numbers."+"\033[0m"
        return '*', pv, cr

    token =  sorted(token_d, key=lambda k: token_d[k], reverse=True)[0]
    if '~' in token:
        if debug_mode:
            print "\033[43;5m"+"WILDCARD because it is a known pattern."+"\033[0m"
        #print "\033[43;5m"+"WILDCARD because it is a known pattern."+"\033[0m"
        return '*', pv, cr

    #if fillup_ratio!=None and fillup_ratio>WILDCARD_THRESHOLD:
    #    print "\033[43;5m"+"WILDCARD because fill-up ratio is reached."+"\033[0m", fillup_ratio
    #    return "*", pv, cr

    if debug_mode:
        print "\033[43;5m"+"STATIC STRING because it did not meet any condition for the wildcard."+"\033[0m"
    #print "\033[43;5m"+"STATIC STRING because it did not meet any condition for the wildcard."+"\033[0m"

    return token, pv, cr


def read_log_files(flist, filter_str):
    logs = []
    for i in range(0,len(flist)):
        for log in flist[i]:
            if len(log.strip())==0:
                continue

            if filter_str != None:
                if filter_str not in log:
                    continue
                print log.strip()

            log = " ".join(log.split())
            logs.append(log)


    print "Total number of logs loaded:",len(logs)
    return logs


def match_and_remove(tmpl,logs):

    global tm011
    tm_checkpt = time.time()

    # first, build a list of index to delete
    to_delete = []
    match_count = 0
    for i in range(0,len(logs)):
        log = logs[i]
        matched = re.match("^"+tmpl+"$",log)

        if matched!=None:
            match_count = match_count + 1
            to_delete.append(i)
            #if debug_mode:
            #    print "DEL:",log

    # delete matched logs
    before_removal = len(logs)
    to_delete = sorted(to_delete)
    for i in reversed(sorted(to_delete)):
        del logs[i]
    del_count = before_removal - len(logs)
    if match_count!=del_count:
        print "Error 263: match count and delete count mismatch!"
        sys.exit(0)

    elapsed = time.time() - tm_checkpt
    tm011+=elapsed
    return del_count


def exist_match(log_template, logs):
    for i in range(0,len(logs)):
        matched = re.match("^"+log_template+"$", logs[i])
        if matched!=None:
            if debug_mode:
                print "\033[0;32mMatch found at "+str(i)+":", logs[i], "\033[0m "
            return i
    return -1


def test_multiple_match(rlogs, vect, log_template):

    cnt = 0
    to_delete = []
    for j in range(0,len(rlogs)):
        if vect[j]==-1:
            continue
        log = rlogs[j]
        matched = re.match("^"+log_template+"$", "".join(log))
        if matched!=None:
            #print "      ", "".join(log)
            to_delete.append(j)
            cnt += 1
    return cnt, to_delete


def mark_matched_logs(logs, vect, rlogs, log_template, i):
    #print "Entering mark_matched_logs() Star_count:",log_template.count(".*")
    marked = 0 # how many logs match to the log template?
    replog_selected = False
    for j in range(0,len(logs)):

        if vect[j]>-1: # skip logs already matched by previous templates
            continue

        log = logs[j]

        # TODO if log is too long, it takes too long to match the regular expression even though the number of wildcard is OK.
        # I am shortening the log
        if len(log)>LOGLEN_THRESHOLD:
            log = log[:LOGLEN_THRESHOLD]

        matched = re.match("^"+log_template+"$",log)
        if matched!=None:
            vect[j] = i
            marked += 1
            
            if not replog_selected:
                replog_selected = True
                rlogs.append(log)

#    if sum(1 for x in vect if x>0) != len(rlogs):
#        print "ERROR: vect("+str(sum(1 for x in vect if x>0))+") and rlogs length("+str(len(rlogs))+") does not match!!"
#        sys.exit(0)

    #print "Leaving mark_matched_logs()"
    return marked


def remove_log_template_matches(logs, logtem):
    sample_logs = []
    # Remove matching logs using pre-filled log templates
    for i in range(0,len(logtem)):
        log_template = logtem[i]
        # Remove any matched logs from the logs.
        before_removal = len(logs)
        alog = None
        for j in reversed(range(0,len(logs))):
            matched = re.match("^"+log_template+"$",logs[j])
            if matched!=None:
                alog = logs[j]
                del logs[j]
        removed_logs = before_removal - len(logs)

        if removed_logs==0:
            print "ERROR: no matching logs found from the given template."
            print "template:", log_template
            sys.exit(0)

        print "\033[0;31m"+"["+format(i,'3d')+"]",format(len(logs),'5d'),format(removed_logs,'4d'),"\033[0m","\033[0;32m\""+logtem[i]+"\",\033[0m"
        #print "["+format(i,'3d')+"]",format(len(logs),'5d'),format(removed_logs,'4d'), logtem[i]

        sample_logs.append(alog)
    print "Done removing logs using pre-populated log templates. Remaining logs:", len(logs)
    return sample_logs


def build_random_index(data_len, sample_len):
    # create random index list
    if data_len<=sample_len:
        numbers = range(0,data_len)
    else:
        numbers = set()
        #while len(numbers)<sample_len:
        #    numbers.add(randint(0,data_len-1))

        # TODO deterministically generate the number list
        num = 0
        while len(numbers)<sample_len:
            numbers.add(num)
            num += 1

    return numbers


def sample_by_length(logs, vect, ssize):
    global tm006
    tm_checkpt = time.time()

    tlen_d = defaultdict()
    for log in logs:
        toklen = len(log.split())
        if toklen not in tlen_d:
            tlen_d[toklen]= []
        tlen_d[toklen].append(log)
    most_popular = sorted(tlen_d, key=lambda k: len(tlen_d[k]), reverse=True)[0]


    for x in tlen_d:
        print x, len(tlen_d[x])

    sys.exit(0)

    selected=[]
    for log in logs:
        if len(log.split())==most_popular:
            selected.append(log)
        if len(selected)>=ssize:
            elapsed = time.time() - tm_checkpt
            #print "{0:.3f}".format(elapsed), "Random log selection"
            tm006+=elapsed
            return selected
    sys.exit(0)


def sample_by_token_length_and_space_count(logs, tlogs, vect):
    global tm006
    tm_checkpt = time.time()

    logscore_d = defaultdict()
    tlog_d = defaultdict()
    for i in range(0,len(logs)):
        if vect[i]>-1: 
            continue
        log = logs[i]
        tlog = tlogs[i]

        toklen = len(log.split())
        space_cnt = log.count(' ') + log.count('\t') + log.count(',') + log.count(':') + log.count(';') + log.count(',')
        score = toklen*1000+space_cnt

        if score not in logscore_d:
            logscore_d[score]= []
        logscore_d[score].append(log)
        if score not in tlog_d:
            tlog_d[score]= []
        tlog_d[score].append(tlog)
    most_popular = sorted(logscore_d, key=lambda k: len(logscore_d[k]), reverse=True)[0]

    if debug_mode:
        print "** Summary of log groups using characters **"
        for x in sorted(logscore_d, key=lambda k: len(logscore_d[k]), reverse=True)[:20]:
            print "  For the key of",format(x,'5d')+",", format(len(logscore_d[x]), '5d'),"logs are grouped."
        print "    ..."

    selected = []
    for log in logscore_d[most_popular]:
        selected.append(log)
        if len(selected)>=1000:
            break
    tselected = []
    for tlog in tlog_d[most_popular]:
        tselected.append(tlog)
        if len(tselected)>=1000:
            break

    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed), "Random log selection"
    tm006+=elapsed

    return selected, tselected


def sample_by_term_correlation(logs, ssize):
    global tm008
    tm_checkpt = time.time()

    numset = build_random_index(len(logs), ssize)

    rand_logs = random_sample_logs(logs, ssize)
    toke_logs = do_tokenization(rand_logs)
    bow,bow_list = select_significant_terms(toke_logs)
    all_vect = build_term_vectors(bow_list,toke_logs)
    corr_d = compute_term_correlation(all_vect, bow, bow_list, max(len(x) for x in bow.keys()))
    while corr_d==None:
        rand_logs = random_sample_logs(logs, ssize)
        toke_logs = do_tokenization(rand_logs)
        bow,bow_list = select_significant_terms(toke_logs)
        all_vect = build_term_vectors(bow_list,toke_logs)
        corr_d = compute_term_correlation(all_vect, bow, bow_list, max(len(x) for x in bow.keys()))
        print "corr_d None. Looping one more time."

    term_groups = determine_term_groups(corr_d,bow,rand_logs)
    tb = sorted(term_groups, key=len, reverse=True)[0]

    selected=[]
    tselected=[]
    for log in logs:
        if all_terms_exist(log, tb):
            selected.append(log)
            tselected.append(custom_split(log))
        if len(selected)>=ssize:
            elapsed = time.time() - tm_checkpt
            #print "{0:.3f}".format(elapsed), "Term correlation-based selection"
            tm008+=elapsed
            return selected,tselected

    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed), "Term correlation-based selection"
    tm008+=elapsed

    return selected,tselected


def sample_by_signature(logs, ssize):
    global tm007
    tm_checkpt = time.time()

    numset = build_random_index(len(logs), ssize)

    # dictionary of signature
    tlen_d = defaultdict()
    for n in numset:
        signature = re.sub("[\d\s\w]","",logs[n])
        if signature not in tlen_d:
            tlen_d[signature]=0
        tlen_d[signature] += 1
    most_popular = sorted(tlen_d, key=lambda k: tlen_d[k], reverse=True)[0]

    #for x in sorted(tlen_d, key=lambda k: tlen_d[k], reverse=True):
    #    print x, tlen_d[x]
    #sys.exit(0)

    selected=[]
    for log in logs:
        if re.sub("[\d\s\w]","",log)==most_popular:
            selected.append(log)
        if len(selected)>=ssize:

            elapsed = time.time() - tm_checkpt
            #print "{0:.3f}".format(elapsed), "Random log selection"
            tm007+=elapsed
            return selected

    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed), "Random log selection"
    tm007+=elapsed

    return selected

def random_sample_logs(logs, n):
    global tm003
    # Randomly select RANDOM_SAMPLE_SIZE logs from the original set
    # To avoid selecting the same log, I first collect the set of RANDOM_SAMPLE_SIZE index and then create a log list.
    if len(logs)<=n:
        return logs
    '''
    tm_checkpt = time.time()
    numset = build_random_index(len(logs), ssize)
    selected = []
    for pos in numset:
        selected.append(logs[pos])
    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed), "Random log selection"
    '''

    tm_checkpt = time.time()
    selected = []
    step_size = len(logs) / RANDOM_SAMPLE_SIZE
    for log in logs[::step_size]:
        if len(selected)>=RANDOM_SAMPLE_SIZE:
            break
        selected.append(log)
    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed), "Random log selection"

    tm003+=elapsed
    return selected


def do_tokenization(logs):
    global tm004
    tm_checkpt = time.time()
    tlogs = []
    for i in range(0,len(logs)):
        log = logs[i]
        tlogs.append(custom_split(log))
    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed), "Tokenizing"
    tm004+=elapsed
    return tlogs


def select_significant_terms(tlog_data):
    # Build bag-of-words with their frequencies
    tm_checkpt = time.time()
    wd = defaultdict()
    for tlog in tlog_data:
        for w in tlog:
            if not w.isalpha():
                continue
            if len(w)<=2:
                continue
            if w not in wd:
                wd[w] = 0
            wd[w] += 1
    # Delete infrequent words
    for w in wd:
        if wd[w] < CUTOFF_COUNT:
            del wd[w]
    wl = sorted(wd, key=lambda k: wd[k], reverse=True)
    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed),"Creating bag-of-words took"
    if debug_mode:
        print "The length of bag-of-words list:",len(wl)
    return wd, wl


def build_term_vectors(wl,tlog_data):
    # Build word vectors for all the words
    # One vector indicates whether the word exists in that log line
    # Thus, the vector length is equal to the length of tlog_data
    tm_checkpt = time.time()
    vectors = []
    for i in range(0,len(wl)):
        w = wl[i]
        ivect = [0]*len(tlog_data)
        for j in range(0,len(tlog_data)):
            tlog = " "+" ".join(tlog_data[j])+" "
            if " "+w+" " in tlog:
                ivect[j]=1
        #print str(bow[w]),w+" "*(40-len(w)),"".join(str(x) for x in ivect)
        vectors.append(ivect)
    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed),"Creating indicator vectors for all words"
    if debug_mode:
        print "Length of vectors:",len(vectors)
    return vectors


def compute_term_correlation(vect, bow, bow_list, max_word_length):
    # Compute the correlations of all word pairs into correlation_dict dictionary
    tm_checkpt = time.time()
    correlation_dict = defaultdict()
    for i in range(0,len(vect)):
        word = bow_list[i]
        if word not in correlation_dict:
            correlation_dict[word] = defaultdict()
        for j in range(0,len(vect)):
            term = bow_list[j]
            if term not in correlation_dict:
                correlation_dict[term] = defaultdict()
            if word==term:
                correlation_dict[word][term]=1.0
                correlation_dict[term][word]=1.0
                continue

            if term not in correlation_dict[word] and word not in correlation_dict[term]:
                val = numpy.corrcoef(vect[i], vect[j])[0][1]
            elif term not in correlation_dict[word] and word in correlation_dict[term]:
                val = correlation_dict[term][word]
            elif term in correlation_dict[word] and word not in correlation_dict[term]:
                val = correlation_dict[word][term]
            else:
                val = correlation_dict[term][word]

            if val>CC_THRESHOLD: # or val<-1.0*CC_THRESHOLD:
                #print bow_list[j]+" "*(50-len(bow_list[j])),"\t","{0:.3f}".format(val)
                #corr_list.append("{0:.3f}".format(val)+":"+bow_list[j])
                correlation_dict[word][term] = val
                correlation_dict[term][word] = val
    # Print word correlations
#    for word in sorted(correlation_dict, key=lambda k: bow[k], reverse=True):
#        cc_list="["
#        for term in sorted(correlation_dict[word], key=lambda k: correlation_dict[word][k], reverse=True):
#            if word==term:
#                continue
#            cc_list = cc_list+term+":"+"{0:.3f}".format(correlation_dict[word][term])+", "
#        cc_list += "]" 
#        print "*",format(bow[word],'4d'),word+" "*(max_word_length-len(word)), cc_list
    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed), "Correlation computation"

    # check if correlation_dict is empty or not
    clen = 0
    for word in sorted(correlation_dict, key=lambda k: bow[k], reverse=True):
        clen += len(correlation_dict[word])
        #print "        [compute_term_correlation] length of "+word, len(correlation_dict[word]), correlation_dict[word]
    if clen==len(correlation_dict):
        print "WARNING 219: correlation_dict is empty"
        #for t in log_templates:
        #    print "\033[0;43m"+str(t["count"])+"\033[0m",t["template"]
        print "Total number of templates:",len(log_templates)
        return None

    #print print_correlation(bow_list,vect,"Stopping","Exit")
    return correlation_dict


def display_term_groups(correlation_dict, word_dict):
    # Determine the correlated term group
    tm_checkpt = time.time()
    tgrp = [] # term group
    circ = []
    max_word_length = max(len(x) for x in word_dict.keys())
    for word in sorted(correlation_dict, key=lambda k: word_dict[k], reverse=True):
        # if word is already part of any of the previous groups, skip it
        is_member = False
        for g in tgrp:
            if word in g:
                is_member = True
        if is_member:
            continue
        for term in sorted(correlation_dict[word], key=lambda k: correlation_dict[word][k], reverse=True):
            # append new term only if its correlation coefficient to all the members are above threshold
            eligible = True
            for member in circ:
                #print "member:",member," term:",term, " member keys:", correlation_dict[member].keys()
                if term not in correlation_dict[member]:
                    eligible = False
                    break
                val = correlation_dict[member][term]
                if val<CC_THRESHOLD:
                    eligible = False
                    break
            if eligible: # the word itself is added automatically as the first member here because correlation_dict has itself
                circ.append(term)
        if len(circ)>=2:
            print word, " "*(max_word_length-len(word)), format(len(correlation_dict[word]),'3d'), format(len(circ),'3d'), circ
            tgrp.append(sorted(circ))
        # circle members determined at this point
        circ = []
    elapsed = time.time() - tm_checkpt
    if debug_mode:
        print "Displaying term groups took", elapsed, "seconds"
    return


def determine_term_groups(correlation_dict, word_dict, logs):
    # Determine the correlated term group
    tm_checkpt = time.time()
    tgrp = [] # term group
    circ = []
    max_word_length = max(len(x) for x in word_dict.keys())
    for word in sorted(correlation_dict, key=lambda k: word_dict[k], reverse=True):
        # if word is already part of any of the previous groups, skip it
        is_member = False
        for g in tgrp:
            if word in g:
                is_member = True
        if is_member:
            continue
        for term in sorted(correlation_dict[word], key=lambda k: correlation_dict[word][k], reverse=True):
            # append new term only if its correlation coefficient to all the members are above threshold
            eligible = True
            for member in circ:
                #print "member:",member," term:",term, " member keys:", correlation_dict[member].keys()
                if term not in correlation_dict[member]:
                    eligible = False
                    break
                val = correlation_dict[member][term]
                if val<CC_THRESHOLD:
                    eligible = False
                    break
            if eligible: # the word itself is added automatically as the first member here because correlation_dict has itself
                circ.append(term)
        if len(circ)>=2:
            # counting inclusion of all terms is very costly
            #print word, " "*(max_word_length-len(word)), format(len(correlation_dict[word]),'3d'), format(len(circ),'3d'), format(multiple_term_inclusion_count(logs,circ),'5d'), circ
            if debug_mode:
                print word, " "*(max_word_length-len(word)), format(len(correlation_dict[word]),'3d'), format(len(circ),'3d'), circ
            tgrp.append(sorted(circ))
        # circle members determined at this point
        circ = []
    elapsed = time.time() - tm_checkpt
    #print "{0:.3f}".format(elapsed), "Determining term groups"
    return tgrp 


def do_filtering(tlogs, valid_vect, flt_words, flt_valid_vect):

    # Filter logs based on the flt_words
    filtered = []
    for i in range(0,len(tlogs)):

        if valid_vect[i]==0:
            continue

        tlog = tlogs[i]

        # go through the tokens for each log and select only the ones that match the filter words in active positions
        add_ok = True
        for j in range(0,len(flt_valid_vect)):

            if flt_valid_vect[j]==0:
                continue
            if flt_words[j]=='*':
                continue

            if j>len(tlog)-1: 
                add_ok = False
                break
            if tlog[j] != flt_words[j]: 
                add_ok = False
                break

        if add_ok:
            filtered.append(tlog)
        else:
            valid_vect[i] = 0

    # TODO: if there is 0 filtered logs, we need to drop some filter words.

    return filtered
 

def generate_log_template(fwords):

    # Transform the discovered log template into a python-ready form
    log_template = "".join(fwords).strip()

    log_template = re.sub("\\\\","\\\\\\\\",log_template )
    log_template = re.sub("\-","\-",log_template)
    log_template = re.sub("\[","\[",log_template)
    log_template = re.sub("\]","\]",log_template)
    log_template = re.sub("\(","\(",log_template)
    log_template = re.sub("\)","\)",log_template)
    log_template = re.sub("\$","\$",log_template)
    log_template = re.sub("\?","\?",log_template)
    log_template = re.sub("\+","\+",log_template)
    log_template = re.sub("\|","\|",log_template)
    #log_template = re.sub("\\\\","~201~",log_template)

    log_template = re.sub("\*","\S+",log_template)

    for it in standalone_patterns:
        log_template = re.sub(it["label"],"\S+",log_template)
    
    log_template = re.sub("~200~","\*",log_template)

    # recover special common word patterns
    #for n in range(0,len(common_patterns)):
    #    if common_patterns[n]["label"]=="~324~" or common_patterns[n]["label"]=="~325~":
    #        log_template = re.sub(common_patterns[n]["label"], "\S+ \S+", log_template)
    #    else:
    #        log_template = re.sub(common_patterns[n]["label"], "\S+", log_template)
    #log_template = re.sub("~300~","\S+",log_template)


    log_template = re.sub("\\S+~","\S+",log_template)

    log_template = re.sub("\\\\\[\\\\\]","\\[\S*\\]",log_template)
    log_template = re.sub("{}","{\S*}",log_template)

    final_template = []
    for t in log_template.split():
        if "\S+" in t and "=\S+" not in t and ":\S+" not in t: 
            final_template.append("\S+")
        else:
            final_template.append(t)
    log_template = " ".join(final_template)

    found = True
    while found:
        diff = 0
        found = False
        matched = re.finditer("(\\\\S\+)[:]\\\\S\+", log_template)
        for m in matched:
            log_template = log_template[0:m.start()-diff]+"\S+"+log_template[m.end()-diff:]
            diff = diff + 4
            found = True

    while "\S+\S+" in log_template:
        log_template = re.sub("\\\\S\+\\\\S\+", "\S+", log_template)


    #log_template = re.sub("\\\\S\+ \\\\S\+",".*",log_template)
    #log_template = re.sub("\\\\S\+,\\\\S\+",".*",log_template)
    #log_template = re.sub("\\\\S\+:\\\\S\+",".*",log_template)

    #print "Log Sample:  ", "\033[0;31m"+log_sample+"\033[0m"
    #print "Log template:", "\033[0;47m"+log_template+"\033[0m"
    return log_template


def generate_log_template_star(fwords):

    # Transform the discovered log template into a python-ready form
    log_template = "".join(fwords).strip()

    for item in number_patterns:

        if item["type"]=="hexa1" or item["type"]=="hexa2": 
            continue

        found = True
        while found:
            diff = 0
            found = False
            matched = re.finditer(item["pattern"], log_template)
            for m in matched:


                log_template = log_template[0:m.start()-diff]+".*"+log_template[(m.end())-diff:]
                diff = diff + ((m.end())-(m.start())) - len(".*")
                found = True


    log_template = re.sub("\\\\","\\\\\\\\",log_template )
    log_template = re.sub("\-","\-",log_template)
    log_template = re.sub("\[","\[",log_template)
    log_template = re.sub("\]","\]",log_template)
    log_template = re.sub("\(","\(",log_template)
    log_template = re.sub("\)","\)",log_template)
    log_template = re.sub("\$","\$",log_template)
    log_template = re.sub("\?","\?",log_template)
    log_template = re.sub("\+","\+",log_template)
    log_template = re.sub("\|","\|",log_template)
    #log_template = re.sub("\\\\","~201~",log_template)

    #log_template = re.sub("\*","\S+",log_template)
    log_template = re.sub("\*",".*",log_template)

    for it in standalone_patterns:
        #log_template = re.sub(it["label"],"\S+",log_template)
        log_template = re.sub(it["label"],".*",log_template)
    
    log_template = re.sub("~200~","\*",log_template)

    # recover special common word patterns
    #for n in range(0,len(common_patterns)):
    #    log_template = re.sub(common_patterns[n]["label"], ".*", log_template)
    #log_template = re.sub("~300~",".*",log_template)

    log_template = re.sub("~\\S+~",".*",log_template)

    for p in common_patterns:
        marker = p["prefix"]+"_"+"\\d{9}"
        log_template = re.sub(marker,".*",log_template)

    log_template = re.sub("\\\\\[\\\\\]","\\[.*\\]",log_template)
    log_template = re.sub("{}","{.*}",log_template)

    final_template = []
    for t in log_template.split():
        if ".*" in t and "=.*" not in t and ":.*" not in t: 
            final_template.append(".*")
        else:
            final_template.append(t)
    log_template = " ".join(final_template)

    found = True
    while found:
        diff = 0
        found = False
        matched = re.finditer("(\.\*)[:]\.\*", log_template)
        for m in matched:
            log_template = log_template[0:m.start()+1-diff]+log_template[(m.end()-1)-diff:]
            diff += 3
            found = True

    while ".* .*" in log_template:
        log_template = re.sub("\.\* \.\*", ".*", log_template)

    while "..*" in log_template:
        log_template = re.sub("\.\.\*", ".*", log_template)

    while ".*.*.*.*" in log_template:
        log_template = re.sub("\.\*\.\*\.\*\.\*", ".*", log_template)
    while ".*.*.*" in log_template:
        log_template = re.sub("\.\*\.\*\.\*", ".*", log_template)
    while ".*.*" in log_template:
        log_template = re.sub("\.\*\.\*", ".*", log_template)

    if log_template.count(".*") > STAR_THRESHOLD:
        parts=log_template.split(".*")
        log_template=""
        for i in range(STAR_THRESHOLD):
            log_template+=parts[i]
            log_template+=".*"

    return log_template


def postprocess_raw_template(fwords):
    logtem= "".join(fwords).strip()
    for it in standalone_patterns:
        logtem= re.sub(it["label"],"*",logtem)
    logtem = re.sub("~200~","*",logtem)
    logtem = re.sub("~300~","*",logtem)
    return logtem 


def compute_uniformity_pvalue(val_dict):
    tensor=[] 
    for n in sorted(val_dict, key=lambda k: val_dict[k], reverse=False):
        if len(tensor)==0:
            cumul = 0
        else:
            cumul = tensor[-1]
        tensor.append(cumul+val_dict[n]) 
    pval = stats.kstest(tensor,'uniform',args=(tensor[0],tensor[-1])).pvalue
    return pval


def exist_partial_match2(rlogs, fwords, fmask):
    survived = []
    for i in range(0,len(rlogs)):
        tlog = rlogs[i]
        # go through the tokens for each log and select only the ones that match the filter words in active positions
        add_ok = True
        for j in range(0,len(fmask)):
            if fmask[j]==0:
                continue
            if j>len(tlog)-1: 
                add_ok = False
                break
            if fwords[j]=='*':
                continue
            if '~' in fwords[j]:
                continue
            if tlog[j] != fwords[j]:
                add_ok = False
                break
        if add_ok:
            survived.append(tlog)
    print "    ",len(survived)
    return survived
 
def exist_partial_match(rlogs, rvect, fwords, fmask):
    print "Entering exist_partial_match."
    print "Filter words:",fwords
    print "Filter masks:",fmask

    for k in range(0,len(rlogs)):
        if rvect[k]==-1: 
            continue
        #print k, rvect[k], rlogs[k]

        #print "+",tokenized
        #print "+", "".join(str(x) for x in fmask)
        tokenized = rlogs[k] 

        #print "Repre log:",tokenized

        match_found = True
        for m in range(0,len(fmask)): 
            if fmask[m]==0:
                continue
            if m >= len(tokenized):
                #print "Difference found!! Length too short."
                match_found = False
                break
            #print "fwords[m]="+fwords[m], "m="+str(m), "len(tokenized)="+str(len(tokenized)), "match_found="+str(match_found)

            if '~' in fwords[m]:
                continue

            if fwords[m]!=tokenized[m]:
                #print "Difference found!!", fwords[m], tokenized[m]
                match_found = False
                break

        if match_found:
            #print "Leaving exist_partial_match() with True! Checked", k, "logs."
            return True
        #    print "===> Match found"
        #    print rlogs[k]
        #    print "".join(fwords)
        #    print "".join(str(x) for x in fmask)
    #print "Leaving exist_partial_match() with False! Checked", k, "logs."
    return False


def finalize_filter_with_star(fword,fmask):
    tfword = copy.deepcopy(fword)
    tfmask = copy.deepcopy(fmask)
    for i in range(0,len(tfmask)):
        if tfmask[i]==0:
            tfword[i]="*"
            tfmask[i]=1
    return tfword, tfmask


def forge_candidate_log_templates(input_logs, rep_logs):

    global tm005
    global debug_mode

    valid_mask = [1]*len(input_logs)

    column_cnt = max(len(x) for x in input_logs)
    filter_words = [""]*column_cnt
    filter_mask = [0]*column_cnt

    tm_checkpt = time.time()
    token_added_order = []
    count_added_order = [] 
    candidate_set = [] 
    add_candidate = False
    while 0 in filter_mask: # only when filter selection vector is not full

        filtered_logs = do_filtering(input_logs, valid_mask, filter_words, filter_mask)

        if debug_mode:
            print "Updating column_cnt from", column_cnt,"to",max(len(x) for x in filtered_logs)
            for x in filtered_logs:
                print "    ",x
        column_cnt = max(len(x) for x in filtered_logs)
        filter_words = filter_words[:column_cnt]
        filter_mask = filter_mask[:column_cnt]

        if 0 not in filter_mask: 
            break

        max_runlen = 0
        max_runlen_pct = 0.0
        max_runlen_pos = -1

        all_column_dict = defaultdict() # dict of dict, column_dict per tpos is saved here
        for tpos in range(0,column_cnt):

            # Skip the token position if it has been added to the filter words already
            if filter_mask[tpos]==1:
                continue

            column_dict = defaultdict()
            for tlog in filtered_logs: # tlog is the filtered and tokenized log lines
                if len(tlog) > tpos: 
                    tok = tlog[tpos]
                    if tok not in column_dict: # if key is not yet created, make one
                        column_dict[tok] = 0
                    column_dict[tok] += 1 # increment count
            all_column_dict[tpos] = column_dict

            if len(column_dict)==0: 
                print "ERROR 832: No token values collected into the dictionary. Perhaps log's length ran out."
                sys.exit(0)

            runlength_token = sorted(column_dict, key=lambda k: column_dict[k], reverse=True)[0]
            runlength = column_dict[runlength_token] 

            runlen_percent = float(runlength*100.0)/float(len(filtered_logs))
            if runlen_percent<100.0 and runlen_percent>max_runlen_pct: 
                max_runlen = runlength
                max_runlen_pct = runlen_percent
                max_runlen_pos = tpos
            #print tpos, max_runlen_pos, "->"+runlength_token+"<-", runlen_percent

            if runlength==len(filtered_logs):
                filter_mask[tpos] = 1
                filter_words[tpos] = runlength_token
                if debug_mode:
                    print "\033[0;36mAdding single-valued column to the filter\033[0m [tpos:"+str(tpos)+"]", "->"+runlength_token+"<-", runlength
                #print "\033[0;36mAdding single-valued column to the filter\033[0m [tpos:"+str(tpos)+"]", "->"+runlength_token+"<-", runlength
                token_added_order.append(runlength_token)
                count_added_order.append(1) 

        if 0 not in filter_mask: 
            if debug_mode:
                print "Exiting loop since all filters are determined.", filter_mask
                print "    filter_words:", filter_words
            break

        if max_runlen_pos==-1: 
            print "ERROR: max column not selected!!!"
            sys.exit(0)

        if max_runlen_pos>=0: 

            target_dict = all_column_dict[max_runlen_pos] 
            filled = float(sum(filter_mask))/float(len(filter_mask)) 
            new_fword, pv, cr = determine_filter_word(target_dict, len(filtered_logs), filled)
            #print "{0:.5f}".format(pv), "{0:.5f}".format(cr), "---->"+new_fword

            if debug_mode:
                print "max_column:"+str(max_runlen_pos)+",\033[35;47mdetermine_filter_word returned:\033[0m", "->"+new_fword+"<-"
                
                print "*** Max token from each column ***"
                for h in all_column_dict: # h is a column position
                    d = all_column_dict[h]
                    for n in sorted(d, key=lambda k: d[k], reverse=True):
                        print "["+str(h)+"]",d[n], "\t","->"+n+"<-"
                        break
            if debug_mode and ' ' not in target_dict and '=' not in target_dict:
                print "------------------------------------------------------------"
                print "[Column:"+str(max_runlen_pos)+"]", "\033[1;91mpval:",pv,"\033[0m", "Cardinality:", "{0:.2f}".format(cr),"%"
                print "------------------------------------------------------------"
                print "num  |   count   |       token      "
                print "------------------------------------------------------------"
                # print each line
                cnt = 1
                for n in sorted(target_dict, key=lambda k: target_dict[k], reverse=True):
                    print "["+str(cnt)+"]\t",target_dict[n], "\t\t","->"+n+"<-"
                    cnt += 1
                    if cnt>40:
                        print "..."
                        break

            if new_fword=="*":
                if pv<UNIFORM_THRESHOLD: 
                    filter_words,filter_mask = finalize_filter_with_star(filter_words,filter_mask)
                else: 
                    
                    tw,tm = finalize_filter_with_star(filter_words,filter_mask)

                    log_template = generate_log_template_star(tw)
                    #print "log template:", log_template
                    if exist_match(log_template, rep_logs)>=0:
                        new_fword = sorted(target_dict, key=lambda k: target_dict[k], reverse=True)[0]
                    else:
                        filter_words,filter_mask = finalize_filter_with_star(filter_words,filter_mask)

            else: 
                if pv>0.9 and pv<UNIFORM_THRESHOLD: 
                    
                    if add_candidate==False: 
                        add_candidate = True
                else: 
                    pass

                if add_candidate:
                    tw,tm = finalize_filter_with_star(filter_words,filter_mask)
                    log_template = generate_log_template_star(tw)
                    #print "log template:", log_template
                    #print "\033[1;95mCandidate:", "\033[0m \033[90;102m", "".join(tw), "\033[0m "
                    if exist_match(log_template, rep_logs)<0:
                        candidate_set.append(log_template)
                        #print "\033[1;94mCandidate ACCEPTED", "\033[0m ", "\033[1;95mCandidate:", "\033[0m \033[90;102m", "".join(tw), "\033[0m "
                    else:
                        #print "\033[1;91mCandidate REJECTED", "\033[0m ","\033[1;95mCandidate:", "\033[0m \033[90;102m", "".join(tw), "\033[0m "
                        pass

            if new_fword!="*": 
                
                filter_mask[max_runlen_pos] = 1
                filter_words[max_runlen_pos] = new_fword
                token_added_order.append(new_fword)
                count_added_order.append(len(target_dict))
            #print "\033[0;35mNew filter word\033[0m [tpos:"+str(max_runlen_pos)+"]", "=>"+new_fword+"<="

        if debug_mode:
            print "Current column_cnt:",column_cnt
            print "Max runlength percent:", "{0:.2f}".format(max_runlen_pct),"%"
            print "Max runlength percent position:", max_runlen_pos
            print "\033[1;94mMax runlength percent word :", "->"+filter_words[max_runlen_pos]+"<-\033[0m"
            print "sum of valid_mask:", sum(valid_mask)
            print "filter_mask(sum:"+str(sum(filter_mask))+"/"+str(len(filter_mask))+"):","".join(str(x) for x in filter_mask)
            print "Filtering logs using filter_words:",filter_words, len(filter_words)

            #print "Added order:", "\033[1;95m|\033[0m".join(token_added_order)
            #print "Cardinality order:", "\033[1;95m|\033[0m".join(count_added_order)
            print "\033[1;95mToken list in the added order:(The number is the count of unique tokens.)\033[0m"
            for w in range(0,len(token_added_order)):
                print "        ", format(count_added_order[w],'3d'),token_added_order[w]

            raw_input("\033[0;35m->Press ENTER to continue filtering ...\033[0m")
            print " "

        #print "filter_vect(sum:"+str(sum(filter_mask))+"/"+str(len(filter_mask))+"):","".join(str(x) for x in filter_mask)
    # END of inner while loop

    elapsed = time.time() - tm_checkpt
    tm005 += elapsed
    #print "{0:.3f}".format(elapsed),"Filter construction going through all columns"

    # All filter_mask is filled now.
    # Generate a log template
    log_template = generate_log_template_star(filter_words)
    #log_template = generate_log_template(filter_words)

    candidate_set.append(log_template)

    return candidate_set
    #return log_template


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


# thelogs will be reduced by matching log templates
# logtem is a list of 'count','template' dict
def compute_slcpl(thelogs, logtem):


    selected = []
    for t in sorted(logtem, key=lambda k: k['count'], reverse=True):
        # first, build a list of index to delete
        to_delete = []
        for i in range(0,len(thelogs)):
            log = thelogs[i]
            matched = re.match("^"+t["template"]+"$",log)

            if matched!=None:
                to_delete.append(i)

        # delete matched logs
        before_removal = len(thelogs)
        to_delete = sorted(to_delete)
        for i in reversed(sorted(to_delete)):
            del thelogs[i]
        del_count = before_removal - len(thelogs)

        #print "Removed",del_count,"logs."
        #if del_count==0:
        #    print "\033[33;31m",format(del_count,'5d'), format(int(t['count']),'5d'), t['template'], "\033[0m"

        if del_count>0:
            #t['template'] = t['template'].replace(".*","")
            t['template'] = do_tokenization([t['template'].replace(".*","")])[0] 
            selected.append(t)
            #print t['template']

    #print len(thelogs),"logs remaining."
    #print "Initial log template count:", len(logtem)
    #print "Selected log template count:", len(selected)

#    selected =  [
#        { 
#            'template': "abcd abcd ",
#            'count': 10
#        },
#        { 
#            'template': "abcd ccc ",
#            'count': 10
#        },
#        { 
#            'template': "abcd ddd",
#            'count': 10
#        },
#    ]

    SL_sum = 0
    total_log_count = 0 
    for t in selected:

        u = []
        for s in t['template']:
            if s==" ":
                continue
            u.append(s)
        t['template'] = u

        static_length = len(t['template'])
        #print "  \033[1;94m", static_length, "\033[0m", "\033[33;36m","".join(t['template']), "\033[0m"
        SL_sum += (static_length * t['count'])
        total_log_count += t['count']
    if total_log_count==0:
        print "Error: total_log_count 0"
        sys.exit(0)
    SL = float(SL_sum)/float(total_log_count) 
    print "\033[0;32mAverage weighted SL:", SL, "\033[0m"

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

            if len(selected[i]['template'])>len(basechar):
                print "\n\n    Number of tokens in the template:",len(selected[i]['template'])
                print "    Length of basechar:",len(basechar)
                print "    "+str(selected[i]['template'])
                print "    \033[1;31mTemplate too long. Truncating to match the basechar length.\033[0m"
                selected[i]['template'] = selected[i]['template'][0:len(basechar)-1]
                #sys.exit(0)

            if len(selected[j]['template'])>len(basechar):
                print "\n\n    Number of tokens in the template:",len(selected[j]['template'])
                print "    Length of basechar:",len(basechar)
                print "    "+str(selected[j]['template'])
                print "    \033[1;31mTemplate too long. Truncating to match the basechar length.\033[0m"
                selected[j]['template'] = selected[j]['template'][0:len(basechar)-1]
                #sys.exit(0)

            for tok in selected[i]['template']:
                if tok not in token_d:
                    token_d[tok]=basechar[n]
                    n+=1
            for tok in selected[j]['template']:
                if tok not in token_d:
                    token_d[tok]=basechar[n]
                    n+=1

            if n>len(basechar):
                print "\033[0;32mERROR: Not enough char set!!", n, len(basechar), "\033[0m"
                sys.exit(0)

            str_i=""
            for tok in selected[i]['template']:
                str_i = str_i + token_d[tok]
            str_j=""
            for tok in selected[j]['template']:
                str_j = str_j + token_d[tok]

            lcs_set = lcs(str_i,str_j) 

            #total_cpl +=  sum(len(x) for x in lcs_set)/len(lcs_set) # not used
            # sum length of all the set elements
            # multiply by the 'count'
            # after the loop, divide by the total log count
            for k in lcs_set:
                set_len_sum += len(k)
        set_len_sum = float(set_len_sum)/float(len(selected))

#        set_len_sum = float(set_len_sum)/float(len(logtem)) # divide by the log template count because template i is compared with all the rest

        #print "\n",i,"set_len_sum=",set_len_sum
        #print "    weighted:", float(set_len_sum*selected[i]['count'])/float(total_log_count)
        weighted_cpl_sum += float(set_len_sum*selected[i]['count'])/float(total_log_count) # weighted sum of cpl

        sys.stdout.write('\r'+"Processed "+"{0:.1f}".format(100.0*float(i)/float(len(selected)))+"%")
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



def mark_matched_logs2(logs, mask, template, verbose=False):
    marked=0
    for i in range(0, len(logs)):
        if mask[i]>-1:
            continue

        if verbose:
            if i%(len(logs)/30)==0:
                sys.stdout.write('\r'+"    \033[1;91m "+"{0:.1f}".format(float(i*100)/float(len(logs)))+"% \033[0m"+template)
                sys.stdout.flush()

        log = logs[i]
        first_len = len(log)
        if len(log)>LOGLEN_THRESHOLD: 
            log = log[:LOGLEN_THRESHOLD-3072] 
            if verbose:
                print "\n\033[0;32m"+template+"\033[0m"
                print "\033[1;91m<"+str(i)+"/"+str(len(logs))+":"+"{0:.1f}".format(float(i*100)/float(len(logs)))+"%>\033[0m",log
                print first_len,"->",len(log)
        else:
            if verbose:
                print first_len

        #if len(log)==LOGLEN_THRESHOLD:
        #    new_tem=template.replace(".*","\S*")
        #    template=new_tem

        matched = re.match('^'+template+'$', log)
        if matched!=None:
            mask[i]=99
            marked+=1
    return marked


########################## ########################## ########################## ##########################
########################## ########################## ########################## ##########################

log_templates = []

prepopulated_log_templates = [
"INFO org.apache.hadoop.yarn.server.resourcemanager.rmcontainer.RMContainerImpl: .* Container Transitioned from .* to .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Memory usage of ProcessTree .* for container\-id .*: .* of .* physical memory used; .* of .* virtual memory used",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.application.ApplicationImpl: Application .* transitioned from .* to .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.container.ContainerImpl: Container .* transitioned from .* to .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.AuxServices: Got event .* for appId .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.ContainerManagerImpl: Start request for .* by user .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.ContainerManagerImpl: Stopping container with container Id: .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.launcher.ContainerLaunch: Cleaning up container .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.localizer.ResourceLocalizationService: Writing credentials to the nmPrivate file .*\. Credentials list:",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Stopping resource-monitoring for .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.DefaultContainerExecutor: Deleting absolute path : .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.NMAuditLogger: USER=.* IP=.* OPERATION=.* TARGET=.* RESULT=.* APPID=.* CONTAINERID=.*",
"INFO org.apache.hadoop.yarn.server.nodemanager.NMAuditLogger: USER=.* OPERATION=.* Finished \- Succeeded TARGET=.* RESULT=.* APPID=.* CONTAINERID=.*",
"INFO org.apache.hadoop.yarn.server.nodemanager.NMAuditLogger: USER=.* OPERATION=.* Finished \- Killed TARGET=.* RESULT=.* APPID=.* CONTAINERID=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.RMAuditLogger: USER=.* OPERATION=.* TARGET=.* RESULT=.* APPID=.* CONTAINERID=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.RMAuditLogger: USER=.* IP=.* OPERATION=.* TARGET=.* RESULT=.* APPID=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.RMAuditLogger: USER=.* OPERATION=.* Finished \- Succeeded TARGET=.* RESULT=.* APPID=.*",
"INFO org.apache.hadoop.yarn.server.nodemanager.NodeStatusUpdaterImpl: Removed completed containers from NM context: \[.*\]",
"INFO org.apache.hadoop.yarn.server.nodemanager.NodeStatusUpdaterImpl: Registering with RM using containers :\[.*\]",
"INFO org.apache.hadoop.yarn.server.nodemanager.NodeStatusUpdaterImpl: Sending out 0 NM container statuses: \[.*\]",
"INFO org.apache.hadoop.yarn.server.nodemanager.security.NMContainerTokenSecretManager: Rolling master-key for container-tokens, got key with id .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.NodeStatusUpdaterImpl: Registered with ResourceManager as .*:.* with total resource of <memory:.*, vCores:.*>",
"INFO org.apache.hadoop.yarn.server.resourcemanager.rmapp.attempt.RMAppAttemptImpl: .* State change from .* to .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.application.ApplicationImpl: Adding .* to application .*",
"INFO org.apache.hadoop.mapred.ShuffleHandler: Added token for .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.AuxServices: Adding auxiliary service httpshuffle, \".*\"",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.AuxServices: Got .* for service mapreduce_shuffle",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.application.ApplicationImpl: Removing .* from application .*",
"INFO SecurityLogger.org.apache.hadoop.ipc.Server: Auth successful for .* \(auth:.*\)",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.ParentQueue: Re\-sorting assigned queue: .* stats: default: capacity=.*, absoluteCapacity=.*, usedResources=<memory:.*, vCores:.*>, usedCapacity=.*, absoluteUsedCapacity=.*, numApps=.*, numContainers=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.LeafQueue: completedContainer container=Container: \[ContainerId: .*, NodeId: .*:.*, NodeHttpAddress: .*:.*, Resource: <memory:.*, vCores:.*>, Priority: .*, Token: Token { kind: .*, service: .*:.*}, \] queue=default: capacity=.*, absoluteCapacity=.*, usedResources=<memory:.*, vCores:.*>, usedCapacity=.*, absoluteUsedCapacity=.*, numApps=.*, numContainers=.* cluster=<memory:.*, vCores:.*>",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.ParentQueue: assignedContainer queue=.* usedCapacity=.* absoluteUsedCapacity=.* used=<memory:.*, vCores:.*> cluster=<memory:.*, vCores:.*>",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Starting resource\-monitoring for .*",
"INFO org.apache.hadoop.mapred.ShuffleHandler: Setting connection close header\.\.\.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.common.fica.FiCaSchedulerApp: Completed container: .* in state: .* event:.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Application attempt .* released container .* on node: host: .* #containers=.* available=<memory:.*, vCores:.*> used=<memory:.*, vCores:.*> with event: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.SchedulerNode: Released container .* of capacity <memory:.*, vCores:.*> on host .*, which currently has .* containers, <memory:.*, vCores:.*> used and <memory:.*, vCores:.*> available, release resources=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.ParentQueue: Re\-sorting completed queue: .* stats: default: capacity=.*, absoluteCapacity=.*, usedResources=<memory:.*, vCores:.*>, usedCapacity=.*, absoluteUsedCapacity=.*, numApps=.*, numContainers=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.ParentQueue: completedContainer queue=.* usedCapacity=.* absoluteUsedCapacity=.* used=<memory:.*, vCores:.*> cluster=<memory:.*, vCores:.*>",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.LeafQueue: assignedContainer application attempt=.* container=.*: \[ContainerId: .* NodeId: .*, NodeHttpAddress: .*, Resource: <memory:.*, vCores:.*>, Priority: .* Token: null, \] queue=.*: capacity=.*, absoluteCapacity=.*, usedResources=<memory:.*, vCores:.*>, usedCapacity=.*, absoluteUsedCapacity=.*, numApps=.*, numContainers=.* clusterResource=<memory:.*, vCores:.*> type=.* requestedPartition=",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.LeafQueue: default used=<memory:.*, vCores:.*> numContainers=.* user=.* user\-resources=<memory:.*, vCores:.*>",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.SchedulerNode: Assigned container .* of capacity <memory:.*, vCores:.*> on host .*, which has .* containers, <memory:.*, vCores:.*> used and <memory:.*, vCores:.*> available after allocation",
"INFO org.apache.hadoop.yarn.server.nodemanager.DefaultContainerExecutor: launchContainer: \[.*\]",
"WARN org.apache.hadoop.yarn.server.nodemanager.DefaultContainerExecutor: Exit code from container .* is : .*",
"INFO org.apache.hadoop.yarn.event.AsyncDispatcher: Registering class .* for class .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetAsyncDiskService: Deleted .* file .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetAsyncDiskService: Scheduling .* file .* for deletion",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Neither virutal\-memory nor physical\-memory monitoring is needed. Not running the monitor\-thread",
"INFO org.apache.hadoop.http.HttpServer2: Jetty bound to port .*",
"INFO org.apache.hadoop.ipc.CallQueueManager: Using callQueue: class .* queueCapacity: .*",
"INFO org.apache.hadoop.ipc.Server: Starting Socket Reader .* for port .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.DefaultContainerExecutor: Copying from .* to .*",
"INFO org.apache.hadoop.metrics2.impl.MetricsSystemImpl: Scheduled snapshot period at .* second\(s\)\.",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.localizer.ResourceLocalizationService: Created localizer for .*",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.BlockManager: .* = .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Opened streaming server at .*",
"INFO org.apache.hadoop.yarn.factories.impl.pb.RpcServerFactoryPBImpl: Adding protocol .* to the server",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: .* = .*",
"INFO org.apache.hadoop.util.GSet: Computing capacity for map .*",
"INFO org.apache.hadoop.hdfs.server.namenode.TransferFsImage: Transfer took .* at .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Balancing bandwith is .* bytes/s",
"INFO org.apache.hadoop.util.GSet: VM type = .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Starting DataNode with maxLockedMemory = .*",
"INFO org.apache.hadoop.hdfs.server.common.Storage: Locking is disabled for .*",
"INFO org.apache.hadoop.metrics2.impl.MetricsConfig: loaded properties from .*",
"INFO org.apache.hadoop.ipc.Server: IPC Server Responder: starting",
"INFO org.apache.hadoop.hdfs.server.datanode.BlockScanner: Initialized block scanner with targetBytesPerSec .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Got finalize command for block pool .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Receiving .* src: .* dest: .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.localizer.ResourceLocalizationService: Localizer started on port .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSEditLog: Starting log segment at .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Scanning block pool .* on volume .*",
"INFO org.apache.hadoop.hdfs.StateChange: STATE\* Safe mode is OFF",
"INFO org.apache.hadoop.hdfs.StateChange: STATE\* Safe mode extension entered\.",
"INFO org.apache.hadoop.hdfs.StateChange: STATE\* Safe mode ON\.",
"INFO org.apache.hadoop.hdfs.StateChange: STATE\* UnderReplicatedBlocks has .* blocks",
"INFO org.apache.hadoop.hdfs.StateChange: STATE\* Network topology has .* racks and .* datanodes",
"INFO org.apache.hadoop.hdfs.StateChange: STATE\* Leaving safe mode after .* secs",
"INFO org.apache.hadoop.http.HttpRequestLog: Http request log for .* is not defined",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImage: Start loading edits file .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.ContainerManagerImpl: ContainerManager started at .*",
"INFO org.apache.hadoop.yarn.webapp.WebApps: Web app node started at .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Opened IPC server at .*",
"INFO org.apache.hadoop.ipc.Server: IPC Server listener on .* starting",
"INFO org.apache.hadoop.hdfs.server.datanode.web.DatanodeHttpServer: Listening HTTP traffic on .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Refresh request received for nameservices: .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImage: Loaded image for txid .* from .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.ContainerManagerImpl: ContainerManager bound to .*",
"INFO org.apache.hadoop.hdfs.server.common.Storage: Lock on .* acquired by nodename .*",
"INFO org.apache.hadoop.hdfs.server.namenode.TransferFsImage: Downloaded file .* size .* bytes\.",
"INFO org.apache.hadoop.hdfs.server.namenode.TransferFsImage: Opening connection to .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Number threads for balancing is .*",
"INFO org.apache.hadoop.net.NetworkTopology: Adding a new node: .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: PacketResponder: .*, type=.* terminating",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Configured hostname is .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.recovery.RMStateStore: Storing info for app: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.rmapp.RMAppImpl: .* State change from .* to .* on event=.*",
"INFO org.apache.hadoop.http.HttpServer2: adding path spec: .*",
"INFO org.mortbay.log: Logging to .* via .*",
"INFO org.mortbay.log: Extract jar:file:.* to .*",
"INFO org.mortbay.log: jetty\-6\.1\.26",
"INFO org.mortbay.log: Started .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImage: Reading .* expecting start txid .*",
"INFO org.apache.hadoop.http.HttpServer2: Added global filter '.*' \(class=.*\)",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.localizer.LocalizedResource: Resource .* transitioned from .* to .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.recovery.RMStateStore: Updating info for app: .*",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.DatanodeDescriptor: Number of failed storage changes from .* to .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.ContainerManagerImpl: Creating a new application reference for app .*",
"INFO org.apache.hadoop.conf.Configuration: found resource .* at file:.*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSDirectory: Maximum size of an xattr: .*",
"INFO org.apache.hadoop.metrics2.impl.MetricsSystemImpl: SecondaryNameNode metrics system started",
"INFO org.apache.hadoop.hdfs.server.namenode.NameNode: Caching file names occuring more than .* times",
"INFO org.apache.hadoop.metrics2.impl.MetricsSystemImpl: ResourceManager metrics system started",
"INFO org.apache.hadoop.mapreduce.v2.hs.JobHistory: Starting scan to move intermediate done files",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Using ResourceCalculatorPlugin : .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Physical memory check enabled: .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Retry cache on namenode is enabled",
"INFO org.apache.hadoop.hdfs.StateChange: DIR\* completeFile: .* is closed by .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode.clienttrace: src: .* dest: .* bytes: .* op: .*, cliID: .* offset: .* srvID: .* blockid: .*, duration: .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.DefaultContainerExecutor: Initializing user .*",
"INFO org.apache.hadoop.yarn.server.webproxy.WebAppProxyServlet: dr.who is accessing unchecked .* which is the app master GUI of .* owned by .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.launcher.ContainerLaunch: Container .* succeeded",
"INFO org.apache.hadoop.http.HttpServer2: Added filter .* \(class=.*\) to context .*",
"INFO org.apache.hadoop.util.GSet: capacity = .* = .* entries",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: dnUserName = .*",
"INFO org.apache.hadoop.util.GSet: .* max memory .* = .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.webapp.WebServer: Instantiating NMWebApp at .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.localizer.ResourceLocalizationService: per directory file limit = .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Starting BPOfferServices for nameservices: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ApplicationMasterService: AM registration .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Adding block pool .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.DefaultContainerExecutor: Localizer CWD set to .* = file:.*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Registered .* MBean",
"INFO org.apache.hadoop.hdfs.server.common.Storage: Analyzing storage directories for bpid .*",
"INFO org.apache.hadoop.hdfs.server.namenode.top.window.RollingWindowManager: topN size for command .* is: .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Cached dfsUsed found for .*",
"INFO org.apache.hadoop.yarn.client.RMProxy: Connecting to ResourceManager at .*",
"INFO org.apache.hadoop.hdfs.server.namenode.NameCache: initialized with .* entries .* lookups",
"INFO org.apache.hadoop.security.token.delegation.AbstractDelegationTokenSecretManager: Updating the current master key for generating delegation tokens",
"INFO org.apache.hadoop.hdfs.server.namenode.FSDirectory: XAttrs enabled\? .*",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.DatanodeDescriptor: Adding new storage ID .* for DN .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ResourceManager: Transitioning to active state",
"INFO org.apache.hadoop.hdfs.server.datanode.DirectoryScanner: Periodic Directory Tree Verification scan starting at .* with interval of .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: No KeyProvider found\.",
"INFO org.apache.hadoop.yarn.webapp.WebApps: Registered webapp guice modules",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: supergroup = supergroup",
"INFO org.apache.hadoop.hdfs.server.namenode.FSDirectory: ACLs enabled\? .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.NMTokenSecretManagerInRM: Rolling master\-key for nm\-tokens",
"INFO org.apache.hadoop.security.authentication.server.AuthenticationFilter: Unable to initialize .*, falling back to use random secrets\.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.AppSchedulingInfo: checking for deactivate of application :.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.AMRMTokenSecretManager: Creating password for .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Added new volume: .*",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.BlockManager: Increasing replication from .* to .* for .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.security.NMContainerTokenSecretManager: Updating node address : .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSEditLog: Ending log segment .*",
"INFO org.apache.hadoop.hdfs.server.namenode.NameNode: NameNode RPC up at: .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImageFormatProtobuf: Loaded FSImage in .* seconds.",
"INFO org.apache.hadoop.mapred.IndexCache: IndexCache created with max memory = .*",
"INFO .*: STARTUP_MSG:",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Application Attempt .* is done. finalState=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.AppSchedulingInfo: Application .* requests cleared",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: fsLock is fair: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.rmnode.RMNodeImpl: .* Node Transitioned from .* to .*",
"INFO org.apache.hadoop.mapred.ShuffleHandler: httpshuffle listening on port .*",
"WARN org.apache.hadoop.hdfs.server.common.Util: Path .* should be specified as a URI in configuration files. Please update hdfs configuration.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacitySchedulerConfiguration: max alloc vcore per queue for root.default is undefined",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ApplicationMasterService: .* unregistered successfully.",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Finished loading FSImage in .*",
"INFO org.apache.hadoop.hdfs.server.namenode.LeaseManager: Number of blocks under construction: .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Adding replicas to map for block pool .* on volume .*",
"INFO .* NameNode metrics system started",
"INFO org.apache.hadoop.yarn.server.nodemanager.security.NMTokenSecretManagerInNM: Rolling master\-key for container\-tokens, got key with id .*",
"INFO org.apache.hadoop.hdfs.server.namenode.top.metrics.TopMetrics: NNTop conf: dfs.namenode.top.num.users = .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacitySchedulerConfiguration: max alloc mb per queue for root.default is undefined",
"INFO org.apache.hadoop.hdfs.server.namenode.Checkpointer: Checkpointer about to load edits from .* stream\(s\).",
"INFO org.apache.hadoop.yarn.server.resourcemanager.recovery.RMStateStore: Storing RMDTMasterKey.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.amlauncher.AMLauncher: Cleaning master .*",
"INFO org.apache.hadoop.metrics2.impl.MetricsSystemImpl: NodeManager metrics system started",
"INFO org.apache.hadoop.metrics2.impl.MetricsSystemImpl: DataNode metrics system started",
"INFO BlockStateChange: BLOCK\* addToInvalidates: .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Virtual memory check enabled: false",
"INFO org.apache.hadoop.yarn.server.nodemanager.NodeStatusUpdaterImpl: Initialized nodemanager for null: physical\-memory=.* virtual\-memory=.* virtual\-cores=.*",
"INFO org.apache.hadoop.hdfs.server.namenode.NameNode: fs.defaultFS is .*",
"INFO org.apache.hadoop.hdfs.server.namenode.top.metrics.TopMetrics: NNTop conf: dfs.namenode.top.window.num.buckets = .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSEditLog: Rolling edit logs",
"INFO org.apache.hadoop.hdfs.server.namenode.SecondaryNameNode: Web server init done",
"INFO .*: registered UNIX signal handlers for \[.*\]",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Block pool <registering> \(Datanode Uuid unassigned\) service to .* starting to offer service",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: HA Enabled: .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Detailed lock hold time metrics enabled: .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Roll Edit Log from .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ApplicationMasterService: Unregistering app attempt : .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Added volume \- .* StorageType: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.LeafQueue: Initializing default",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: initializing replication queues",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.RMDelegationTokenSecretManager: storing master key with keyID .*",
"INFO org.apache.hadoop.yarn.util.RackResolver: Resolved .* to /default\-rack",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImageFormatPBINode: Loading .* INodes.",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImage: Edits file .* of size .* edits # .* loaded in .* seconds",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.CacheReplicationMonitor: Starting CacheReplicationMonitor with interval .* milliseconds",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.DatanodeManager: dfs.block.invalidate.limit=.*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Virtual memory check enabled: .*",
"INFO org.apache.hadoop.hdfs.server.namenode.NNStorageRetentionManager: Purging old image FSImageFile\(file=.*, cpktTxId=.*\)",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.localizer.ResourceLocalizationService: usercache path : file:.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.RMNMInfo: Registered .* MBean",
"ERROR org.apache.hadoop.hdfs.server.datanode.DirectoryScanner: dfs.datanode.directoryscan.throttle.limit.ms.per.sec set to value below .* ms/sec. Assuming default value of .*",
"INFO org.apache.hadoop.yarn.webapp.WebApps: Web app cluster started at .*",
"INFO org.apache.hadoop.hdfs.DFSUtil: Starting Web\-server for hdfs at: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ResourceManager: Transitioned to active state",
"INFO org.apache.hadoop.hdfs.server.namenode.FSEditLog: Number of transactions: .* Total time for transactions\(ms\): .* Number of transactions batched in Syncs: .* Number of syncs: .* SyncTimes\(ms\): .*",
"INFO org.apache.hadoop.hdfs.server.namenode.NNStorageRetentionManager: Going to retain .* images with txid >= .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Registered .* MBean",
"INFO org.apache.hadoop.hdfs.server.namenode.SecondaryNameNode: Image has changed. Downloading updated image from NN.",
"INFO org.apache.hadoop.hdfs.server.namenode.FileJournalManager: Finalizing edits file .* \-> .*",
"WARN org.apache.hadoop.hdfs.server.namenode.SecondaryNameNode: Checkpoint done. New Image Size: .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.ContainerManagerImpl: Blocking new container\-requests as container manager rpc server is still starting.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.RMContainerTokenSecretManager: Rolling master\-key for container\-tokens",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: Using ResourceCalculatorProcessTree : null",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Total time to add all replicas to map: .*",
"INFO org.apache.hadoop.hdfs.server.common.Storage: Using .* threads to upgrade data directories \(dfs.datanode.parallel.volumes.load.threads.num=.*, dataDirs=.*\)",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.BlockManager: dfs.block.access.token.enable=.*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Namenode Block pool .* \(Datanode Uuid .* service to .* trying to claim ACTIVE state with txid=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ApplicationMasterService: Registering app attempt : .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.NMTokenSecretManagerInRM: Sending NMToken for nodeId : .* for container : .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.rmapp.RMAppImpl: Storing application with id .*",
"INFO org.apache.hadoop.hdfs.StateChange: .* allocate .*{UCState=.*, truncateBlock=.*, primaryNodeIndex=.*, replicas=\[ReplicaUC\[\[DISK\].*:NORMAL:.*\], ReplicaUC\[\[DISK\].*:NORMAL:.*\], ReplicaUC\[\[DISK\].*:NORMAL:.*\]\]} for .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Acknowledging ACTIVE Namenode Block pool .* \(Datanode Uuid .*\) service to .*/.*:.*",
"INFO org.apache.hadoop.yarn.server.nodemanager.NodeStatusUpdaterImpl: Notifying ContainerManager to unblock new container\-requests",
"INFO BlockStateChange: .* processReport .* from storage .* node .* datanodeUuid=.*, infoPort=.*, infoSecurePort=.*, ipcPort=.*, storageInfo=lv=.*;cid=.*;nsid=.*;c=.*\), blocks: .* hasStaleStorage: .*, processing time: .* msecs",
"INFO org.apache.hadoop.yarn.server.resourcemanager.amlauncher.AMLauncher: Launching .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Successfully sent block report .* containing .* storage report\(s\), of which we sent .*\. The reports had .* total blocks and used .* RPC\(s\)\. This took .* msec to generate and .* msecs for RPC and NN processing\. Got back one command: .*\.",
"INFO org.apache.hadoop.security.Groups: clearing .* cache",
"INFO org.apache.hadoop.yarn.server.resourcemanager.metrics.SystemMetricsPublisher: YARN system metrics publishing service is not enabled",
"INFO org.apache.hadoop.security.token.delegation.AbstractDelegationTokenSecretManager: Starting expired delegation token remover thread, tokenRemoverScanInterval=.* min\(s\)",
"INFO org.apache.hadoop.hdfs.server.datanode.VolumeScanner: .* no suitable block pools found to scan. Waiting .* ms.",
"INFO org.apache.hadoop.hdfs.server.namenode.NameNode: RPC server is binding to .*:.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Added node .* clusterResource: <memory:.*, vCores:.*>",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Total time to scan all replicas for block pool .*: .*ms",
"INFO org.apache.hadoop.yarn.server.nodemanager.containermanager.loghandler.NonAggregatingLogHandler: Scheduling Log Deletion for application: .* with delay of .* seconds",
"INFO org.apache.hadoop.hdfs.server.namenode.SecondaryNameNode: Log Size Trigger :.* txns",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.NMTokenSecretManagerInRM: NMTokenKeyRollingInterval: .* and NMTokenKeyActivationDelay: .*",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.BlockManager: dfs.namenode.startup.delay.block.deletion.sec is set to .*",
"INFO org.apache.hadoop.hdfs.StateChange: BLOCK\* registerDatanode: from DatanodeRegistration\(.*:.*, datanodeUuid=.*, infoPort=.*, infoSecurePort=.*, ipcPort=.*, storageInfo=lv=.*;cid=.*;nsid=.*;c=.*\) storage .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Time taken to scan block pool .* on .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Block pool Block pool .* \(Datanode Uuid null\) service to .* successfully registered with NN",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Setting up storage: nsid=.*;bpid=.*;lv=.*;nsInfo=lv=.*;cid=.*;nsid=.*;c=.*;bpid=.*;dnuuid=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.NMTokenSecretManagerInRM: Clear node set for .*",
"INFO org.apache.hadoop.yarn.server.nodemanager.DefaultContainerExecutor: Deleting path : file:.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacitySchedulerConfiguration: max alloc vcore per queue for root is undefined",
"INFO org.apache.hadoop.yarn.security.YarnAuthorizationProvider: org.apache.hadoop.yarn.security.ConfiguredYarnAuthorizer is instiantiated.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ClientRMService: Application with id .* submitted by user .*",
"INFO org.apache.hadoop.hdfs.server.namenode.TransferFsImage: Uploaded image with txid .* to namenode at .* in .* seconds",
"INFO org.apache.hadoop.yarn.server.resourcemanager.amlauncher.AMLauncher: Setting up container Container: \[ContainerId: .* NodeId: .*, NodeHttpAddress: .*, Resource: <memory:.*, vCores:.*>, Priority: .* Token: Token { kind: .*, service: .* }, \] for AM .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Accepted application .* from user: .*, in queue: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.ParentQueue: Application removed \- appId: .* user: .* leaf\-queue of parent: root #applications: .*",
"INFO org.apache.hadoop.hdfs.server.namenode.EditLogInputStream: Fast\-forwarding stream .* to transaction ID .*",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.DatanodeManager: dfs.namenode.datanode.registration.ip\-hostname\-check=.*",
"INFO org.apache.hadoop.mapreduce.v2.jobhistory.JobHistoryUtils: Default file system \[hdfs://.*:.*\]",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Append Enabled: true",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImageFormatProtobuf: Saving image file .* using no compression",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Starting services required for active state",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.AMRMTokenSecretManager: Application finished, removing password for .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.ParentQueue: Initialized parent\-queue root name=.*, fullname=.*",
"INFO org.apache.hadoop.hdfs.server.namenode.top.metrics.TopMetrics: NNTop conf: dfs.namenode.top.windows.minutes = .*",
"WARN org.apache.hadoop.hdfs.server.datanode.DataNode: Slow BlockReceiver write data to disk cost:.* \(threshold=.*\)",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImage: Planning to load image: FSImageFile\(file=.*, cpktTxId=.*\)",
"INFO org.apache.hadoop.hdfs.server.namenode.NameNode: createNameNode .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacitySchedulerConfiguration: max alloc mb per queue for root is undefined",
"INFO org.apache.hadoop.hdfs.server.namenode.NameNode: Clients are to use .*:.* to access this namenode/service.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.LeafQueue: Application .* from user: .* activated in queue: .*",
"INFO org.apache.hadoop.mapreduce.v2.hs.JobHistory: History Cleaner complete",
"INFO org.apache.hadoop.mapreduce.v2.hs.JobHistory: History Cleaner started",
"INFO org.apache.hadoop.hdfs.server.namenode.TransferFsImage: Image Transfer timeout configured to .* milliseconds",
"INFO org.apache.hadoop.yarn.server.resourcemanager.rmapp.RMAppImpl: Updating application .* with final state: .*",
"WARN org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Only one image storage directory \(dfs.namenode.name.dir\) configured. Beware of data loss due to lack of redundant storage directories!",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Added Application Attempt .* to scheduler from user .* in queue .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.rmapp.attempt.RMAppAttemptImpl: Updating application attempt .* with final state: .*, and exit status: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.AbstractYarnScheduler: Release request cache is cleaned up",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.ParentQueue: Application added \- appId: .* user: .* leaf\-queue of parent: root #applications: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.ParentQueue: root, capacity=.*, asboluteCapacity=.*, maxCapacity=.*, asboluteMaxCapacity=.*, state=.*, acls=.*, labels=.*,",
"INFO org.apache.hadoop.hdfs.DFSUtil: Starting Web\-server for secondary at: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.AMRMTokenSecretManager: AMRMTokenKeyRollingInterval: .* and AMRMTokenKeyActivationDelay: .* ms",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Retry cache will use .* of total heap and retry cache entry expiry time is .* millis",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.LeafQueue: Application removed \- appId: .* user: .* queue: .* #user\-pending\-applications: .* #user\-active\-applications: .* #queue\-pending\-applications: .* #queue\-active\-applications: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.LeafQueue: Application added \- appId: .* user: .* leaf\-queue: default #user\-pending\-applications: .* #user\-active\-applications: .* #queue\-pending\-applications: .* #queue\-active\-applications: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.AMRMTokenSecretManager: Create .* for ApplicationAttempt: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.security.RMContainerTokenSecretManager: ContainerTokenKeyRollingInterval: .* and ContainerTokenKeyActivationDelay: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Initialized queue mappings, override: .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Initialized CapacityScheduler with calculator=.* org.apache.hadoop.yarn.util.resource.DefaultResourceCalculator, minimumAllocation=<<memory:.*, vCores:.*>>, maximumAllocation=<<memory:.*, vCores:.*>>, asynchronousScheduling=.*, asyncScheduleInterval=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Initialized queue: default: capacity=.*, absoluteCapacity=.*, usedResources=<memory:.*, vCores:.*>, usedCapacity=.*, absoluteUsedCapacity=.*, numApps=.*, numContainers=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.rmapp.attempt.RMAppAttemptImpl: Storing attempt: AppId: .* AttemptId: .* MasterContainer: Container: \[ContainerId: .* NodeId: .*, NodeHttpAddress: .*, Resource: <memory:.*, vCores:.*>, Priority: .* Token: Token { kind: ContainerToken, service: .* }, \]",
"INFO org.apache.hadoop.hdfs.StateChange: .* Replication Queue initialization scan for invalid, over\- and under\-replicated blocks completed in .* msec",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Null container completed\.\.\.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ResourceTrackerService: NodeManager from node .* httpPort: .* registered with capability: <memory:.*, vCores:.*>, assigned nodeId .*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: For namenode .* using .* of .*msec .* of .*msec Initial delay: .* heartBeatInterval=.*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Block pool .* \(Datanode Uuid null\) service to .* beginning handshake with NN",
"WARN org.apache.hadoop.yarn.server.nodemanager.containermanager.AuxServices: The Auxilurary Service named .* in the configuration is for class .* which has a name of .*\. Because these are not the same tools trying to send ServiceData and read Service Meta Data may have issues unless the refer to the name in the config\.",
"INFO BlockStateChange: BLOCK\* BlockManager: ask .* to delete \[.*\]",
"INFO BlockStateChange: .* addStoredBlock: blockMap updated: .* is added to .* size .*",
"INFO org.apache.hadoop.hdfs.server.datanode.fsdataset.impl.FsDatasetImpl: Time to add replicas to map for block pool .* on volume .*",
"INFO org.apache.hadoop.http.HttpServer2: Added filter .* \(class=.*\)",
"INFO org.apache.hadoop.hdfs.server.blockmanagement.BlockManager: The block deletion will start around .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FileJournalManager: Recovering unfinalized segments in .*",
"INFO org.apache.hadoop.util.HostsFileReader: Refreshing hosts \(include/exclude\) list",
"WARN org.apache.hadoop.yarn.server.nodemanager.containermanager.monitor.ContainersMonitorImpl: NodeManager configured with .* physical memory allocated to containers, which is more than .* of the total physical memory available .*\)\. Thrashing might happen\.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.recovery.RMStateStore: Updating AMRMToken",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Initialized root queue root: numChildQueue= .* capacity=.*, absoluteCapacity=.*, usedResources=<memory:.*, vCores:.*>usedCapacity=.*, numApps=.*, numContainers=.*",
"INFO org.apache.hadoop.yarn.webapp.View: Getting list of all Jobs.",
"INFO org.apache.hadoop.util.JvmPauseMonitor: Detected pause in JVM or host machine \(eg GC\): pause of approximately .*",
"INFO org.apache.hadoop.yarn.webapp.WebApps: Web app jobhistory started at .*",
"INFO org.apache.hadoop.hdfs.server.namenode.EditLogFileOutputStream: Nothing to flush",
"WARN org.apache.hadoop.hdfs.server.datanode.DataNode: .* datanodeUuid=.*, infoPort=.*, infoSecurePort=.*, ipcPort=.*, storageInfo=lv=.*\):Failed to transfer .* to .* got",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: opWriteBlock .* received exception org.apache.hadoop.hdfs.server.datanode.ReplicaAlreadyExistsException: Block .* already exists in state .* and thus cannot be created.",
"INFO org.apache.hadoop.mapreduce.v2.hs.HistoryFileManager: Explicitly setting permissions to : .*",
"INFO org.apache.hadoop.hdfs.server.namenode.SecondaryNameNode: Image has not changed. Will not download image.",
"INFO org.apache.hadoop.hdfs.server.datanode.DirectoryScanner: BlockPool .* Total blocks: .* missing metadata .* missing block .* missing blocks in .* mismatched .*",
"WARN org.apache.hadoop.yarn.server.nodemanager.containermanager.ContainerManagerImpl: Event EventType: .* sent to absent container .*",
"INFO org.apache.hadoop.mapreduce.jobhistory.JobSummary: jobId=.* mean",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ResourceManager: Using Scheduler: org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler",
"INFO org.apache.hadoop.hdfs.StateChange: STATE\* Safe mode ON, in safe mode extension.",
"INFO org.apache.hadoop.hdfs.server.namenode.FSImageFormatProtobuf: Image file .* of size .* bytes saved in .* seconds.",
"INFO org.apache.hadoop.yarn.server.resourcemanager.ClientRMService: Allocated new applicationId: .*",
"INFO org.apache.hadoop.hdfs.StateChange: .* fsync: .* for .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.RMAppManager\$ApplicationSummary: appId=.*,name=.*,user=.*,queue=.*,state=.*,trackingUrl=.*,appMasterHost=.*,startTime=.*,finishTime=.*,finalStatus=.*,memorySeconds=.*,vcoreSeconds=.*,preemptedAMContainers=.*,preemptedNonAMContainers=.*,preemptedResources=<memory:.*\\, vCores:.*>,applicationType=.*",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: Starting CheckDiskError Thread",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: .* error processing .* operation src: .* dst: .* org.apache.hadoop.hdfs.server.datanode.ReplicaAlreadyExistsException: Block .* already exists in state .* and thus cannot be created.",
"INFO org.apache.hadoop.hdfs.server.datanode.DataNode: .* datanodeUuid=.*, infoPort=.*, infoSecurePort=.*, ipcPort=.*, storageInfo=lv=.*\) Starting thread to transfer .* to .*",
"INFO org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Need to save fs image\? false \(staleImage=.*, haEnabled=.*, isRollingUpgrade=.*\)",
"WARN org.apache.hadoop.hdfs.server.namenode.FSNamesystem: Only one namespace edits storage directory \(dfs.namenode.edits.dir\) configured. Beware of data loss due to lack of redundant storage directories!",
"INFO org.apache.hadoop.mapreduce.v2.hs.HistoryFileManager: Perms after creating .*, Expected: .*",
"INFO org.apache.hadoop.mapreduce.v2.hs.JobHistory: JobHistory Init",
"INFO org.apache.hadoop.mapreduce.v2.hs.HistoryFileManager: Initializing Existing Jobs\.\.\.",
"INFO org.apache.hadoop.metrics2.impl.MetricsSystemImpl: JobHistoryServer metrics system started",
"INFO org.apache.hadoop.mapreduce.v2.hs.HistoryFileManager: Deleting JobSummary file: \[.*\]",
"INFO org.apache.hadoop.mapreduce.v2.hs.HistoryFileManager: Moving .* to .*",
"INFO org.apache.hadoop.hdfs.server.namenode.SecondaryNameNode: Checkpoint Period :.* secs \(.* min\)",
"INFO org.apache.hadoop.mapreduce.v2.hs.CachedHistoryStorage: CachedHistoryStorage Init",
"INFO org.apache.hadoop.mapreduce.v2.hs.HistoryClientService: Instantiated HistoryClientService at .*/.*:.*",
"INFO BlockStateChange: BLOCK\* ask .* to replicate .* to datanode\(s\) .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.scheduler.capacity.CapacityScheduler: Initialized queue: root: numChildQueue= .* capacity=.*, absoluteCapacity=.*, usedResources=<memory:.*, vCores:.*>usedCapacity=.*, numApps=.*, numContainers=.*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.amlauncher.AMLauncher: Done launching container Container: \[ContainerId: .* NodeId: .*, NodeHttpAddress: .*, Resource: <memory:.*, vCores:.*>, Priority: .* Token: Token { kind: .*, service: .* }, \] for AM .*",
"INFO org.apache.hadoop.yarn.server.resourcemanager.amlauncher.AMLauncher: Command to launch container .*: .*",
"INFO org.apache.hadoop.http.HttpServer2: addJerseyResourcePackage: packageName=.*, pathSpec=.*",
]

#rep_logs = [] 

reuse_filename="CODE32_REUSE.p"

if __name__ == '__main__':

    global debug_mode
    global new_pattern_added
    debug_mode = False

    openfile_list = []
    try:
        parser = argparse.ArgumentParser(description="")
        parser.add_argument('--logfile', type=argparse.FileType('r'), nargs='+', required=True, help='List of one or more input log files')
        parser.add_argument('--debug', action='store_true', required=False, help='It walks through each log processing and print out messages.')
        parser.add_argument('--linear', action='store_true', required=False, help='This option forces it to run in linear execution path without branching off. Suitable for quick results.')

        args = parser.parse_args()
        openfile_list = args.logfile

        debug_mode = False
        if args.debug==None:
            debug_mode = False
        else:
            debug_mode = bool(args.debug)

        linear_mode = False
        if args.linear==None:
            linear_mode = False
        else:
            linear_mode = bool(args.linear)

    except Exception, e:
        print('Error: %s' % str(e))

    print "Loading all logs into memory."
    raw_logs = read_log_files( openfile_list, None ) 


#    cur_node = pickle.load(open("cur_node.bin","rb"))
#    #for aa in cur_node.log_templates[:50]:
#    #    print aa
#    #sys.exit(0)
#    mylogs = copy.deepcopy(raw_logs)
#    SL,CPL = compute_slcpl(mylogs, cur_node.log_templates)
#    sys.exit(0)


    old_log_count = len(raw_logs)
    rep_logs = remove_log_template_matches(raw_logs, prepopulated_log_templates)
    print "Old Log count:",old_log_count
    print "New Log count:",len(raw_logs)
    print "Prepopulated log template count:", len(prepopulated_log_templates)
    print "Representative logs count:", len(rep_logs)
    log_templates = copy.deepcopy(prepopulated_log_templates)
    if len(raw_logs)==0:
        print("No logs to process.")
        sys.exit(0)


    print "Preprocessing logs..."
    all_logs = preprocess_known_patterns(raw_logs) 
                                                   
    print "Sequence number:", seqnum
    print "Tokenizing all logs."
    all_tlogs = do_tokenization(all_logs)
    print "Done tokenizing.", len(all_logs)

    #apply_all_patterns(all_tlogs)
    print "Rearranging numbers of known patterns to make values unique ..."
    uniquify_numbers(all_tlogs)
    print "Done rearranging values."
    print "len(all_tlogs):",len(all_tlogs)

    tree = Tree()
    tree.create_node("TOP", len(raw_logs), "top")
    cur_node = tree.find_inprogress_node()

    runtime_checkpt = time.time()
    #while all_vect.count(-1)>0:
    while -1 in cur_node.all_vect:

        sampled_logs, tokenized_logs = sample_by_token_length_and_space_count(all_logs, all_tlogs, cur_node.all_vect)
        #sampled_logs = random_sample_logs(all_logs, RANDOM_SAMPLE_SIZE)
        #sampled_logs = sample_by_length(all_logs, all_vect, RANDOM_SAMPLE_SIZE)
        #sampled_logs = sample_by_signature(all_logs, RANDOM_SAMPLE_SIZE)
        if len(sampled_logs)==0:
            break
        if debug_mode:
            print "Sampled_logs size:", len(sampled_logs)

        #tokenized_logs = do_tokenization(sampled_logs)
        #sampled_logs,tokenized_logs = sample_by_term_correlation(all_logs, RANDOM_SAMPLE_SIZE)

        if debug_mode:
            print "Number of logs:",len(sampled_logs)
            print "Number of tokenized logs:",len(tokenized_logs)

        #apply_all_patterns(tokenized_logs)
        replace_known_patterns(tokenized_logs)

        candidate_set = forge_candidate_log_templates(tokenized_logs,cur_node.rep_logs) # returns a list of candidate log templates

#        log_template = candidate_set[0]
#        removed_count = mark_matched_logs(all_logs, cur_node.all_vect, cur_node.rep_logs, log_template, len(cur_node.log_templates))
#        print "["+format(len(cur_node.log_templates),'3d')+"]", format(cur_node.all_vect.count(-1),'5d'), format(removed_count,'4d'), log_template
#        cur_node.log_templates.append({"count":removed_count,"template":log_template})
#        #raw_templates.append(postprocess_raw_template(filter_words))

        if (linear_mode):
            tval = 99999
        else:
            tval = 1

        if len(candidate_set)>tval:
            for log_template in candidate_set:
                print "\033[93;100mCandidate:"+"\033[0m \033[90;103m", log_template, "\033[0m "
            
                conflicted = exist_match(log_template, cur_node.rep_logs)
                if conflicted>=0:
                    print "ERROR<1>: log_template overlaps with one of the already-found log templates."
                    sys.exit(0)

                new_name = 'n'+str(tree.serial)
                new_identifier  = 'i'+str(tree.serial)
                tree.create_node( new_name, len(raw_logs), new_identifier, parent=cur_node.identifier)
                tree.serial += 1
    
                new_node = tree.find_node(new_identifier)
                # copy internal state
                new_node.log_templates = copy.deepcopy(cur_node.log_templates)
                new_node.rep_logs = copy.deepcopy(cur_node.rep_logs)
                new_node.all_vect = copy.deepcopy(cur_node.all_vect)

                new_node.print_node()
    
                # Mark matched logs from all_logs using new log template. 
                removed_count = mark_matched_logs(all_logs, new_node.all_vect, new_node.rep_logs, log_template, len(new_node.log_templates))
                if removed_count==0:
                    print "WARNING[1]: no logs removed from the template!"
                    print "TEMPLATE->",log_template
                    print "Remaining logs:", len(all_logs)-sum(1 for x in cur_node.all_vect if x>0)
                    print "Log template count:", len(cur_node.log_templates)
                    for p in standalone_patterns:
                        print p['label'],p['pattern']
                    sys.exit(0)
                    break
                print "["+format(len(new_node.log_templates),'3d')+"]", format(new_node.all_vect.count(-1),'5d'), format(removed_count,'4d'), log_template
                # attach new candidate log templates
                new_node.log_templates.append({"count":removed_count,"template":log_template})
                print "\033[1;34m",tree.show("top"),"...\033[0m"
        else: 

            log_template = candidate_set[0]

            #conflicted = exist_match(log_template, cur_node.rep_logs)
            #if conflicted>=0:
            #    print "\033[32;101m ERROR<2> \033[0m ","\033[1;91mlog_template overlaps with one of the already-found log templates.\033[0m"
            #    print "   \033[4;42mConflicting log:\033[0m \033[1;32m", cur_node.rep_logs[conflicted], "\033[0m"
            #    print "   \033[4;43mExisting log template:\033[0m", cur_node.log_templates[conflicted]['template']
            #    print "   \033[4;43mCurrent log template: \033[0m", log_template
            #    sys.exit(0)

            removed_count = mark_matched_logs(all_logs, cur_node.all_vect, cur_node.rep_logs, log_template, len(cur_node.log_templates))
            if removed_count==0:
                print "WARNING[2]: No logs removed from the template!"
                print "TEMPLATE->",log_template
                print "Remaining logs:", len(all_logs)-sum(1 for x in cur_node.all_vect if x>0)
                print "Log template count:", len(cur_node.log_templates)
                for p in standalone_patterns:
                    print p['label'],p['pattern']
                sys.exit(0)
                break
            #print "["+format(len(cur_node.log_templates),'3d')+"]", format(cur_node.all_vect.count(-1),'5d'), format(removed_count,'4d'), log_template
            print log_template
            cur_node.log_templates.append({"count":removed_count,"template":log_template})
            #raw_templates.append(postprocess_raw_template(filter_words))

        if cur_node.is_leaf_node() and -1 in cur_node.all_vect: 
            pass
        else:
            cur_node = tree.find_inprogress_node()
            if cur_node==None:
                break
            print "\033[1;94m=> Switching to ",cur_node.name,"\033[0m"

        if debug_mode:
            raw_input("\033[1;94m->Press ENTER to continue ...\033[0m")

    # END of outer while loop - done removing all the logs
    runtime_elapsed = time.time() - runtime_checkpt

    print "{0:8.3f}".format(tm001), "Apply all patterns"
    print "{0:8.3f}".format(tm002), "Apply new patterns"
    print "{0:8.3f}".format(tm003), "Random sampling logs"
    print "{0:8.3f}".format(tm006), "Sampling by split token length"
    print "{0:8.3f}".format(tm007), "Sampling by signature made of special characters"
    print "{0:8.3f}".format(tm008), "Sampling by term correlation analysis"
    print "{0:8.3f}".format(tm009), "Sampling by term filtering"
    print "{0:8.3f}".format(tm010), "Replenish term band"
    print "{0:8.3f}".format(tm004), "Tokenizing logs"
    print "{0:8.3f}".format(tm005), "Filtering all columns"
    print "{0:8.3f}".format(tm011), "Match and remove logs"
    print "{0:8.3f}".format(runtime_elapsed-tm001-tm002-tm003-tm004-tm005-tm006-tm007-tm008-tm009-tm010-tm011), "Unaccounted"
    print "{0:8.3f}".format(runtime_elapsed), "\033[0;103mEnd-to-end runtime\033[0m"
