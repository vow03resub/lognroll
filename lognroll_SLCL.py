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
import logutils
from scipy import stats
from random import randint
from collections import defaultdict

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
LOGLEN_THRESHOLD=(1024*5)
seqnum=1

MOD_FACTOR = 32

MIN_NUM_OF_TERMS=1 # minimum required number of terms in the term band, initially I introduced this to filter out one word term band.

RANDOM_SAMPLE_SIZE=4096
# frequent words must be above this cut-off line
CUTOFF_COUNT=20
CC_THRESHOLD=0.1 # select only the word-term pairs that have correlation above this value

WILDCARD_THRESHOLD=0.9 # if almost all filter vectors are filled (above this threshold), then the rest of them will be * regardless of value distribution.
# TODO below not used for now
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

    # if input string is just one delimiter character, just return it
    if str_data==delim:
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

    # process each token
    for i in reversed(range(0,len(tokenized))):

        tok = tokenized[i]
        if tok==delim:
            continue

        # if any delimiter in lower priority than current one exists, recursively call the split_by_delimiter 
        low_delim_found = False
        cur = delimiter_list.index(delim)
        for j in range(cur+1,len(delimiter_list)):
            if delimiter_list[j] in tok:

                tokenized = tokenized[:i] + split_by_delimiter(lvl+4, tok, delimiter_list[j]) + tokenized[i+1:]
                low_delim_found = True
                break
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

    # if there are less than 3 values, it is too less to say whether there is a pattern or not
    if len(klist)<PATTERN_THRESHOLD:
        return False

    # check for string length
    slen = len(klist[0])
    for w in klist[1:]:
        if slen!=len(w):
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

    smask = re.sub("[0-9]","\d",smask)

    if fixed_count==0:
        return False

    pat = { "pattern": smask, "label":"~"+str(400+len(discovered_patterns))+"~"}

    pattern_exists = False
    for p in discovered_patterns:
        if p["pattern"]==smask:
            pattern_exists = True
            pat = p
            break
    if not pattern_exists:
        discovered_patterns.append(pat)

    return True


def all_terms_exist(s, keywords_list):
    for term in keywords_list:
        #if term not in custom_split(s):
        if term not in s:
            return False
    return True


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

]

def preprocess_known_patterns(logs):
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


def exist_match(log_template, logs):
    for i in range(0,len(logs)):
        matched = re.match("^"+log_template+"$", logs[i])
        if matched!=None:
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
 
def finalize_filter_with_star(fword,fmask):
    tfword = copy.deepcopy(fword)
    tfmask = copy.deepcopy(fmask)
    for i in range(0,len(tfmask)):
        if tfmask[i]==0:
            tfword[i]="*"
            tfmask[i]=1
    return tfword, tfmask

def mark_matched_logs2(logs, mask, template):
    marked=0
    for i in range(0, len(logs)):
        if mask[i]>-1:
            continue

        log = logs[i]
        if len(log)>LOGLEN_THRESHOLD:
            log = log[:LOGLEN_THRESHOLD]

        if len(log)==LOGLEN_THRESHOLD:
            tttt=template.replace(".*","\S*")
            template=tttt

        matched = re.match('^'+template+'$', log)
        if matched!=None:
            mask[i]=99
            marked+=1
    return marked


log_templates = []

prepopulated_log_templates = [ ]

if __name__ == '__main__':

    openfile_list = []
    try:
        parser = argparse.ArgumentParser(description="")
        parser.add_argument('--logname', required=True, help='One of hadoop, openstack, and cassandra')
        args = parser.parse_args()
        #openfile_list = args.logfile
        log_name=args.logname

    except Exception, e:
        print('Error: %s' % str(e))

    if log_name=="hadoop":
        openfile_list=[open("log_sample/hadoop_clean.log","r")]
        reuse_file="./CODE32_REUSE_hadoop_lognroll.p"
    elif log_name=="openstack":
        openfile_list=[open("log_sample/openstack_small_clean.log","r")]
        reuse_file="./CODE32_REUSE_openstack_lognroll.p"
    elif log_name=="cassandra":
        openfile_list=[open("log_sample/cassandra_debug.log","r")]
        reuse_file="./CODE32_REUSE_cassandra_lognroll.p"
    else:
        print "Unrecognized log name:", log_name
        sys.exit(0)

    xlog_templates = pickle.load(open(reuse_file,"rb"))

    print "Loading all logs into memory."
    raw_logs = read_log_files( openfile_list, None ) 
    mylogs = copy.deepcopy(raw_logs)

    sum_matched=0
    for t in sorted(xlog_templates,reverse=True):
        sum_matched += int(t[0])
    print "    Sum of matched logs:",sum_matched

    selected=[]
    for t in sorted(xlog_templates, reverse=True):
        # first, build a list of index to delete
        to_delete = []
        for i in range(0,len(mylogs)):
            log = mylogs[i]
            matched = re.match("^"+t[1]+"$",log)

            if matched!=None:
                to_delete.append(i)

        # delete matched logs
        before_removal = len(mylogs)
        to_delete = sorted(to_delete)
        for i in reversed(sorted(to_delete)):
            del mylogs[i]
        del_count = before_removal - len(mylogs)
        #print "Removed",del_count,"logs.",t[1]
        #if del_count>0:
        #    print "\033[33;31m",format(del_count,'5d'), format(int(t[0]),'5d'), t[1], "\033[0m"

        if del_count>0:
            #selected.append(t[1])
            selected.append({"count":del_count,"template":str(t[1])})
        #if len(selected)==1:
        #    break
    
    print "   ",len(mylogs),"logs remaining."
    print "    Initial template count:", len(xlog_templates)
    print "    Selected template count:", len(selected)

    mylogs2 = copy.deepcopy(raw_logs)
    SL,CPL = logutils.compute_slcpl(mylogs2, selected)

    print "    SL= "+str(SL)
    print "    CPL= "+str(CPL)
    #print "    \033[1;94mscore= "+str(SL*1.0/(1.0+CPL)), "\033[0m"

    if log_name=="hadoop":
        log_name_color= "91;106m"
    elif log_name=="openstack":
        log_name_color= "91;103m"
    elif log_name=="cassandra":
        log_name_color= "1;106m"
    else:
        log_name_color= "32;31m"

    mname_color= "1;94m"

    print "    \033[1;94mscore= "+str(SL-CPL), "\033[0m", "\033["+log_name_color+log_name+"\033[0m", "\033["+mname_color+"Lognroll\033[0m"

    print sum_matched, len(mylogs), len(xlog_templates), len(selected), SL, CPL
