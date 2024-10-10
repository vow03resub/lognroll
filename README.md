# Lognroll
**source paper**: https://dl.acm.org/doi/abs/10.1145/3464298.3493400

Lognroll is a log parser that generates a log template through iterative filtering. 
It continuously learns one rule that occupies the largest part of the log data set, and generates a log template through iterative filtering that removes the logs targeted by the rule. 


### Running environment

python 3.7+
regex 2022.3.2 (or previous version)
numpy
scipy
argparse
pickle



### How to run
$ python ./lognroll_actual.py --linear --logfile logs/hadoop_clean.log

### Execution Parameter
* --linear: Whether to follow linear execution path along the tree or not.
* --clean: When specified, it deletes intermediate pickle files of tokenized log data and reprocess them. It takes longer.
* --debug: When specified, it walks through each log processing and print out messages.
* --logfile: List of one or more input log files

### Default log files
By default, the logs folder has log files extracted from Hadoop and Cassandra, respectively. If you want to input the log files you extracted yourself, you'd better put the timestamp-removed log file.

