#!/usr/bin/python
#-*- coding: utf-8 -*-

import os
import re
import ast
import sys
import copy
import time
import uuid
import argparse
import pickle
from random import randint
from collections import defaultdict

prepopulated_log_templates = [
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* ColumnFamilyStore.java:899 \- Enqueuing flush of .*: .* \(.*%\) on\-heap, .* \(.*%\) off\-heap",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* Memtable.java:490 \- Completed flushing .* \(.*\) for commitlog position CommitLogPosition\(segmentId=.*, position=.*\)",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* Memtable.java:461 \- Writing .*\(.* serialized bytes, .* ops, .* of on/off\-heap limit\), flushed range = .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* ColumnFamilyStore.java:1197 \- Flushed to \[BigTableReader\(path='.*'\)\] \(.* sstables, .*\), biggest .*, smallest .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CompactionTask.java:155 \- Compacting \(.*\) \[.*\]",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CompactionTask.java:255 \- Compacted \(.*\) .* sstables to \[.*\] to level=.*\. .* to .* \(.*\ of original\) in .*\. Read Throughput = .*, Write Throughput = .*, Row Throughput = .*\. .* total partitions merged to .*. Partition merge counts were {.*}",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* IndexSummaryRedistribution.java:75 \- Redistributing index summaries",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* Schema.java:425 \- Adding .* to cfIdMap",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* ColumnFamilyStore.java:406 \- Initializing .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* SSTableReader.java:506 \- Opening .* \(.*\)",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* ColumnFamilyStore.java:954 \- forceFlush requested but everything is clean in .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* AutoSavingCache.java:394 \- Saved KeyCache \(.* items\) in .* ms",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:2322 \- New node .*/.* at token .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* AbstractCommitLogSegmentManager.java:107 \- No segments in reserve; creating a fresh one",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:484 \- .* memory: init = .*\(.*\) used = .*\(.*\) committed = .*\(.*\) max = .*\(.*\)",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* ViewManager.java:137 \- Not submitting build tasks for views in keyspace .* as storage service is not initialized",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:290 \- opening keyspace .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* YamlConfigurationLoader.java:108 \- Loading settings from file:.*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* AbstractCommitLogSegmentManager.java:329 \- Segment CommitLogSegment\(.*\) is no longer active and will be deleted now",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CommitLog.java:317 \- Commit log segment CommitLogSegment\(.*\) is unused",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StartupChecks.java:280 \- Checking directory .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:101 \- .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:1441 \- DRAINING: starting drain process",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:1441 \- DRAINED",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:1441 \- NORMAL",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:606 \- Token metadata: Normal Tokens:",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:606 \- Token metadata:",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:599 \- Populating token metadata from system tables",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* HintsService.java:220 \- Paused hints dispatch",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* TokenMetadata.java:479 \- Updating topology for localhost/127.0.0.1",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:51 \- .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MigrationManager.java:572 \- Gossiping my schema version .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:2252 \- Node .* state .*, token \[.*\]",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:2255 \- Node .* state jump to .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CommitLogReader.java:223 \- Reading .* \(CL version 6, messaging version 11, compression null\)",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CommitLogReader.java:214 \- Finished reading .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* YamlConfigurationLoader.java:89 \- Configuration location: file:.*",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StartupChecks.java:197 \- OpenJDK is not recommended. Please upgrade to the newest Oracle Java release",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StartupChecks.java:160 \- JMX is not enabled to receive remote connections. Please see cassandra\-env.sh for more info\.",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* SigarLibrary.java:44 \- Initializing SIGAR library",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* SigarLibrary.java:174 \- Cassandra server running in degraded mode\. Is swap disabled\? : .*, Address space adequate\? : .*, nofile limit adequate\? : .*, nproc limit adequate\? : .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* RateBasedBackPressure.java:123 \- Initialized back\-pressure with high ratio: .*, factor: .*, flow: .*, window size: .*\.",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* QueryProcessor.java:115 \- Initialized prepared statement caches with .* \(.*\) and .* \(Thrift\)",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* NativeLibrary.java:174 \- JNA mlockall successful",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* JMXServerUtils.java:249 \- Configured JMX server at: service:jmx:rmi:.*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* DatabaseDescriptor.java:710 \- Back\-pressure is disabled with strategy org.apache.cassandra.net.RateBasedBackPressure{high_ratio=.*, factor=.*, flow=.*}.",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* DatabaseDescriptor.java:424 \- Global memtable off\-heap threshold is enabled at .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* DatabaseDescriptor.java:420 \- Global memtable on\-heap threshold is enabled at .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* DatabaseDescriptor.java:366 \- DiskAccessMode '.*' determined to be mmap, indexAccessMode is mmap",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* DatabaseDescriptor.java:358 \- Syncing log with a period of .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* Config.java:481 \- Node configuration:\[.*\]",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:488 \- JVM Arguments: \[.*\]",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:486 \- Classpath: .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:479 \- Heap size: .*/.*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:478 \- JVM vendor/version: OpenJDK 64-Bit Server VM/1.8.0_131",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:471 \- Hostname: .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CacheService.java:174 \- Scheduling counter cache save to every .* seconds \(going to save all keys\).",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CacheService.java:163 \- Initializing counter cache with capacity of .* .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CacheService.java:134 \- Initializing row cache with capacity of .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CacheService.java:112 \- Initializing key cache with capacity of .*\.",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* BufferPool.java:230 \- Global buffer pool is enabled, when pool is exhausted \(max is .*\) it will allocate on heap",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* ApproximateTime.java:44 \- Scheduling approximate time-check task with a precision of .* milliseconds",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:849 \- Bootstrap variables: .* .* .* .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:819 \- Starting up server gossip",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:706 \- Loading persisted ring state",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:621 \- Native protocol supported versions: 3/v3, 4/v4, 5/v5-beta \(default: 4/v4\)",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:619 \- CQL supported versions: 3.4.4 \(default: 3.4.4\)",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:618 \- Thrift API version: 20.1.0",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:617 \- Cassandra version: 3.11.0",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:549 \- Starting shadow gossip round to check for endpoint collision",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* QueryProcessor.java:162 \- Preloaded .* prepared statements",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MessagingService.java:753 \- Starting Messaging Service on .* \(lo\)",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* IndexSummaryManager.java:85 \- Initializing index summary manager with a memory pool size of .* and a resize interval of .* minutes",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* AutoSavingCache.java:173 \- Completed loading \(.*;.*\) KeyCache cache",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StartupChecks.java:127 \- jemalloc shared library could not be preloaded to speed up memory allocations",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MessagingService.java:984 \- Waiting for messaging service to quiesce",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MessagingService.java:1338 \- MessagingService has terminated the accept\(\) thread",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* Gossiper.java:1530 \- Announcing shutdown",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CommitLogReplayer.java:129 \- Global replay position is CommitLogPosition\(segmentId=-1, position=0\) from columnfamilies .*{.*}",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CommitLog.java:159 \- Log replay complete, .* replayed mutations",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CommitLog.java:157 \- Replaying .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* AutoSavingCache.java:197 \- reading saved cache .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* OutboundTcpConnection.java:545 \- Unable to connect to .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* SecondaryIndexManager.java:508 \- Executing pre\-join tasks for: CFS\(Keyspace='.*', ColumnFamily='.*'\)",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:231 \- Setting tokens to \[.*\]",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:1439 \- JOINING: Finish joining ring",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* Server.java:156 \- Starting listening for CQL clients on .*:.* \(unencrypted\)\.\.\.",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* Server.java:155 \- Using Netty Version: \[.*\]",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* NativeTransportService.java:70 \- Netty using native Epoll event loop",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:527 \- Not starting RPC server as requested. Use JMX \(StorageService->startRPCServer\(\)\) or nodetool \(enablethrift\) to start it",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:992 \- Using saved tokens \[.*\]",
"ERROR \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:706 \- Fatal configuration error",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraDaemon.java:408 \- Completed submission of build tasks for any materialized views defined at startup",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* Server.java:176 \- Stop listening for CQL clients",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* GCInspector.java:284 - ParNew GC in .*ms. CMS Old Gen: .* -> .*; Par Eden Space: .* -> .*; Par Survivor Space: .* -> .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* OutputHandler.java:42 \- Verifying BigTableReader\(path='.*'\) \(.*\)",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* OutputHandler.java:42 \- Checking computed hash of BigTableReader\(path='.*'\)",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MigrationManager.java:310 \- Create new Keyspace: KeyspaceMetadata{.*}",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StartupChecks.java:286 \- Directory .* doesn't exist",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StartupChecks.java:131 \- jemalloc seems to be preloaded from .*",
# https://github.com/hobinyoon/apache-cassandra-3.0.5-src/blob/master/src/java/org/apache/cassandra/service/MigrationManager.java:413
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MigrationManager.java:426 \- Update table '.*' From .*\[.*\] To .*\[.*\]",
# logger.debug("application result is {}", this);
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CFMetaData.java:793 \- application result is .*\[.*\]",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CFMetaData.java:768 \- applying .*\[.*\] to .*\[.*\]",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MigrationManager.java:355 \- Create new table: .*\[.*\]",
#- logger.trace("ScheduledThreadPoolExecutor has shut down as part of C* shutdown");
#+ logger.debug("ScheduledThreadPoolExecutor has shut down as part of C* shutdown");
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* DebuggableScheduledThreadPoolExecutor.java:57 \- ScheduledThreadPoolExecutor has shut down as part of C\* shutdown",
# 예) INFO  [Service Thread] 2017-12-25 18:29:56,343 StatusLogger.java:98 - Table                       Memtable ops,data
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:98 \- .*",
# INFO  [Service Thread] 2017-12-25 18:29:56,343 StatusLogger.java:91 - RowCache                          0                        0                      all
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:91 \- .*",
# INFO  [Service Thread] 2017-12-25 18:29:56,342 StatusLogger.java:85 - KeyCache                      42704                101711872                      all
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:85 \- .*",
# INFO  [Service Thread] 2017-12-25 18:29:56,342 StatusLogger.java:83 - Cache Type                     Size                 Capacity               KeysToSave
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:83 \- .*",
# INFO  [Service Thread] 2017-12-25 18:29:56,342 StatusLogger.java:73 - MessagingService                n/a       0/0
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:73 \- .*",
# INFO  [Service Thread] 2017-12-25 18:29:56,341 StatusLogger.java:61 - CompactionManager                 0         0
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:61 \- .*",
# INFO  [Service Thread] 2017-12-25 18:29:55,652 StatusLogger.java:47 - Pool Name                    Active   Pending      Completed   Blocked  All Time Blocked
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:47 \- .*",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* GCInspector.java:282 - ParNew GC in .*ms. CMS Old Gen: .* -> .*; Par Eden Space: .* -> .*; Par Survivor Space: .* -> .*",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* SystemKeyspace.java:1083 \- No host ID found, created .* \(Note: This should happen exactly once per node\).",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:977 \- Generated random tokens. tokens are \[.*\]",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageService.java:857 \- This node will not auto bootstrap because it is configured to be a seed node.",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* OutboundTcpConnection.java:424 \- Attempting to connect to .*",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* NoSpamLogger.java:94 \- Out of .* commit log syncs over the past .* with average duration of .*, .* have exceeded the configured commit interval by an average of .*",
# https://github.com/PaytmLabs/cassandra/blob/master/src/java/org/apache/cassandra/service/MigrationManager.java
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MigrationManager.java:406 \- Update Keyspace '.*' From .* To .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CommitLog.java:152 \- No commitlog files found; skipping replay",
# https://fossies.org/linux/apache-cassandra/src/java/org/apache/cassandra/auth/CassandraRoleManager.java
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraRoleManager.java:355 \- Created default superuser role .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraIndex.java:719 \- Index build of .* complete",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraIndex.java:707 \- Submitting index build of .* for data in .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* ReadCallback.java:132 \- Timed out; received .* of .* responses",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* OutboundTcpConnection.java:560 \- Handshaking version with .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* OutboundTcpConnection.java:532 \- Done connecting to .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* NoSpamLogger.java:91 \- Some operations were slow, details available at debug level \(debug.log\)",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MonitoringTask.java:93 \- Scheduling monitoring task with report interval of .* ms, max operations .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MonitoringTask.java:173 \- .* operations were slow in the last .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MigrationManager.java:478 \- Drop table '.*'",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* MigrationManager.java:462 \- Drop Keyspace '.*'",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraRoleManager.java:399 \- Setup task failed with error, rescheduling",
"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CassandraRoleManager.java:360 \- CassandraRoleManager skipped default role setup: some nodes were not ready",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CompactionManager.java:300 \- No sstables to .* for .*\..*",
# https://www.javatips.net/api/cassandra-master/cassandra-trunk/src/java/org/apache/cassandra/db/compaction/CompactionManager.java
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* CompactionManager.java:1813 \- Executor has been shut down, not submitting .*",
"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* OutboundTcpConnection.java:108 \- OutboundTcpConnection using coalescing strategy .*",
"DEBUG \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StorageProxy.java:2361 \- Schemas are in agreement.",
## https://github.com/hobinyoon/apache-cassandra-3.0.5-src/blob/fc6710f9ce117e22286b9f42955b3e7632844160/src/java/org/apache/cassandra/utils/StatusLogger.java
##"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* StatusLogger.java:.* \- .*",
##"INFO \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* GCInspector.java:284 \- ParNew GC in .*\. CMS Old Gen: .* \-> .*; Par Eden Space: .* \-> .*; Par Survivor Space: .* \-> .*",
##"WARN \[.*\] \d{4}-\d\d\-\d\d \d\d:\d\d:\d\d,\d* GCInspector.java:282 \- ParNew GC in .*\. CMS Old Gen: .* \-> .*; Par Eden Space: .* \-> .*; Par Survivor Space: .* \-> .*",
]

def remove_log_template_matches(logs, logtem):
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

        print "\033[0;34m"+"["+format(i,'3d')+"]",format(len(logs),'5d'),format(removed_logs,'4d'),"\033[0m","\033[0;35m\""+logtem[i]+"\",\033[0m"


def read_a_logfile(infile):
    logs = []

    lines = [line.rstrip() for line in infile]
    for log in lines:
        if len(log.strip())==0:
            continue
        log = " ".join(log.split())
        logs.append(log)

    print "Total number of logs loaded:",len(logs)
    return logs


if __name__ == '__main__':

    global new_pattern_added

    try:
        parser = argparse.ArgumentParser(description="")
        parser.add_argument('--logfile', type=argparse.FileType('r'), required=True, help='One logfile name.')

        args = parser.parse_args()
        input_logfile = args.logfile

    except Exception, e:
        print('Error: %s' % str(e))

    print "Loading all logs into memory."
    raw_logs = read_a_logfile(input_logfile)

    loaded_log_count = len(raw_logs)
    remove_log_template_matches(raw_logs, prepopulated_log_templates)
    print "Loaded Log count:", loaded_log_count
    print "Remaining Log count:",len(raw_logs)
    print "Prepopulated log template count:", len(prepopulated_log_templates)

    for i in range(0,len(raw_logs)):
        print raw_logs[i]

    log_templates = []
    for t in prepopulated_log_templates:
        log_templates.append({"count":1,"template":t}) 
    pickle.dump(log_templates,open("TMP1029348.bin","wb"))

