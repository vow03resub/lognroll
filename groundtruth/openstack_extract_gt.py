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
# https://github.com/openvswitch/ovs
"DEBUG ovsdbapp.backend.ovs_idl.vlog \[.*\] \[.*\] on fd .* {{\(pid=.*\) __log_wakeup /usr/local/lib/python2.7/dist\-packages/ovs/poller.py:246}}",
"DEBUG neutron.plugins.ml2.drivers.openvswitch.agent.ovs_neutron_agent \[.*\] Agent rpc_loop \- iteration:.* started {{\(pid=.*\) rpc_loop /opt/stack/neutron/neutron/plugins/ml2/drivers/openvswitch/agent/ovs_neutron_agent.py:1949}}",
"DEBUG neutron.plugins.ml2.drivers.openvswitch.agent.ovs_neutron_agent \[.*\] Agent rpc_loop \- iteration:.* completed. Processed ports statistics: {'regular': {'updated': .*, 'added': .*, 'removed': .*}}. Elapsed:.* {{\(pid=.*\) loop_count_and_wait /opt/stack/neutron/neutron/plugins/ml2/drivers/openvswitch/agent/ovs_neutron_agent.py:1769}}",
"DEBUG neutron.plugins.ml2.drivers.openvswitch.agent.openflow.native.ofswitch \[.*\] ofctl request version=.*,msg_type=.*,msg_len=.*,xid=.*,OFPFlowStatsRequest\(cookie=.*,cookie_mask=.*,flags=.*,match=.*,out_group=.*,out_port=.*,table_id=.*,type=.*\) result \[OFPFlowStatsReply\(body=\[OFPFlowStats\(byte_count=.*,cookie=.*,duration_nsec=.*,duration_sec=.*,flags=.*,hard_timeout=.*,idle_timeout=.*,instructions=\[\],length=.*,match=.*,packet_count=.*,priority=.*,table_id=.*\)\],flags=.*,type=.*\)\] {{\(pid=.*\) _send_msg /opt/stack/neutron/neutron/plugins/ml2/drivers/openvswitch/agent/openflow/native/ofswitch.py:96}}",
"DEBUG oslo_service.periodic_task \[.*\] Running periodic task .* {{\(pid=.*\) run_periodic_tasks /usr/local/lib/python2.7/dist-packages/oslo_service/periodic_task.py:215}}",
"DEBUG neutron_lib.callbacks.manager \[.*\] Notify callbacks .* for agent, after_update {{\(pid=.*\) _notify_loop /usr/local/lib/python2.7/dist-packages/neutron_lib/callbacks/manager.py:167}}",
"DEBUG keystone.middleware.auth \[.*\] Authenticating user token {{\(pid=.*\) process_request /usr/local/lib/python2.7/dist-packages/keystonemiddleware/auth_token/__init__.py:401}}",
"DEBUG nova.api.openstack.placement.requestlog \[.*\] Starting request: .* \".*\" {{\(pid=.*\) __call__ /opt/stack/nova/nova/api/openstack/placement/requestlog.py:38}}",
"DEBUG keystone.policy.backends.rules \[.*\] enforce identity:validate_token: {'is_delegated_auth': .*, 'access_token_id': .*, 'user_id': .*, 'roles': .*, 'user_domain_id': .*, 'consumer_id': .*, 'trustee_id': .*, 'is_domain': .*, 'is_admin_project': .*, 'trustor_id': .*, 'token': .*, 'project_id': .*, 'trust_id': .*, 'project_domain_id': .*} {{\(pid=.*\) enforce /opt/stack/keystone/keystone/policy/backends/rules.py:33}}",
"DEBUG keystone.middleware.auth \[.*\] RBAC: auth_context: {'is_delegated_auth': .*, 'access_token_id': .*, 'user_id': .*, 'roles': .*, 'user_domain_id': .*, 'consumer_id': .*, 'trustee_id': .*, 'is_domain': .*, 'is_admin_project': .*, 'trustor_id': .*, 'token': .*, 'project_id': .*, 'trust_id': .*, 'project_domain_id': .*} {{\(pid=.*\) fill_context /opt/stack/keystone/keystone/middleware/auth.py:239}}",
"DEBUG keystone.common.authorization \[.*\] RBAC: Authorization granted {{\(pid=.*\) check_policy /opt/stack/keystone/keystone/common/authorization.py:240}}",
"DEBUG keystone.common.authorization \[.*\] RBAC: Authorizing identity:validate_token\(\) {{\(pid=.*\) _build_policy_check_credentials /opt/stack/keystone/keystone/common/authorization.py:136}}",
"DEBUG oslo_concurrency.lockutils \[.*\] Lock .* released by .* :: held .* {{\(pid=.*\) inner /usr/local/lib/python2.7/dist-packages/oslo_concurrency/lockutils.py:282}}",
"DEBUG oslo_concurrency.lockutils \[.*\] Lock .* acquired by .* :: waited .* {{\(pid=.*\) inner /usr/local/lib/python2.7/dist-packages/oslo_concurrency/lockutils.py:270}}",
"DEBUG oslo_concurrency.processutils \[.*\] CMD \".*\" returned: .* in .*s {{\(pid=.*\) execute /usr/local/lib/python2.7/dist-packages/oslo_concurrency/processutils.py:385}}",
"DEBUG oslo_concurrency.processutils \[.*\] Running cmd \(subprocess\): sudo cinder-rootwrap /etc/cinder/rootwrap.conf env LC_ALL=C .* {{\(pid=.*\) execute /usr/local/lib/python2.7/dist\-packages/oslo_concurrency/processutils.py:355}}",
"DEBUG nova.scheduler.client.report \[.*\] Refreshing aggregate associations for resource provider .* {{\(pid=.*\) _ensure_resource_provider /opt/stack/nova/nova/scheduler/client/report.py:438}}",
"DEBUG neutron.db.agents_db \[.*\] Agent healthcheck: found .* active agents {{\(pid=.*\) agent_health_check /opt/stack/neutron/neutron/db/agents_db.py:281}}",
"DEBUG cinder.volume.drivers.lvm \[.*\] Updating volume stats {{\(pid=.*\) _update_volume_stats /opt/stack/cinder/cinder/volume/drivers/lvm.py:204}}",
"DEBUG cinder.scheduler.host_manager \[.*\] Received volume service update from .*: {u'filter_function': .*, u'goodness_function': .*, u'volume_backend_name': .*, u'driver_version': .*, u'sparse_copy_volume': .*, u'pools': \[{u'pool_name': .*, u'filter_function': .*, u'goodness_function': .*, u'total_volumes': .*, u'multiattach': .*, u'provisioned_capacity_gb': .*, u'allocated_capacity_gb': .*, u'thin_provisioning_support': .*, u'free_capacity_gb': .*, u'location_info': .*, u'total_capacity_gb': .*, u'thick_provisioning_support': .*, u'reserved_percentage': .*, u'QoS_support': .*, u'max_over_subscription_ratio': .*}\], u'vendor_name': .*, u'storage_protocol': .*} {{\(pid=.*\) update_service_capabilities /opt/stack/cinder/cinder/scheduler/host_manager.py:510}}",
"DEBUG cinder.manager \[.*\] Notifying Schedulers of capabilities \.\.\. {{\(pid=.*\) _publish_service_capabilities /opt/stack/cinder/cinder/manager.py:193}}",
"DEBUG nova.compute.resource_tracker \[.*\] Total usable vcpus: .*, total allocated vcpus: .* {{\(pid=.*\) _report_final_resource_view /opt/stack/nova/nova/compute/resource_tracker.py:763}}",
"DEBUG nova.compute.resource_tracker \[.*\] Hypervisor: free VCPUs: .* {{\(pid=.*\) _report_hypervisor_resource_view /opt/stack/nova/nova/compute/resource_tracker.py:730}}",
"DEBUG nova.compute.resource_tracker \[.*\] Compute_service record updated for .*:.* {{\(pid=.*\) _update_available_resource /opt/stack/nova/nova/compute/resource_tracker.py:701}}",
"DEBUG nova.compute.resource_tracker \[.*\] Auditing locally available compute resources for .* \(node: .*\) {{\(pid=.*\) update_available_resource /opt/stack/nova/nova/compute/resource_tracker.py:609}}",
"DEBUG nova.compute.manager \[.*\] CONF.reclaim_instance_interval <= 0, skipping\.\.\. {{\(pid=.*\) _reclaim_queued_deletes /opt/stack/nova/nova/compute/manager.py:6560}}",
"DEBUG nova.compute.manager \[.*\] Didn't find any instances for network info cache update. {{\(pid=.*\) _heal_instance_info_cache /opt/stack/nova/nova/compute/manager.py:5964}}",
"DEBUG nova.compute.manager \[.*\] Rebuilding the list of instances to heal {{\(pid=.*\) _heal_instance_info_cache /opt/stack/nova/nova/compute/manager.py:5892}}",
"DEBUG nova.compute.manager \[.*\] Starting heal instance info cache {{\(pid=.*\) _heal_instance_info_cache /opt/stack/nova/nova/compute/manager.py:5888}}",
# https://github.com/openstack/keystone/blob/master/keystone/common/fernet_utils.py
"DEBUG keystone.common.fernet_utils \[.*\] Loaded .* Fernet keys from .*, but \`\[fernet_tokens\] max_active_keys = .*\`; perhaps there have not been enough key rotations to reach \`max_active_keys\` yet\? {{\(pid=.*\) load_keys /opt/stack/keystone/keystone/common/fernet_utils.py:307}}",
"DEBUG nova.compute.manager \[.*\] There are .* instances to clean {{\(pid=.*\) _run_pending_deletes /opt/stack/nova/nova/compute/manager.py:6939}}",
"DEBUG nova.compute.manager \[.*\] Cleaning up deleted instances {{\(pid=.*\) _run_pending_deletes /opt/stack/nova/nova/compute/manager.py:6930}}",
"DEBUG nova.compute.manager \[.*\] Cleaning up deleted instances with incomplete migration {{\(pid=.*\) _cleanup_incomplete_migrations /opt/stack/nova/nova/compute/manager.py:6973}}",
"DEBUG keystone.middleware.auth \[.*\] There is either no auth token in the request or the certificate issuer is not trusted\. No auth context will be set\. {{\(pid=.*\) fill_context /opt/stack/keystone/keystone/middleware/auth.py:203}}",
"DEBUG keystone.auth.core \[.*\] MFA Rules not processed for user \`.*\`. Rule list: \`\[.*\]\` \(Enabled: \`.*\`\)\. {{\(pid=.*\) check_auth_methods_against_rules /opt/stack/keystone/keystone/auth/core.py:388}}",
"DEBUG ovsdbapp.backend.ovs_idl.vlog \[.*\] tcp:.*:.*: entering .* {{\(pid=.*\) _transition /usr/local/lib/python2.7/dist-packages/ovs/reconnect.py:468}}",
"DEBUG ovsdbapp.backend.ovs_idl.vlog \[.*\] .*\-ms timeout {{\(pid=.*\) __log_wakeup /usr/local/lib/python2.7/dist-packages/ovs/poller.py:231}}",
"DEBUG nova.virt.libvirt.imagecache \[.*\] Skipping verification, no base directory at .* {{\(pid=.*\) _get_base /opt/stack/nova/nova/virt/libvirt/imagecache.py:406}}",
"DEBUG ovsdbapp.backend.ovs_idl.vlog \[.*\] tcp:.*:.*: idle .* ms, sending inactivity probe {{\(pid=.*\) run /usr/local/lib/python2.7/dist-packages/ovs/reconnect.py:105}}",
# https://github.com/openstack-archive/deb-python-eventlet/blob/master/eventlet/wsgi.py
"DEBUG eventlet.wsgi.server \[.*\] \(.*\) accepted \('.*', .*\) {{\(pid=.*\) server /usr/local/lib/python2.7/dist-packages/eventlet/wsgi.py:883}}",

"INFO nova.api.openstack.placement.requestlog \[.*\] .* \"GET .*\" status: .* len: .* microversion: .*",
"INFO keystone.common.wsgi \[.*\] GET http://.*",
"INFO nova.scheduler.host_manager \[.*\] Successfully synced instances from host .*\.",
"INFO keystone.common.wsgi \[.*\] POST http://.*",
"INFO nova.compute.resource_tracker \[.*\] Final resource view: name=.* phys_ram=.* used_ram=.* phys_disk=.* used_disk=.* total_vcpus=.* used_vcpus=.* pci_stats=\[.*\]",
"DEBUG nova.compute.resource_tracker \[.*\] Hypervisor/Node resource view: name=.* free_ram=.* free_disk=.* free_vcpus=.* pci_devices=\[{\"dev_id\": .* \"product_id\": .* \"dev_type\": .*, \"numa_node\": .*, \"vendor_id\": .*, \"label\": .*, \"address\": .*}, {\"dev_id\": .*, \"product_id\": .*, \"dev_type\": .*, \"numa_node\": .*, \"vendor_id\": .*, \"label\": .*, \"address\": .*}, {\"dev_id\": .*, \"product_id\": .*, \"dev_type\": .*, \"numa_node\": .*, \"vendor_id\": .*, \"label\": .*, \"address\": .*}, {\"dev_id\": .*, \"product_id\": .*, \"dev_type\": .*, \"numa_node\": .*, \"vendor_id\": .*, \"label\": .*, \"address\": .*",
"INFO eventlet.wsgi.server \[.*\] .* \- \- \[.*/.*/.* .*:.*:.*\] \".*\" .* .* .*",

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
