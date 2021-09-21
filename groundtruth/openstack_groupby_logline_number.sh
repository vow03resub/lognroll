#!/bin/bash

#cat openstack_small.log | awk '{print $(NF)}' | sort | uniq -c | sort -nr


LLN="/usr/local/lib/python2.7/dist-packages/ovs/poller.py:246 /opt/stack/neutron/neutron/plugins/ml2/drivers/openvswitch/agent/ovs_neutron_agent.py:1949 /opt/stack/neutron/neutron/plugins/ml2/drivers/openvswitch/agent/ovs_neutron_agent.py:1769 /opt/stack/neutron/neutron/plugins/ml2/drivers/openvswitch/agent/openflow/native/ofswitch.py:96 /usr/local/lib/python2.7/dist-packages/oslo_service/periodic_task.py:215 /usr/local/lib/python2.7/dist-packages/neutron_lib/callbacks/manager.py:167 /usr/local/lib/python2.7/dist-packages/keystonemiddleware/auth_token/__init__.py:401 /opt/stack/nova/nova/api/openstack/placement/requestlog.py:38 /opt/stack/keystone/keystone/policy/backends/rules.py:33 /opt/stack/keystone/keystone/middleware/auth.py:239 /opt/stack/keystone/keystone/common/authorization.py:240 /opt/stack/keystone/keystone/common/authorization.py:136 /usr/local/lib/python2.7/dist-packages/oslo_concurrency/lockutils.py:282 /usr/local/lib/python2.7/dist-packages/oslo_concurrency/lockutils.py:270 /usr/local/lib/python2.7/dist-packages/oslo_concurrency/processutils.py:385 /usr/local/lib/python2.7/dist-packages/oslo_concurrency/processutils.py:355 /opt/stack/nova/nova/scheduler/client/report.py:438 /opt/stack/neutron/neutron/db/agents_db.py:281 /opt/stack/cinder/cinder/volume/drivers/lvm.py:204 /opt/stack/cinder/cinder/scheduler/host_manager.py:510 /opt/stack/cinder/cinder/manager.py:193 /opt/stack/nova/nova/compute/resource_tracker.py:763 /opt/stack/nova/nova/compute/resource_tracker.py:747 /opt/stack/nova/nova/compute/resource_tracker.py:730 /opt/stack/nova/nova/compute/resource_tracker.py:701 /opt/stack/nova/nova/compute/resource_tracker.py:609 /opt/stack/nova/nova/compute/manager.py:6560 /opt/stack/nova/nova/compute/manager.py:5964 /opt/stack/nova/nova/compute/manager.py:5892 /opt/stack/nova/nova/compute/manager.py:5888 /opt/stack/keystone/keystone/common/fernet_utils.py:307 /opt/stack/nova/nova/compute/manager.py:6939 /opt/stack/nova/nova/compute/manager.py:6930 /opt/stack/nova/nova/compute/manager.py:6973 /opt/stack/keystone/keystone/middleware/auth.py:203 /opt/stack/keystone/keystone/auth/core.py:388 /usr/local/lib/python2.7/dist-packages/ovs/reconnect.py:468 /usr/local/lib/python2.7/dist-packages/ovs/poller.py:231 /opt/stack/nova/nova/virt/libvirt/imagecache.py:406 /usr/local/lib/python2.7/dist-packages/ovs/reconnect.py:105 /usr/local/lib/python2.7/dist-packages/eventlet/wsgi.py:883"

for lln in $LLN
do
    clear
    #echo $lln
    cat ../log_sample/openstack_small_clean.log | grep $lln | head -50
    read -p "Press Enter to continue" Status
done

