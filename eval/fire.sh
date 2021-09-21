#!/bin/bash

./go.sh abnormal_cinder_volume.1.log 1
./go.sh abnormal_large_instance.1.log 1
./go.sh abnormal_large_instance.2.log 1
./go.sh abnormal_large_instance.3.log 1
./go.sh abnormal_neutron-dhcp_down.1.log 1
./go.sh abnormal_neutron-dhcp_down.2.log 1
./go.sh abnormal_nova-compute_down.1.log 1
./go.sh abnormal_nova-compute_down.2.log 1
./go.sh abnormal_volume_limit_exceed.1.log 1
./go.sh abnormal_volume_limit_exceed.2.log 1

