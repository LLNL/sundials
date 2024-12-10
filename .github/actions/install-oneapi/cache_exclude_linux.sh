#!/bin/bash

# SPDX-FileCopyrightText: 2020 Intel Corporation
#
# SPDX-License-Identifier: MIT

#shellcheck disable=SC2010
LATEST_VERSION=$(ls -1 /opt/intel/oneapi/compiler/ | grep -v latest | sort | tail -1)

sudo rm -rf /opt/intel/oneapi/compiler/"$LATEST_VERSION"/linux/compiler/lib/ia32_lin
sudo rm -rf /opt/intel/oneapi/compiler/"$LATEST_VERSION"/linux/bin/ia32
sudo rm -rf /opt/intel/oneapi/compiler/"$LATEST_VERSION"/linux/lib/emu
sudo rm -rf /opt/intel/oneapi/compiler/"$LATEST_VERSION"/linux/lib/oclfpga
