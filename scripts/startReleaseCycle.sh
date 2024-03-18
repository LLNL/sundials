#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# Update files to start new release cycle
# ------------------------------------------------------------------------------

set -e
set -o pipefail

# Wrapper for editing inplace with different sed implementations
sedi() {
    case $(uname) in
        Darwin*) sedi=('-i' '') ;;
        *) sedi='-i' ;;
    esac
    sed "${sedi[@]}" "$@"
}

# ------------------------------------------------------------------------------
# Update versions
# ------------------------------------------------------------------------------

fn="../doc/shared/versions.py"
sedi "s/doc_version =.*/doc_version = \'develop\'/" $fn

# ------------------------------------------------------------------------------
# Update Markdown changelog
# ------------------------------------------------------------------------------

fn="../CHANGELOG.md"

# Add new entry to changelog
cat > tmp.txt <<HEREDOC
# SUNDIALS Changelog

## Changes to SUNDIALS in release X.Y.Z

### New Features

### Bug Fixes
HEREDOC

sedi -e '/SUNDIALS Changelog/ {' \
     -e 'r tmp.txt' \
     -e 'd' \
     -e '}' \
     $fn

rm -f tmp.txt

# ------------------------------------------------------------------------------
# Update RST changelog
# ------------------------------------------------------------------------------

fn="../doc/shared/Changelog.rst"

# Move recent changes to changelog
sedi -e '/RecentChanges_link.rst/ {' \
     -e 'r ../doc/shared/RecentChanges.rst' \
     -e 'd' \
     -e '}' \
     $fn

# Clear recent changes file
cat > ../doc/shared/RecentChanges.rst <<HEREDOC
**New Features**

**Bug Fixes**
HEREDOC

# Add new entry to changelog
cat > tmp.txt <<HEREDOC
.. SED_REPLACEMENT_KEY

Changes to SUNDIALS in release X.Y.Z
====================================

.. include:: RecentChanges_link.rst
HEREDOC

sedi -e '/SED_REPLACEMENT_KEY/ {' \
     -e 'r tmp.txt' \
     -e 'd' \
     -e '}' \
     $fn

rm -f tmp.txt
