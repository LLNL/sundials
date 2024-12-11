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

fn="../doc/shared/sundials_vars.py"
sedi "s/doc_version =.*/doc_version = \'develop\'/" $fn

# ------------------------------------------------------------------------------
# Update Markdown changelog
# ------------------------------------------------------------------------------

fn="../CHANGELOG.md"

# Add new entry to changelog
cat > tmp.txt <<HEREDOC
# SUNDIALS Changelog

## Changes to SUNDIALS in release X.Y.Z

### Major Features

### New Features and Enhancements

### Bug Fixes

### Deprecation Notices
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

# Replace line containing RecentChanges_link.rst in Changelog.rst with the
# contents of RecentChanges.rst
sedi -e '/RecentChanges_link.rst/ {' \
     -e 'r ../doc/shared/RecentChanges.rst' \
     -e 'd' \
     -e '}' \
     $fn

# Clear recent changes file
cat > ../doc/shared/RecentChanges.rst <<HEREDOC
**Major Features**

**New Features and Enhancements**

**Bug Fixes**

**Deprecation Notices**
HEREDOC

# Create temporary file with new changelog entry
cat > tmp.txt <<HEREDOC
.. SED_REPLACEMENT_KEY

Changes to SUNDIALS in release X.Y.Z
====================================

.. include:: RecentChanges_link.rst
HEREDOC

# Replace the line containing SED_REPLACEMENT with the new changelog entry
sedi -e '/SED_REPLACEMENT_KEY/ {' \
     -e 'r tmp.txt' \
     -e 'd' \
     -e '}' \
     $fn

rm -f tmp.txt

# ------------------------------------------------------------------------------
# Update package introductions
# ------------------------------------------------------------------------------

# Replace section titles
for pkg in arkode cvode cvodes ida idas kinsol
do
    sedi 's/Changes to SUNDIALS.*/Changes to SUNDIALS in release X.Y.Z/I' \
         "../doc/${pkg}/guide/source/Introduction.rst"
done
