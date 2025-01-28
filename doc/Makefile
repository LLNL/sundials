# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# Makefile to build or clean all RST package docs at once
# ------------------------------------------------------------------------------

DIRS = arkode/guide arkode/examples \
       cvode/guide \
       cvodes/guide \
       ida/guide \
       idas/guide \
       kinsol/guide

# prefix with desired command to create unique dependencies/targets
LATEXDIRS = $(DIRS:%=latexpdf-%) latexpdf-install_guide
HTMLDIRS = $(DIRS:%=html-%) html-superbuild
CLEANDIRS = $(DIRS:%=clean-%) clean-superbuild clean-install_guide

latexpdf: $(LATEXDIRS)
$(LATEXDIRS):
	$(MAKE) -C $(@:latexpdf-%=%) latexpdf

html: $(HTMLDIRS)
$(HTMLDIRS):
	$(MAKE) -C $(@:html-%=%) html

clean: $(CLEANDIRS)
$(CLEANDIRS):
	$(MAKE) -C $(@:clean-%=%) clean

.PHONY: $(LATEXDIRS)
.PHONY: $(HTMLDIRS)
.PHONY: $(CLEANDIRS)
.PHONY: latexpdf html clean
