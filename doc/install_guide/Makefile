# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?= -W --keep-going
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make all".
all:
	@$(SPHINXBUILD) -M latexpdf "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

.PHONY: all Makefile

latexpdf: all

clean: Makefile
	@$(SPHINXBUILD) -M clean "$(SOURCEDIR)" "$(BUILDDIR)"