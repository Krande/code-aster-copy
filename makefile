#:This makefile gives convenient shortcuts to build and test code_aster.
#:
#:Targets:
#:  default         Default task, selected by DEFAULT environment variable: 'safe' by default.
#:  bootstrap       Run 'configure', 'safe' and 'doc'.
#:  all             Run 'bootstrap' for 'mpi' and 'debug' configurations.
#:
#:  configure       Configure the project: 'waf configure'
#:  install         Alias for 'safe'
#:  safe            Build and install in safe mode: 'waf install --safe'
#:  fast            Build and install in fast mode: 'waf install --fast'
#:                  The difference between safe and fast modes is the algorithm used
#:                  to check changes of dependencies.
#:  doc             Build the embedded Html documentation: 'waf doc'
#:  distclean       Perform a distclean of the build directory: 'waf distclean'
#:  install-tests   Same as 'fast' by adding the '--install-tests' option
#:  test            Execute a testcase, use 'n' variable: 'waf test -n xxx'
#:
#:The same targets exist with '_debug' suffix (configure_debug, install_debug, test_debug...)
#:which use 'waf_debug' instead of 'waf'.
#:
#:  help            Show this help message
#:
#:  <testname>      An unknown target is treated as a testname, same 'make test n=testname'
#:
#:Environment variables:
#:  ASTER_CONFIG    Environment file to be sourced for building (default: automatic detection)
#:  ASTER_BUILD     Build variant 'mpi', 'debug' (default: %ASTER_BUILD%)
#:  DEFAULT         Default selected target (default: %DEFAULT%)
#:  OPTS            Options passed to waf commands, example OPTS='-p'
#:
#:With all prerequisites well configured (example in a up-to-date container) you may run:
#:      ./configure
#:      make
#:or
#:      make bootstrap
#:
#:Build both optimized and debug versions:
#:      make all
#:
#:To build a sequential version, you must explicitly set ASTER_BUILD=std (but you
#:can not build std+debug using this makefile).
#:
#:You may add options to the 'waf' commands by using the OPTS environment variable
#:on the command line (example: with a progress bar):
#:      make safe OPTS='-p'
#:The number of jobs is directly passed to make (example: sequential build, limit to 4 tasks):
#:      export ASTER_BUILD=std
#:      make -j 16
#:
#:Execute a testcase:
#:      make test n=ssll112a
#:or:
#:      make ssll112a

# for compatibility: use value of BUILD if ASTER_BUILD is not defined
ifeq ($(origin ASTER_BUILD), undefined)
	BUILD ?= mpi
	ifneq (,$(filter $(BUILD),mpi debug std))
		ASTER_BUILD := $(BUILD)
	else
		ASTER_BUILD := mpi
	endif
endif

OPTS ?=
# extract '-j' option to be passed to waf
JOBS ?= $(shell \
	j="-j"; \
	if grep -q -- "-j" <<< "$(MAKEFLAGS)"; then \
		j=-j$$( sed -e 's/.*-j\([0-9]\+\).*/\1/' <<< "$(MAKEFLAGS)" ) ; \
	fi; \
	if [ "$$j" = "-j" ]; then \
		jobs=$$(nproc) ; \
		[ -z "$$jobs" ] && jobs=12 ; \
		[ $$jobs -gt 12 ] && jobs=12 ; \
		j="-j$$jobs" ; \
	fi ; \
	echo $$j )
DEFAULT ?= safe

SHELL = /bin/bash

.PHONY: help default bootstrap bootstrap_debug all
# targets for ASTER_BUILD configuration
.PHONY: configure install safe fast build doc distclean install-tests test
# same targets for 'debug' configuration
.PHONY: configure_debug install_debug safe_debug fast_debug build_debug doc_debug distclean_debug
.PHONY: install-tests_debug test_debug
.PHONY: clean_cache

default: $(DEFAULT)

all:
	$(MAKE) ASTER_BUILD=mpi bootstrap
	$(MAKE) ASTER_BUILD=debug bootstrap

bootstrap: configure safe doc

configure:
	./waf_$(ASTER_BUILD) configure $(OPTS)

install: safe

safe:
	./waf_$(ASTER_BUILD) install $(OPTS) --safe $(JOBS)

fast:
	./waf_$(ASTER_BUILD) install $(OPTS) --fast $(JOBS)

build:
	./waf_$(ASTER_BUILD) build $(OPTS) --safe $(JOBS)

doc:
	@( \
		if [ $(ASTER_BUILD) = "std" ]; then \
			echo "doc skipped, only available in parallel" ; \
			exit 0 ; \
		fi ; \
		./waf_$(ASTER_BUILD) doc $(OPTS) ; \
	)

distclean: ##- perform a distclean of the build directory.
	./waf_$(ASTER_BUILD) distclean

install-tests:
	$(MAKE) fast OPTS="$(OPTS) --install-tests"

n ?=
test:
	@( \
		if [ -z "$(n)" ]; then \
			echo "usage: make test n=testname" ; \
			exit 1 ; \
		fi ; \
		./waf_$(ASTER_BUILD) test $(OPTS) -n $(n) ; \
	)

bootstrap_debug: configure_debug safe_debug doc_debug

configure_debug:
	$(MAKE) ASTER_BUILD=debug configure

install_debug:
	$(MAKE) ASTER_BUILD=debug install

safe_debug:
	$(MAKE) ASTER_BUILD=debug safe

fast_debug:
	$(MAKE) ASTER_BUILD=debug fast

build_debug:
	$(MAKE) ASTER_BUILD=debug build

doc_debug:
	$(MAKE) ASTER_BUILD=debug doc

distclean_debug:
	$(MAKE) ASTER_BUILD=debug distclean

clean_cache:
	@rm -rf $$(find * -type d -name __pycache__ ! -path 'build/*')

install-tests_debug:
	$(MAKE) ASTER_BUILD=debug install-tests

test_debug:
	$(MAKE) ASTER_BUILD=debug test n=$(n)

help: makefile
	@sed -n 's/^#://p' $< | \
		sed -e 's/%ASTER_BUILD%/$(ASTER_BUILD)/g' -e 's/%DEFAULT%/$(DEFAULT)/g'

%:
	@( \
		if [ ! -f ../src/astest/"$@".export ] && [ ! -f ../validation/astest/"$@".export ]; then \
			echo "not a testcase: $@"; \
			exit 1; \
		fi; \
	)
	$(MAKE) --no-print-directory test n="$@"

# parallelism is enabled by waf
.NOTPARALLEL:
