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
#:  BUILD           Build variant 'mpi', 'debug', 'std'... (default: %BUILD%)
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
#:To build a sequential version, you must explicitly set BUILD=std.
#:
#:You may add options to the 'waf' commands by using the OPTS environment variable
#:on the command line (example: with a progress bar):
#:      make safe OPTS='-p'
#:or by setting them before (example: sequential build, limit to 4 tasks):
#:      export BUILD=std OPTS='-j 4'
#:      make
#:
#:Execute a testcase:
#:      make test n=ssll112a
#:or:
#:      make ssll112a

BUILD ?= mpi
OPTS ?=
DEFAULT ?= safe

SHELL = /bin/bash

.PHONY: help default bootstrap bootstrap_debug all
# targets for BUILD configuration
.PHONY: configure install safe fast doc distclean install-tests test
# same targets for 'debug' configuration
.PHONY: configure_debug install_debug safe_debug fast_debug doc_debug distclean_debug install-tests_debug test_debug

default: $(DEFAULT)

all:
	@make BUILD=mpi bootstrap
	@make BUILD=debug bootstrap

bootstrap: configure safe doc

configure:
	./waf_$(BUILD) configure $(OPTS)

install: safe

safe:
	./waf_$(BUILD) install $(OPTS) --safe -j $$(nproc)

fast:
	./waf_$(BUILD) install $(OPTS) --fast -j $$(nproc)

doc:
	@( \
		if [ $(BUILD) = "std" ]; then \
			echo "doc skipped, only available in parallel" ; \
			exit 0 ; \
		fi ; \
		./waf_$(BUILD) doc $(OPTS) ; \
	)

distclean: ##- perform a distclean of the build directory.
	./waf_$(BUILD) distclean

install-tests:
	@make fast OPTS="$(OPTS) --install-tests"

n ?=
test:
	@( \
		if [ -z "$(n)" ]; then \
			echo "usage: make test n=testname" ; \
			exit 1 ; \
		fi ; \
		./waf_$(BUILD) test $(OPTS) -n $(n) ; \
	)

bootstrap_debug: configure_debug safe_debug doc_debug

configure_debug:
	@make BUILD=debug configure

install_debug:
	@make BUILD=debug install

safe_debug:
	@make BUILD=debug safe

fast_debug:
	@make BUILD=debug fast

doc_debug:
	@make BUILD=debug doc

distclean_debug:
	@make BUILD=debug distclean

install-tests_debug:
	@make BUILD=debug install-tests"

test_debug:
	@make BUILD=debug test n=$(n)

help : makefile
	@sed -n 's/^#://p' $< | \
		sed -e 's/%BUILD%/$(BUILD)/g' -e 's/%DEFAULT%/$(DEFAULT)/g'

%:
	@make --no-print-directory test n="$@"
