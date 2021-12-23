#:This makefile gives convenient shortcuts to build and check code_aster.
#:
#:Targets:
#:  configure       Configure the project: 'waf configure'
#:  default         Default task, depending on DEFAULT environment variable: 'safe' or 'fast'
#:  install         Alias for 'default'
#:  bootstrap       Run 'configure', 'safe' and 'doc'.
#:  all             Run 'bootstrap' for 'std' and 'mpi' variants.
#:
#:  safe            Build and install in safe mode: 'waf install --safe'
#:  fast            Build and install in fast mode: 'waf install --fast'
#:  dbg_safe        Build and install in safe mode with debug variant: 'waf install_debug --safe'
#:  dbg             Build and install in fast mode with debug variant: 'waf install_debug --fast'
#:  install-tests   Same as 'fast' by adding the '--install-tests' option
#:
#:  doc             Build the embedded Html documentation: 'waf doc'
#:
#:  test            Execute a testcase, use 'n' variable: 'waf test -n xxx'
#:  check           Execute a list of testcases, use 'LIST' variable.
#:
#:  clean_doc       Remove the build directory of the documentation (debug & release)
#:  distclean       Perform a distclean of the build directory
#:
#:  help            Show this help message
#:
#:  <testname>      An unknown target is treated as a testname, same 'make test n=testname'
#:
#:Environment variables:
#:  BUILD           Build variant 'std' or 'mpi' (default: %BUILD%)
#:  DEFAULT         Algorithm to check changes on dependencies 'safe', 'fast'
#:                  or 'dbg_safe', 'dbg' with debugging informations (default: %DEFAULT%)
#:  OPTS            Options passed to waf commands, example OPTS='-p'
#:
#:With all prerequisites well configured (example in a up-to-date container) you may run:
#:      make
#:or
#:      make bootstrap
#:
#:Build both std & mpi versions:
#:      make all
#:
#:You may add options to the 'waf' commands by using the OPTS environment variable:
#:      make safe BUILD=std OPTS='-p'
#:or:
#:      export BUILD=std OPTS='--progress'
#:      make
#:
#:Execute a testcase:
#:      make test n=ssll112a
#:or:
#:      make ssll112a
#:
#:Check testcases (submit list by default):
#:      make check
#:
#:Check testcases (all 'verification' ones, '-L' is automatically added
#:if LIST does not start with an hyphen):
#:      make check LIST=verification
#:
#:Check testcases (using testlist):
#:      make check LIST="--testlist list.snl"
#:
#:Check testcases (using regexp on testnames):
#:      make check LIST="-R ssnv230[ac]"

BUILD ?= mpi
OPTS ?=
DEFAULT ?= safe

SHELL = /bin/bash

.PHONY: help default configure bootstrap all install-tests test check clean_doc distclean
.PHONY: safe fast dbg_safe dbg doc install

default: $(DEFAULT)
bootstrap: configure safe doc

all:
	@make bootstrap BUILD=std
	@make bootstrap BUILD=mpi

configure:
	./waf_$(BUILD) configure $(OPTS)

install: default

safe:
	./waf_$(BUILD) install $(OPTS) --safe -j $$(nproc)

fast:
	./waf_$(BUILD) install $(OPTS) --fast -j $$(nproc)

dbg_safe:
	./waf_$(BUILD) install_debug $(OPTS) --safe -j $$(nproc)

dbg:
	./waf_$(BUILD) install_debug $(OPTS) --fast -j $$(nproc)

doc:
	./waf_$(BUILD) doc $(OPTS)

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

LIST ?= submit
check:
	@( \
		bindir=$$(python3 -c "ctx = {} ; \
			exec(open('build/$(BUILD)/c4che/release_cache.py').read(), ctx) ; \
			print(ctx['BINDIR'])") ; \
		tmpres=$$(mktemp -d -t run_ctest.XXXXXX) ; rmdir $${tmpres} ; \
		opt=$$(egrep '^-' <<< "$(LIST)" > /dev/null || echo "-L") ; \
		$${bindir}/run_ctest --resutest $${tmpres} $${opt} $(LIST) ; \
		echo "result directory: $${tmpres}" \
	)

clean_doc:
	rm -rf ./build/$(BUILD)/{release,debug}/doc

distclean: ##- perform a distclean of the build directory.
	./waf_$(BUILD) distclean

help : makefile
	@sed -n 's/^#://p' $< | \
		sed -e 's/%BUILD%/$(BUILD)/g' -e 's/%DEFAULT%/$(DEFAULT)/g'

%:
	@make --no-print-directory test n="$@"
