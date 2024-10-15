LOGS = /tmp/coupling.log
ENV = env_ple.sh

default: $(LOGS)

initialize: .initialized

$(LOGS): run

run: .initialized clean
	mpiexec --tag-output -n 1 run_aster coupled_mechanics.py : \
	                     -n 1 run_aster coupled_thermics.py | tee $@

.initialized:
	export CA_COUPLING_LOGDIR=$(PWD)

mesh/solid_new_0.med: mesh/solid.med
	mpiexec --tag-output -n 4 run_aster mesh/split.py

split: mesh/solid_new_0.med

clean:
	rm -f *.log *.log.*

distclean: clean
	rm -f mesh/*_new_*.med

.PHONY: default clean distclean split $(LOGS)
