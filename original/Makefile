FC = gfortran
FFLAGS ?= -O2 -std=f2008 -ffree-line-length-none -Wall -Wextra -Wno-compare-reals

.PHONY: all debug test study clean
all: electrospinning

electrospinning: electrospinning.f90
	$(FC) $(FFLAGS) $< -o $@

debug: electrospinning.f90
	$(FC) $(FFLAGS) -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace $< -o electrospinning_debug

test: debug
	python3 tests/verify.py --exe ./electrospinning_debug
	python3 tests/verify_plasma.py --exe ./electrospinning_debug

study: electrospinning
	python3 mesh_study.py --exe ./electrospinning --config config_spinning.nml

clean:
	rm -f electrospinning electrospinning_debug *.o *.mod
