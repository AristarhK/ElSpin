FC = ./scripts/fortran.sh
PYTHON ?= python3
FFLAGS ?= -O2 -std=f2008 -ffree-line-length-none -Wall -Wextra
DEBUG_FLAGS = -O1 -g -std=f2008 -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace

.PHONY: all debug test test-original demo clean
all: build/jet3d build/electrospinning
build:
	mkdir -p build
build/jet3d: src/jet3d.f90 | build
	$(FC) $(FFLAGS) -Jbuild $< -o $@
build/electrospinning: original/electrospinning.f90 | build
	$(FC) $(FFLAGS) -Jbuild $< -o $@
build/jet3d_debug: src/jet3d.f90 | build
	$(FC) $(DEBUG_FLAGS) -Jbuild $< -o $@
build/electrospinning_debug: original/electrospinning.f90 | build
	$(FC) $(DEBUG_FLAGS) -Jbuild $< -o $@
debug: build/jet3d_debug build/electrospinning_debug
test: all build/jet3d_debug
	build/jet3d_debug --self-test
	MPLCONFIGDIR=/tmp/electrospin-mpl $(PYTHON) -m unittest discover -s tests -v
test-original: build/electrospinning_debug
	$(PYTHON) original/tests/verify.py --exe build/electrospinning_debug
	$(PYTHON) original/tests/verify_plasma.py --exe build/electrospinning_debug
demo: all
	MPLCONFIGDIR=/tmp/electrospin-mpl $(PYTHON) -m electrospin run --config configs/demo.json --out results/demo
clean:
	rm -rf build
