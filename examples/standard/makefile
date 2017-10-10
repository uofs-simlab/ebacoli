# Makefile for eBACOLI examples that could have been solved with BACOLI

# Fortran compiler to be used
FC=gfortran
# Compile flags, code must be aware of ebacoli library and module location
FCFLAGS=-I../../lib
# Linker flags
FLFLAGS=-Wl,-R,'$$ORIGIN/../../lib'

# Recompile examples if libbacoli.so has changed.
# Ensure module gets linked as well.
BACOLI=../../src/libebacoli.so

EXECUTABLES:= burg1-simple         \
		burg2-simple       \
		steady-simple      \
		cahn_allen-simple  \
		rcdsys-simple      \
		burg1-curve        \
		burg2-curve        \
		steady-curve       \
		cahn_allen-curve   \
		rcdsys-curve       \
		burg1-trimesh      \
		burg2-trimesh      \
		steady-trimesh     \
		cahn_allen-trimesh \
		rcdsys-trimesh

.PHONY: all debug opt clean-obj clean test

all: ${EXECUTABLES}

debug: FCFLAGS+=-g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
debug: all

opt: FCFLAGS+=-O2
opt: all

# From the simple driver
burg1-simple: driver95-simple.o burg1.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

burg2-simple: driver95-simple.o burg2.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

steady-simple: driver95-simple.o steady.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

cahn_allen-simple: driver95-simple.o cahn_allen.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

rcdsys-simple: driver95-simple.o rcdsys.o ${BACOLI}
	rm -rf driver95-simple.o
	${FC} ${FCFLAGS} -DMAC_NPDE=4 -o driver95-simple.o -c driver95-simple.F95
	${FC} ${FLFLAGS} -DMAC_NPDE=4 -o $@ $^
	rm -rf driver95-simple.o

# From the curve driver
burg1-curve: driver95-curve.o burg1.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

burg2-curve: driver95-curve.o burg2.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

steady-curve: driver95-curve.o steady.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

cahn_allen-curve: driver95-curve.o cahn_allen.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

rcdsys-curve: driver95-curve.o rcdsys.o ${BACOLI}
	rm -rf driver95-curve.o
	${FC} ${FCFLAGS} -DMAC_NPDE=4 -o driver95-curve.o -c driver95-curve.F95
	${FC} ${FLFLAGS} -DMAC_NPDE=4 -o $@ $^
	rm -rf driver95-curve.o

# From the trimesh driver
burg1-trimesh: driver95-trimesh.o burg1.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

burg2-trimesh: driver95-trimesh.o burg2.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

steady-trimesh: driver95-trimesh.o steady.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

cahn_allen-trimesh: driver95-trimesh.o cahn_allen.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

rcdsys-trimesh: driver95-trimesh.o rcdsys.o ${BACOLI}
	rm -rf driver95-trimesh.o
	${FC} ${FCFLAGS} -DMAC_NPDE=4 -o driver95-trimesh.o -c driver95-trimesh.F95
	${FC} ${FLFLAGS} -DMAC_NPDE=4 -o $@ $^
	rm -rf driver95-trimesh.o


%.o: %.F
	${FC} ${FCFLAGS} -o $@ -c $<

%.o: %.F95
	${FC} ${FCFLAGS} -o $@ -c $<

clean-obj:
	rm -rf *.o

clean: clean-obj
	rm -rf ${EXECUTABLES}

########################################################################
# TESTING
#
test: ${EXECUTABLES}
	@echo "########################################################################"
	@echo "## STARTING TEST RUNS ##################################################"
	@echo "########################################################################"
	for e in ${EXECUTABLES}; do echo $${e}; ./runtest.sh $${e}; done;
