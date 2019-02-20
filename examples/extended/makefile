# Makefile for eBACOLI exclusive examples

# Fortran compiler to be used
FC=gfortran
# Compile flags, code must be aware of ebacoli library and module location
FCFLAGS=-I../../lib
# Linker flags
FLFLAGS=-Wl,-R,'$$ORIGIN/../../lib'

# Recompile examples if libbacoli.so has changed.
# Ensure module gets linked as well.
BACOLI=../../lib/libebacoli.so

EXECUTABLES:=   FHN_monodomain-trimesh  	\
		FHN2_monodomain-trimesh  	\
		toy-trimesh			\
		toy2-trimesh			\
		CDM-trimesh		\
		rabid-trimesh			\
		FHN_monodomain-simple		\
		FHN2_monodomain-simple		\
		toy-simple			\
		toy2-simple			\
		toy3-simple			\
		CDM-simple			\
		test-app

.PHONY: all debug opt clean-obj clean test

all: ${EXECUTABLES}

debug: FCFLAGS+=-g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
debug: all

opt: FCFLAGS+=-O2
opt: all

FHN_monodomain-simple: driver-FHN_monodomain-simple.o FHN_monodomain.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

FHN2_monodomain-simple: driver-FHN2_monodomain-simple.o FHN2_monodomain.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

toy-simple: driver-toy-simple.o toy.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

toy2-simple: driver-toy2-simple.o toy2.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

toy3-simple: driver-toy3-simple.o toy2.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

CDM-simple: driver-CDM-simple.o CDM.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^


FHN_monodomain-trimesh: driver-FHN_monodomain-trimesh.o FHN_monodomain.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

FHN2_monodomain-trimesh: driver-FHN2_monodomain-trimesh.o FHN2_monodomain.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

toy-trimesh: driver-toy-trimesh.o toy.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

toy2-trimesh: driver-toy2-trimesh.o toy2.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

CDM-trimesh: driver-CDM-trimesh.o CDM.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

rabid-trimesh: driver-rabid-trimesh.o rabid.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

test-app: test-app-driver.o test-app.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^


%.o: %.f
	${FC} ${FCFLAGS} -o $@ -c $<

%.o: %.f95
	${FC} ${FCFLAGS} -o $@ -c $<

clean-obj:
	rm -rf *.o

clean: clean-obj
	rm -rf ${EXECUTABLES}
