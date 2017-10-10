# Makefile for ebacoli exclusive examples

# Fortran compiler to be used
FC=gfortran
# Compile flags, code must be aware of ebacoli library and module location
FCFLAGS=-I../../../src
# Linker flags
FLFLAGS=-Wl,-R,'$$ORIGIN/../../../src'

# Recompile examples if libbacoli.so has changed.
# Ensure module gets linked as well.
BACOLI=../../../src/libebacoli.so

EXECUTABLES:=   TT_monodomain-simple	\
		TT_monodomain-trimesh

SRC_FILES:=TenTusscher.f TT_ebacoli.f f_d_COL10.f f_d_COL12.f f_d_COL14.f f_d_COL16.f f_d_COL18.f f_d_COL2.f f_d_COL4.f f_d_COL6.f f_d_COL8.f f_d_COL1.f f_d_COL11.f f_d_COL13.f f_d_COL15.f f_d_COL17.f f_d_COL19.f f_d_COL3.f f_d_COL5.f f_d_COL7.f f_d_COL9.f
OBJ_FILES:=${SRC_FILES:.f=.o}


.PHONY: all debug opt clean-obj clean test

all: ${EXECUTABLES}

debug: FCFLAGS+=-g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
debug: all

opt: FCFLAGS+=-O2
opt: all

TT_monodomain-simple: driver-TT_monodomain-simple.o ${OBJ_FILES} ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

TT_monodomain-trimesh: driver-TT_monodomain-trimesh.o ${OBJ_FILES} ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

%.o: %.f
	${FC} ${FCFLAGS} -o $@ -c $<

%.o: %.f95
	${FC} ${FCFLAGS} -o $@ -c $<

clean-obj:
	rm -rf *.o

clean: clean-obj
	rm -rf ${EXECUTABLES}
