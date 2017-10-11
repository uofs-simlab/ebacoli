# Makefile for ebacoli library

LIBNAME=libebacoli.so
MODNAME=ebacoli95_mod.mod

FC=gfortran
FCFLAGS=-fPIC
FLFLAGS=-Wl,-soname,${LIBNAME}

EBACOLI_F_SRC=ebacoli.f ebacoli-aux.f
EBACOLI_F95_SRC=ebacoli95.f95 d1mach_i1mach.f95

EBACOLI_F_OBJ=${EBACOLI_F_SRC:.f=.o}
EBACOLI_F95_OBJ=${EBACOLI_F95_SRC:.f95=.o}

.PHONY: all debug opt clean-obj clean

all: ${LIBNAME}

debug: FCFLAGS+=-g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
debug: ${LIBNAME}

opt: FCFLAGS+=-O2
opt: ${LIBNAME}

%.o: %.f
	${FC} ${FCFLAGS} -o $@ -c $<

%.o: %.f95
	${FC} ${FCFLAGS} -o $@ -c $<

${LIBNAME}: ${EBACOLI_F_OBJ} ${EBACOLI_F95_OBJ}
	${FC} -shared ${FLFLAGS} -o $@ $^
	cp ${LIBNAME} ../lib/${LIBNAME}
	cp ${MODNAME} ../lib/${MODNAME}

clean-obj:
	rm -rf *.o

clean: clean-obj
	rm -rf ../lib/${LIBNAME} ../lib/${MODNAME}
	rm -rf ${LIBNAME} ${MODNAME}
