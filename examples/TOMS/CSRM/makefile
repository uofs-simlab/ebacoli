# Makefile for Competition-diffusion model

# Fortran compiler to be used
FC=gfortran
# Compile flags, code must be aware of bacoli library and module location
FCFLAGS=-I../../../lib
# Linker flags
FLFLAGS=-Wl,-R,'$$ORIGIN/../../../lib'

# Recompile examples if libbacoli.so has changed.
# Ensure module gets linked as well.
BACOLI=../../../lib/libebacoli.so

EXECUTABLES:= CSRM-frameoutput

.PHONY: all CSRM-frameoutput heatmaps frame debug opt clean-obj clean test

all: run-CSRM-frameoutput heatmap frame

debug: FCFLAGS+=-g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
debug: all

opt: FCFLAGS+=-O2
opt: all

CSRM-frameoutput: driver-CSRM-frameoutput.o CSRM.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

run-CSRM-frameoutput: CSRM-frameoutput
	./CSRM-frameoutput
%.o: %.f
	${FC} ${FCFLAGS} -o $@ -c $<

%.o: %.f95
	${FC} ${FCFLAGS} -o $@ -c $<

clean-obj:
	rm -rf *.o

heatmap: plot_heatmap_from_Bsplines.py
	for ii in `seq 0 3`; do ./plot_heatmap_from_Bsplines.py u $${ii} Bsplines000???; done;

frame: plot_all_bsplines.py
	./plot_all_bsplines.py png Bsplines000???

movie:
	@avconv -i 'Bsplines%06d.png' -r 24 -c:v libx264 Bsplines_movie.mp4


clean-frames:
	rm -rf *.png

clean: clean-obj clean-frames
	rm -rf Bsplines??????
	rm -rf ${EXECUTABLES}
