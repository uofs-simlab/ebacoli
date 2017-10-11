# Makefile for Competition-diffusion model

# Fortran compiler to be used
FC=gfortran
# Compile flags, code must be aware of bacoli library and module location
FCFLAGS=-I../../lib
# Linker flags
FLFLAGS=-Wl,-R,'$$ORIGIN/../../lib'

# Recompile examples if libbacoli.so has changed.
# Ensure module gets linked as well.
BACOLI=../../lib/libebacoli.so

EXECUTABLES:= CDM-frameoutput

.PHONY: all run-CDM-frameoutput heatmap movie frames debug opt clean-obj clean-frames clean

all: run-CDM-frameoutput heatmap frames movie

debug: FCFLAGS+=-g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
debug: all

opt: FCFLAGS+=-O2
opt: all

CDM-frameoutput: driver-CDM-frameoutput.o CDM.o ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

run-CDM-frameoutput: CDM-frameoutput
	./CDM-frameoutput
%.o: %.f
	${FC} ${FCFLAGS} -o $@ -c $<

%.o: %.f95
	${FC} ${FCFLAGS} -o $@ -c $<

clean-obj:
	rm -rf *.o

heatmap: plot_heatmap_from_Bsplines.py
	./plot_heatmap_from_Bsplines.py v 0 Bsplines000???

frames: plot_all_bsplines.py
	./plot_all_bsplines.py Bsplines000???

movie:
	@avconv -i 'Bsplines%06d.png' -r 24 -c:v libx264 Bsplines_movie.mp4

clean-frames:
	rm -rf *.png

clean: clean-obj clean-frames
	rm -rf Bsplines??????
	rm -rf Bsplines_movie.mp4
	rm -rf ${EXECUTABLES}
