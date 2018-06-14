# Makefile for ebacoli exclusive examples

# Fortran compiler to be used
FC=gfortran
# Compile flags, code must be aware of bacoli library and module location
FCFLAGS=-I../../../lib
# Linker flags
FLFLAGS=-Wl,-R,'$$ORIGIN/../../../lib'

# Recompile examples if libbacoli.so has changed.
# Ensure module gets linked as well.
BACOLI=../../../lib/libebacoli.so

EXECUTABLES:= TTP_monodomain-frameoutput

SRC_FILES:=TenTusscher.f TTP_ebacoli.f f_d_COL10.f f_d_COL12.f f_d_COL14.f f_d_COL16.f f_d_COL18.f f_d_COL2.f f_d_COL4.f f_d_COL6.f f_d_COL8.f f_d_COL1.f f_d_COL11.f f_d_COL13.f f_d_COL15.f f_d_COL17.f f_d_COL19.f f_d_COL3.f f_d_COL5.f f_d_COL7.f f_d_COL9.f
OBJ_FILES:=$(addprefix TTP_src/,${SRC_FILES:.f=.o})


.PHONY: all TTP_monodomain-frameoutput heatmap frame debug opt clean-obj clean test

all: run-TTP_monodomain-frameoutput heatmap frame

debug: FCFLAGS+=-g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
debug: all

opt: FCFLAGS+=-O2
opt: all

TTP_monodomain-frameoutput: driver-TTP_monodomain-frameoutput.o ${OBJ_FILES} ${BACOLI}
	${FC} ${FLFLAGS} -o $@ $^

run-TTP_monodomain-frameoutput: TTP_monodomain-frameoutput
	./TTP_monodomain-frameoutput
%.o: %.f
	${FC} ${FCFLAGS} -o $@ -c $<

%.o: %.f95
	${FC} ${FCFLAGS} -o $@ -c $<

clean-obj:
	rm -rf ${OBJ_FILES}

heatmap: plot_heatmap_from_Bsplines.py
	./plot_heatmap_from_Bsplines.py u 0 Bsplines000???; done
	for ii in `seq 0 17`; do ./plot_heatmap_from_Bsplines.py v $${ii} Bsplines000???; done;

frame: plot_all_bsplines.py
	./plot_all_bsplines.py png Bsplines000???

movie:
	@avconv -i 'Bsplines%06d.png' -r 24 -c:v libx264 Bsplines_movie.mp4

clean-frames:
	rm -rf *.png

clean: clean-obj clean-frames
	rm -rf Bsplines??????
	rm -rf ${EXECUTABLES}
