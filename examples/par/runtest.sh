#!/usr/bin/env sh

# run test for a given executable

EXECUTABLE=$1
OUTPUTDIR=output/exact_jac

./${EXECUTABLE} > ${EXECUTABLE}.out
# test standard output (from each type of driver)
diff ${OUTPUTDIR}/${EXECUTABLE}.out ${EXECUTABLE}.out || printf "${PWD}\nPossible problem with ${EXECUTABLE}, diffs above\n=========================================\n"; rm ${EXECUTABLE}.out

# test Point outputs (from trimesh driver)
if [ -e Points1 ]
then
   for p in `ls Points?`
   do
       diff ${OUTPUTDIR}/${EXECUTABLE}_${p}.out ${p} || printf "${PWD}\nPossible problem with ${EXECUTABLE}, diffs above\n=========================================\n"; rm ${p}
   done
fi

# test Bspline output (from curve driver)
if [ -e Bsplines ]
then
   diff ${OUTPUTDIR}/${EXECUTABLE}_Bsplines.out Bsplines || printf "${PWD}\nPossible problem with ${EXECUTABLE}, diffs above\n=========================================\n"; rm Bsplines
fi
