#!/usr/bin/env python

import sys, os
import numpy as np
import csv
import glob

# EBacoli
from ebacoli import EbacoliData

space_points_file_name = sys.argv[1]

bspline_file_names = sys.argv[2:]

with open(space_points_file_name) as space_points_file:
    reader = csv.reader(space_points_file, delimiter=" ")
    space_points = np.array(list(zip(*reader)),dtype=float)[0]

# Time and space sizes
Nt = len(bspline_file_names)
Nx = len(space_points)

sol_v = np.zeros((Nt*Nx,1))
sol_s = np.zeros((Nt*Nx,1))

# Get all solution points packed into the sol array
i=0
for bspline_file in bspline_file_names:
    ed = EbacoliData(bspline_file)
    # print(len(ref_space_points))
    sol_v[i*Nx:(i+1)*Nx] = ed.u_at_points(space_points,0).reshape(Nx,1)
    sol_s[i*Nx:(i+1)*Nx] = ed.v_at_points(space_points,0).reshape(Nx,1)
    i = i+1


# write solution to stdout
for v,s in zip(sol_v,sol_s):
    print(v[0],s[0])
