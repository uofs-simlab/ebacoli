#!/usr/bin/env python

import numpy as np
import sys, os

# lower, upper, npts as input args
x_lower = float(sys.argv[1])
x_upper = float(sys.argv[2])
num_pts = int(sys.argv[3])

x = np.linspace(x_lower,x_upper,num_pts)

for i in x:
    print(i)
