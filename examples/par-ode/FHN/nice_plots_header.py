"""
Header to handle consistent figure styling for our ebacoli plots.
"""

# *NOTE* colors sets rc things, load it first so that you can overwrite
# afterwards if needed
# the python libs of ebacoli must be accessible within python (i.e., by PYTHONPATH or other)
import colors
from ebacoli import *

# MATPLOTLIBRC stuff
import matplotlib as mpl
# mpl.use('AGG')  # for systems not running a GUI

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}
mpl.rc('font', **font)

# PYPLOTRC stuff
import matplotlib.pyplot as plt
plt.rc('lines', linewidth=2)
# Latex formatting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

from matplotlib import cm
