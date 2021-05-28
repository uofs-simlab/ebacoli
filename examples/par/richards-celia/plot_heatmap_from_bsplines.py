#!/usr/bin/env python

# All we need is the header file in the parent directory
import sys, os
from nice_plots_header import *

############################################################################
##
## End of import/settings
##
############################################################################

# Customization for this figure
cmap = colors.cmap_widesep_r

############################################################################
##
## End of setup
##
############################################################################

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: plot_heatmap_from_Bsplines.py <list of Bspline filenames>")
        sys.exit(1)

    fnames = sys.argv[1:]
    n_files = len(fnames)

    # read all data files
    ed_list = []
    for ii in range(n_files):
        ed_list.append(EbacoliData(fnames[ii]))

    title = "$\psi(x,t)$"
    plot_heatmap_of_field(ed_list, "u", 0, cmap=cmap, title=title)
