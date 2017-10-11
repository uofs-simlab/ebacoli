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
        print "Usage: plot_heatmap_from_Bsplines.py <u or v> <component> <list of Bspline filenames>"
        sys.exit(1)

    plot_var_str = sys.argv[1]   # u or v
    plot_var_num = int(sys.argv[2])
    fnames = sys.argv[3:]
    n_files = len(fnames)

    # read all data files
    ed_list = []
    for ii in range(n_files):
        ed_list.append(EbacoliData(fnames[ii]))

    # Edit title based on component
    if ("%s%d" % (plot_var_str, plot_var_num+1) == "u1"):
        title = "$u(x,t)$"
    elif ("%s%d" % (plot_var_str, plot_var_num+1) == "u2"):
        title = "$v(x,t)$"
    elif ("%s%d" % (plot_var_str, plot_var_num+1) == "v1"):
        title = "$w(x,t)$"


    plot_heatmap_of_field(ed_list, plot_var_str, plot_var_num, cmap=cmap,
                          title=title)
