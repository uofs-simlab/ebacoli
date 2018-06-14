#!/usr/bin/env python

# All we need is the header file in the parent directory
import sys, os
from plot_header import *

############################################################################
##
## End of import/settings
##
############################################################################


# cmap = uofs.cmap_sep # default color map for all but u1
cmap = cm.bone_r # default color map for all but u1

# separator location - unscaled from original data
sep_bottom  = 11  # chosen to resolve the tip of TT front in variable u1
sep_top     = 17

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

    if ("%s%d" % (plot_var_str, plot_var_num+1) == "u1"):
        title = "$v(x,t)$"
    else:
        title = "$s_{%d}(x,t)$" % (plot_var_num+1)


    ########################################################################
    # cmap to Resolve the tip of TT wavefront (in u1 variable)
    # if (plot_var_str == "u" and plot_var_num == 0) :

    #     # get continuous representation to determine highest and lowest
    #     z = np.zeros([n_files,n_files])
    #     for ii in range(n_files):
    #         good, garb = ed_list[ii].continuous_representation(n_cont=n_files,
    #                                     u_num=[plot_var_num], v_num = [])
    #         z[:,ii] = np.array(good[0])
    #     span = np.max(z) - np.min(z)
    #     low = np.min(z)
    #     high = np.min(z) + span

    #     rel_bot = (sep_bottom-low)/span # relative coordinates for bottom of sep
    #     rel_top = (sep_top-low)/span    # relative coordinates for top of sep

    #     p = [0,0.0,rel_bot,rel_top,1.0,1.0]
    #     cdict = uofs.create_cdict_for_sep_LSC(uofs.green,uofs.black,uofs.yellow,p)
    #     TT_custom_cmap = uofs.LinearSegmentedColormap('ten Tusscher tip resolution', cdict)
    #     plt.register_cmap(cmap=TT_custom_cmap)
    #     cmap = TT_custom_cmap
    ########################################################################

    plot_heatmap_of_field(ed_list, plot_var_str, plot_var_num, cmap=cmap,
                          title=title, extension="eps")
