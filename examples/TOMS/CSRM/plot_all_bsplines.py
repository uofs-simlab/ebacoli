#!/usr/bin/env python

# All we need is the header file in the parent directory
import sys, os
from plot_header import *

############################################################################
##
## End of import/settings
##
############################################################################

# Customization for this figure
fig_bot = -0.2
fig_top = 1.2
leg_loc = "upper right"
nint_text_loc = [0.8,0.9]

u_to_plot = range(4)
v_to_plot = []

############################################################################
##
## End of setup
##
############################################################################

if __name__=="__main__":

    if len(sys.argv) < 3:
        print "Usage: plot_all_bsplines.py <extension> <list of Bspline filenames>"
        sys.exit(1)

    extension = sys.argv[1]
    fnames = sys.argv[2:]

    for ii in range(len(fnames)):

        ed = EbacoliData(fnames[ii])

        plot_interpolated_fields(ed, u_to_plot, v_to_plot,
                                 knot_show=True, knot_bottom=fig_bot)

        ########################################################################
        ## Annotating and additional formatting here
        ########################################################################
        plt.text(nint_text_loc[0], nint_text_loc[1], r'$N=%d$' % ed.nint)

        axes = plt.gca()
        axes.set_ylim([fig_bot,fig_top])
        plt.xlabel('$x$')
        plt.title('$t=%.2f$' % ed.time)
        leg = plt.legend(loc="center right",fancybox=True)
        # set the alpha value of the legend: it will be translucent
        leg.get_frame().set_alpha(0.8)
        plt.tight_layout()
        plt.savefig(fnames[ii]+"."+extension)
        plt.gcf().clear()
