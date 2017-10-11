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
fig_bot = -0.2
fig_top = 1.2
leg_loc = "upper right"
nint_text_loc = [78,0.1]

u_to_plot = range(2)
v_to_plot = range(1)
vline_positions = [50]

############################################################################
##
## End of setup
##
############################################################################

if __name__=="__main__":

    if len(sys.argv) < 2:
        print "Usage: plot_all_bsplines.py <list of Bspline filenames>"
        sys.exit(1)

    fnames = sys.argv[1:]

    for ii in range(len(fnames)):

        ed = EbacoliData(fnames[ii])

        # mark the vlines you want
        for v in vline_positions:
            plt.axvline(x=v,color='k',linestyle='--',linewidth=0.5)

        plot_interpolated_fields(ed, u_to_plot, v_to_plot,
                                 u_labels=["$u$","$v$"],
                                 u_colors = [colors.presen_green,colors.yellow],
                                 v_labels=["$w$"], v_colors = ["k"],
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
        plt.savefig(fnames[ii]+'.png')
        plt.gcf().clear()
