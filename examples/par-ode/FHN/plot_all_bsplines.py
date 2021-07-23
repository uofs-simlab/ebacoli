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
fig_bot = -80
fig_top = 20
leg_loc = "lower left"

nint_text_loc = [0.5,-75]
p_text_loc = [0.85,10]

u_to_plot = range(1)
u_labels  = ["$V$"]
u_colors = [colors.presen_green]

v_to_plot = range(1)
v_labels = ["$S_1$"]
v_colors = ["k"]

w_to_plot = []
w_labels = []
w_colors = []

vline_positions = []

############################################################################
##
## End of setup
##
############################################################################

def print_usage_exit():
    print("Usage: plot_all_bsplines.py <type (png or eps)> <list of Bspline filenames>")
    sys.exit(1)

if __name__=="__main__":

    if len(sys.argv) < 2:
        print_usage_exit()

    extension = sys.argv[1]
    if (extension == "png"):
        pass
    elif (extension == "eps"):
        pass
    else:
        print("Extension must be either 'png' or 'eps'")
        print_usage_exit()

    fnames = sys.argv[2:]

    for ii in range(len(fnames)):

        ed = EbacoliData(fnames[ii])

        # mark the vlines you want
        for v in vline_positions:
            plt.axvline(x=v,color='k',linestyle='--',linewidth=0.5)

        plot_interpolated_fields(ed, u_to_plot, v_to_plot, w_to_plot,
                                 u_labels=u_labels,
                                 u_colors=u_colors,
                                 v_labels=v_labels, v_colors=v_colors,
                                 w_labels=w_labels, w_colors=w_colors,
                                 knot_show=False, knot_bottom=fig_bot)

        # Add ticks with no labels to interior
        plt.xticks(ed.xbs)
        plt.gcf().canvas.draw()
        ax = plt.gca()
        labels = ax.get_xticklabels()
        print(labels[0],labels[-1])
        print(len(labels))
        for i in range(1,len(labels)-1):
            labels[i] = ""
        labels[0] = "0"
        labels[-1] = "1"
        print(labels[0],labels[-1])
        ax.set_xticklabels(labels)

        ########################################################################
        ## Annotating and additional formatting here
        ########################################################################
        plt.text(p_text_loc[0], p_text_loc[1], r'$p=%d$' % ed.p)
        plt.text(nint_text_loc[0], nint_text_loc[1], r'$N=%d$' % ed.nint)

        axes = plt.gca()
        axes.set_ylim([fig_bot,fig_top])
        plt.xlabel(r'$x$')
        plt.ylabel(r'$\{V,S_1\}(x,t)$')
        plt.title(r'$t=%.2f$\qquad\qquad\qquad {\fontsize{%d}{%f}$\texttt{atol}=%1.0e$}' % (ed.time,20,24,ed.atol))
        leg = plt.legend(loc=leg_loc,fancybox=True)
        # set the alpha value of the legend: it will be translucent
        leg.get_frame().set_alpha(0.8)
        plt.tight_layout()
        plt.savefig(fnames[ii]+"."+extension)
        plt.gcf().clear()
