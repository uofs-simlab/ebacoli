# define colors and some colormaps

from cycler import cycler

# This stuff is used for matplotlib plots
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

############################################################################
# rgb lists
green  = [0.0, 105.0/255, 62.0/255]
yellow = [1.0, 203.0/255, 0.0]
white  = [1.0, 1.0, 1.0]
black  = [0.0, 0.0, 0.0]
grey   = [77.0/255, 78.0/255, 83.0/255]

def create_cdict_for_sep_LSC(bottom, middle, top, sep_points):
    """ Generates dictionary for separated color maps.

    Visually appealing results tend to have the middle color much darker
    than top and bottom

    """
    bc = bottom
    mc = middle
    tc = top
    # 6 important points for transitions
    # The 5 regions go:
    #   0-1 constant
    #   1-2 linear change
    #   2-3 constant (with discontinuity)
    #   3-4 linear change
    #   4-5 constant
    p  = sep_points
    cdict = {'red':   ((p[0],bc[0],bc[0]),
                       (p[1],bc[0],bc[0]),
                       (p[2],(bc[0]+mc[0])/2.0,(bc[0]+mc[0])/2.0),
                       (p[2],mc[0],mc[0]),
                       (p[3],mc[0],mc[0]),
                       (p[3],(mc[0]+tc[0])/2.0,(mc[0]+tc[0])/2.0),
                       (p[4],tc[0],tc[0]),
                       (p[5],tc[0],tc[0])),

             'green': ((p[0],bc[1],bc[1]),
                       (p[1],bc[1],bc[1]),
                       (p[2],(bc[1]+mc[1])/2.0,(bc[1]+mc[1])/2.0),
                       (p[2],mc[1],mc[1]),
                       (p[3],mc[1],mc[1]),
                       (p[3],(mc[1]+tc[1])/2.0,(mc[1]+tc[1])/2.0),
                       (p[4],tc[1],tc[1]),
                       (p[5],tc[1],tc[1])),

             'blue':  ((p[0],bc[2],bc[2]),
                       (p[1],bc[2],bc[2]),
                       (p[2],(bc[2]+mc[2])/2.0,(bc[2]+mc[2])/2.0),
                       (p[2],mc[2],mc[2]),
                       (p[3],mc[2],mc[2]),
                       (p[3],(mc[2]+tc[2])/2.0,(mc[2]+tc[2])/2.0),
                       (p[4],tc[2],tc[2]),
                       (p[5],tc[2],tc[2])),
    }
    return cdict

############################################################################
# Standard cmap_sep has green on bottom, yellow on top, black separator
p = [0,0.2,0.48,0.52,0.8,1.0]
cdict_sep = create_cdict_for_sep_LSC(green,black,yellow,p)
cmap_sep = LinearSegmentedColormap('UofS separated', cdict_sep)
plt.register_cmap(cmap=cmap_sep)
# reversed sep colormap
cdict_sep_r = create_cdict_for_sep_LSC(yellow,black,green,p)
cmap_sep_r = LinearSegmentedColormap('UofS separated - reversed', cdict_sep_r)
plt.register_cmap(cmap=cmap_sep_r)

############################################################################
# Wider separator
p = [0,0.0,0.3,0.7,1.0,1.0]
cdict_widesep = create_cdict_for_sep_LSC(green,black,yellow,p)
cmap_widesep = LinearSegmentedColormap('UofS separated', cdict_widesep)
plt.register_cmap(cmap=cmap_widesep)
# and reversed
cdict_widesep_r = create_cdict_for_sep_LSC(yellow,black,green,p)
cmap_widesep_r = LinearSegmentedColormap('UofS separated', cdict_widesep_r)
plt.register_cmap(cmap=cmap_widesep_r)

############################################################################
# No separator, though its color still influences the surrounding colors
p = [0,0.0,0.5,0.5,1.0,1.0]
cdict_nosep = create_cdict_for_sep_LSC(green,black,yellow,p)
cmap_nosep = LinearSegmentedColormap('UofS separated', cdict_nosep)
plt.register_cmap(cmap=cmap_nosep)
# and reversed
cdict_nosep_r = create_cdict_for_sep_LSC(yellow,black,green,p)
cmap_nosep_r = LinearSegmentedColormap('UofS separated', cdict_nosep_r)
plt.register_cmap(cmap=cmap_nosep_r)

############################################################################
# Lower position of separator
p = [0,0.0,0.18,0.22,1.0,1.0]
cdict_lowsep = create_cdict_for_sep_LSC(green,black,yellow,p)
cmap_lowsep = LinearSegmentedColormap('UofS separated', cdict_lowsep)
plt.register_cmap(cmap=cmap_lowsep)
# reversed
cdict_lowsep_r = create_cdict_for_sep_LSC(yellow,black,green,p)
cmap_lowsep_r = LinearSegmentedColormap('UofS separated', cdict_lowsep_r)
plt.register_cmap(cmap=cmap_lowsep_r)

############################################################################
# Higher position of separator
p = [0,0.0,0.78,0.82,1.0,1.0]
cdict_highsep = create_cdict_for_sep_LSC(green,black,yellow,p)
cmap_highsep = LinearSegmentedColormap('UofS separated', cdict_highsep)
plt.register_cmap(cmap=cmap_highsep)
# reversed
cdict_highsep_r = create_cdict_for_sep_LSC(yellow,black,green,p)
cmap_highsep_r = LinearSegmentedColormap('UofS separated', cdict_highsep_r)
plt.register_cmap(cmap=cmap_highsep_r)


############################################################################
############################################################################
##   Lines
############################################################################
############################################################################

# lighten up the green for lines in presentations
#   -> (weighted average with white)
nwhite = 1.0
ngreen = 2.0
presen_green = [ (nwhite+ngreen*x)/(nwhite+ngreen) for x in green]
# Set the line cycler to use defined colors!
# plt.rc('axes', prop_cycle=(cycler('color', [presen_green, yellow, black, grey]*4) +
#                            cycler('linestyle', ['-', '-', '--', '-',
#                                                 '--','--','-','--',
#                                                 ':',':',':',':',
#                                                 '-.','-.','-.','-.'])))
