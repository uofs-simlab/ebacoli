"""
Python data structures and methods to be used with eBACOLI data.
"""
import colors as c

import sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from scipy.interpolate import splev
from scipy.integrate import quadrature

NOTHING = object()

class EbacoliData:
    """
    Structure for storing eBACOLI output:
    """
    def __init__(self,fname=None):
        if fname:
            self.read_bspline_file(fname)
        else:
            self.has_data = False

    def read_bspline_file(self,fname=None):
        """Reads eBACOLI bspline output files."""
        self.filename = fname
        with open(fname) as f:
            lines     = f.readlines()
            firstline = lines[0].split(',')
            self.npde = int(firstline[0])
            self.nu   = int(firstline[1])
            self.nint = int(firstline[2])
            self.time = float(firstline[3])
            self.atol = float(firstline[4])
            self.rtol = float(firstline[5])
            self.xbs  = np.fromstring(lines[1], sep=' ')
            y         = np.fromstring(lines[2], sep=' ').reshape((self.npde,-1))
            self.p    = int(lines[3])

        # Derived data
        self.nv = self.npde - self.nu
        self.u_coeff = y[:self.nu,:]
        self.v_coeff = y[self.nu:,:]
        self.xrange = [min(self.xbs),max(self.xbs)]

        self.has_data = True

    def clear_data(self):
        """Nullify object's data."""
        self.npde     = None
        self.nu       = None
        self.nv       = None
        self.nint     = None
        self.time     = None
        self.xbs      = None
        self.p        = None
        self.u_coeff  = None
        self.v_coeff  = None
        self.xrange   = None
        self.has_data = False

    def u_at_points(self, x, i):
        """
        Returns u[i] value at input points x.  NOTE that if i > nu,
        this will return v[j], where j = i-nu.
        """
        if i < self.nu:
            u = splev(x, (self.xbs,self.u_coeff[i],self.p))
        else:
            u = splev(x, (self.xbs,self.v_coeff[i-self.nu],self.p))
        return u

    def v_at_points(self, x, i):
        """ Returns v[i] at input points x."""
        v = splev(x, (self.xbs,self.v_coeff[i],self.p))
        return v

    def continuous_representation(self,n_cont=1000,u_num=NOTHING,v_num=NOTHING):
        """
        Returns a 'continuous' representation of the data on self.xrange

        INPUT:
          n_cont  - number of points in the 'continuous' independent variable
          [u_num] - list of which u fields to return
          [v_num] - list of which v fields to return

        OUTPUT:
          u_list  - list of u continuous variables
          v_list  - list of v continuous variables
        """
        # Default behaviour is to return all fields
        if u_num is NOTHING:
            u_num = range(self.nu)
        if v_num is NOTHING:
            v_num = range(self.nv)

        u_list = []
        v_list = []

        # continuous x variable
        x = np.linspace(self.xrange[0],self.xrange[1],n_cont)

        for i in u_num:
            u_list.append(splev(x, (self.xbs,self.u_coeff[i],self.p)))
        for i in v_num:
            v_list.append(splev(x, (self.xbs,self.v_coeff[i],self.p)))

        return u_list, v_list

## END OF class EbacoliData

############################################################################
############################################################################

#### Functions that use EbacoliData objects

def plot_heatmap_of_field(ed_list, field_str, field_num, cmap=c.cmap_sep,
                          title=NOTHING, extension="png"):
    """
    Plots the heatmap of a field from a list of EbacoliData structures.

    INPUT:
      ed_list   - list of EbacoliData structures, assumed to be ordered in time,
                  with each element having the same xranges, first element
                  considered to be at time 0
      field_str - string for the field, must be either "u" or "v"
      field_num - which component of u or v to plot
      [cmap]    - optional colormap
    """

    # Default title is field name
    if title is NOTHING:
        title = "$%s_{%d}(x,t)$" % (field_str, field_num+1)

    # lambda functions for getting the right coeffs
    if field_str == "u":
        field_coeff = lambda ed: ed.u_coeff[field_num]
    elif field_str == "v":
        field_coeff = lambda ed: ed.v_coeff[field_num]
    else:
        error("Must use 'u' or 'v' as first argument")
        sys.exit(1) # can probably be handled more nicely with an exception

    final_time = ed_list[-1].time
    n_ed = len(ed_list)

    # square grid for space and time
    x = np.linspace(ed_list[0].xrange[0],ed_list[0].xrange[1],n_ed)
    t = np.linspace(0.0,final_time,n_ed)

    tr_x,tr_t = np.meshgrid(x,t)
    z = np.zeros([n_ed,n_ed])

    plt.figure()

    # pass through the data, now reading into our square grid
    counter = -1
    for ed in ed_list:
        counter = counter+1
        z[:,counter] = splev(x, (ed.xbs,field_coeff(ed),ed.p))

    # Plot the heatmap.
    # NOTES: transpose to get the x on the horiz axis
    #        origin set to lower
    span = np.max(z) - np.min(z)
    low = np.min(z)
    high = np.min(z) + span
    plt.imshow(np.transpose(z),cmap=cmap,
               vmax=high,vmin=low,
               origin='lower',extent=[ed.xrange[0],ed.xrange[1],0.0,final_time],aspect="auto")
    # plt.contourf(tr_x,tr_t,np.transpose(z),cmap=cmap)
    plt.colorbar()

    plt.xlabel('$x$')
    plt.ylabel('$t$')
    plt.title(title)

    plt.tight_layout()

    plt.savefig("heatmap_%s%d." % (field_str, field_num+1) + extension)
    # plt.show()
    plt.gcf().clear()

def plot_interpolated_fields(ed, u_num=NOTHING, v_num=NOTHING,
                             u_labels=NOTHING, v_labels=NOTHING,
                             u_linestyles=NOTHING, v_linestyles=NOTHING,
                             u_colors=NOTHING, v_colors=NOTHING, n_cont=1000,
                             knot_show=False, knot_bottom=None):
    """
    Plot interpolated fields for eBACOLI data.

    INPUT:
      ed             - EbacoliData struct
      [u_num]        - list of u fields to plot
      [v_num]        - list of v fields to plot
      [u_labels]     - list of u field names for the legend
      [v_labels]     - list of v field names for the legend
      [u_linestyles] - list of u linestyles
      [v_linestyles] - list of v linestyles
      [u_colors]     - list of colors to plot u fields
      [v_colors]     - list of colors to plot v fields
      [n_cont]       - size of 'continuous' variable
      [knot_show]    - show the knots (currently from knot_bottom to min plotted value)
      [knot_bottom]  - how far down do the knots' lines go
    """
    # Default behaviour is to plot all fields
    if u_num is NOTHING:
        u_num = range(self.nu)
    if v_num is NOTHING:
        v_num = range(self.nv)

    # Default field labels
    if u_labels is NOTHING:
        u_labels = ["$u_{%d}$" % (i+1) for i in u_num]
    if v_labels is NOTHING:
        v_labels = ["$v_{%d}$" % (i+1) for i in v_num]

    assert(len(u_num) == len(u_labels))
    assert(len(v_num) == len(v_labels))

    # default linestyles
    if u_linestyles is NOTHING:
        u_linestyles = ["-" for i in u_num]
    if v_linestyles is NOTHING:
        v_linestyles = ["--" for i in v_num]

    assert(len(u_num) == len(u_linestyles))
    assert(len(v_num) == len(v_linestyles))

    # ASSERT that there are enough colors for all u_num
    if u_colors is NOTHING:
        u_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][:len(u_num)]
    else:
        assert(len(u_colors) == len(u_num))
    # and v_num
    if v_colors is NOTHING:
        v_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][
            len(u_num):len(u_num)+len(v_num) ]
    else:
        assert(len(v_colors) == len(v_num))

    ########################################################################
    ## Plotting Lines here
    ########################################################################

    # 'continuous' x variable
    x = np.linspace(ed.xrange[0], ed.xrange[1], n_cont)
    # plot fields

    ii = -1
    for i in u_num:
        ii += 1
        u = splev(x, (ed.xbs,ed.u_coeff[i],ed.p))
        plt.plot(x, u, label=u_labels[ii], linestyle=u_linestyles[ii], color=u_colors[ii])

    ii = -1
    for i in v_num:
        ii += 1
        v = splev(x, (ed.xbs,ed.v_coeff[i],ed.p))
        plt.plot(x, v, label=v_labels[ii], linestyle=v_linestyles[ii], color=v_colors[ii])

    # axes manipulations
    axes = plt.gca()
    axes.set_xlim(ed.xrange[0],ed.xrange[1])
    [fig_bot, fig_top] = axes.get_ylim()
    [fig_left, fig_right] = axes.get_xlim()

    # Show the knot locations
    if knot_show:
        for i in range(len(ed.xbs)):
            uv_at_knot = []
            for j in u_num:
                uv_at_knot.append(splev(ed.xbs[i],(ed.xbs,ed.u_coeff[j],ed.p)))
            for j in v_num:
                uv_at_knot.append(splev(ed.xbs[i],(ed.xbs,ed.v_coeff[j],ed.p)))
            # draw line only up to PDE component with smallest value
            knot_top = min(uv_at_knot)
            plt.plot([ed.xbs[i],ed.xbs[i]],[knot_bottom,knot_top],'k-',linewidth=0.5)

        # reset axes to where they were before knots drawn
        axes.set_ylim([fig_bot,fig_top])
        axes.set_xlim([fig_left,fig_right])

    return axes

def plot_interpolated_error(ed, true_solution, u_num=NOTHING,
                            v_num=NOTHING, ylim=NOTHING, xlim=NOTHING,
                            u_colors=None, v_colors=None, n_cont=1000,
                            knot_show=False):
    """
    Plot the difference between an eBACOLI solution and a true solution over xrange.

    INPUT:
      ed            - EbacoliData structure
      true_solution - function of one variable, must return 2 lists, 1 of size
                      ed.nu and 1 of size ed.nv
      [u_num]       - list of u fields to plot
      [v_num]       - list of v fields to plot
      true_solution - function of one variable that returns
      [u_colors]    - list of colors to plot u fields
      [v_colors]    - list of colors to plot v fields
      [n_cont]      - size of 'continuous' variable
      [knot_show]   - show the knots (currently from knot_bottom to min plotted value)
    """
    # Default behaviour is to plot all fields
    if u_num is NOTHING:
        u_num = range(ed.nu)
    if v_num is NOTHING:
        v_num = range(ed.nv)

    # 'continuous' x variable
    x = np.linspace(ed.xrange[0], ed.xrange[1], 1000)
    u,v = ed.continuous_representation(1000)
    sol_u, sol_v = true_solution(x)

    for i in u_num:
        plt.plot(x, sol_u[i]-u[i], label='$u_{%d}$'%(i+1),linestyle="-")
    for i in v_num:
        plt.plot(x, sol_v[i]-v[i], label='$v_{%d}$'%(i+1),linestyle="-")

    axes = plt.gca()
    if ylim is not NOTHING:
        axes.set_ylim(ylim)
    if xlim is not NOTHING:
        axes.set_xlim(xlim)

    [fig_bot, fig_top] = axes.get_ylim()
    [fig_left, fig_right] = axes.get_xlim()

    # Show the knot locations
    if knot_show:
        for i in range(len(ed.xbs)):
            uv_at_knot = []
            sol_u, sol_v = true_solution(ed.xbs[i])
            for j in range(ed.nu):
                u = ed.u_at_points(ed.xbs[i],j)
                uv_at_knot.append(sol_u[j]-u)
            for j in range(ed.nv):
                v = ed.v_at_points(ed.xbs[i],j)
                uv_at_knot.append(sol_v[j]-v)
            # draw line between max and min values
            plt.plot([ed.xbs[i],ed.xbs[i]],[np.max(uv_at_knot),fig_top],'k-',linewidth=0.5)
            plt.plot([ed.xbs[i],ed.xbs[i]],[fig_bot,np.min(uv_at_knot)],'k-',linewidth=0.5)

        # reset axes to where they were before knots drawn
        axes.set_ylim([fig_bot,fig_top])
        axes.set_xlim([fig_left,fig_right])

    return axes

def compute_error_table_dict(ed_list, true_solution):
    """
    Returns a DataFrame of error table columns.

    INPUT:
      ed_list       - list of EbacoliData structures
      true_solution - function of one variable, must return a list of size ed.npde

    OUTPUT:
      pandas DataFrame with keys: 'atol', 'rtol', 'err', 'nint'
    """

    err_array = []
    atol_array = []
    rtol_array = []
    nint_array = []

    for ed in ed_list:

        err = compute_error(ed, true_solution)
        err_array.append(err)

        atol_array.append(ed.atol)
        rtol_array.append(ed.rtol)
        nint_array.append(ed.nint)

    return pd.DataFrame({'atol': atol_array, 'rtol': rtol_array,
                         'err': err_array, 'nint': nint_array})

def compute_error(ed, true_solution):
    """
    Error computed as
    $$
    E = \sqrt{\sum_k^{npde}\int_{xR}^{xL}\left(Y_k(x,1)-y_k(x,1)\right)^2\mathrm{d}x}
    $$
    by Gaussian quadrature (scipy.integrate.quadrature)

    INPUT:
      ed            - an EbacoliData structure
      true_solution - function of one variable, must return a list of size ed.npde

    OUTPUT:
      err    - error computed as above
    """

    err = 0.0
    for i in range(ed.npde):
        integrand = lambda x: (ed.u_at_points(x,i)-true_solution(x)[i])**2
        [q, quad_err] = quadrature(integrand, ed.xrange[0], ed.xrange[1],
                                   args=(), tol=1e-12, rtol=1e-12, maxiter=100)
        # print "Integral computed to %g" % quad_err
        err += q

    err = np.sqrt(err)
    return err

def print_latex_error_table(error_table,filename=NOTHING):
    """
    Writes an error table from a dictionary computed by
    compute_error_table_dict.

    INPUT:
      error_table - pandas DataFrame holding error table columns
      [filename]  - filename to write to (sys.stdout is default)
    """

    # redirect stdout to file
    if filename is not NOTHING:
        stdout = sys.stdout
        sys.stdout = open(filename, "w")

    # write conents
    print "\\begin{tabular}{r|r|r|r}"
    print error_table.to_csv(
        columns=["atol","rtol","err","nint"],
        header=["{\\tt ATOL}", "{\\tt RTOL}", "$E$", "$N$"],
        sep="&", line_terminator="\\\\\n", float_format="%.1e",
        index=False
    ).replace("\n","\n\hline\n",1
    ).replace("&"," & "
    ).replace("\\\\"," \\\\"
    ).replace("\"",""
    )[:-4]  # slice off last newline and LaTeX newline
    print "\\end{tabular}"

    # restore stdout
    if filename is not NOTHING:
        sys.stdout = stdout

def print_csv_error_table(error_table,filename=NOTHING):
    """
    Writes an error table from a dictionary computed by
    compute_error_table_dict.

    INPUT:
      error_table - pandas DataFrame holding error table columns
      [filename]  - filename to write to (sys.stdout is default)
    """

    # redirect stdout to file
    if filename is not NOTHING:
        stdout = sys.stdout
        sys.stdout = open(filename, "w")

    error_table.to_csv(sys.stdout,index=False)

    # restore stdout
    if filename is not NOTHING:
        sys.stdout = stdout

def read_csv_error_table(filename):
    """
    Reads a csv error table from file.

    INPUT:
      filename    - name of csv file to read

    OUTPUT:
      pandas DataFrame of the table
    """

    error_table = pd.read_csv(filename)

    return error_table
