
#-------------------------------------------------------------------------------
# function definition of Plotting utilities
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

import numpy as np

import matplotlib.pyplot as plt

def set_imshow_ticks_(axe, arr, axis, points=None, fmt='%1.3f', rot=0, fontsize=None):
    r"""
    customize ticks and ticklabels for a specific axe.

    Parameters
    ----------
    axe : matplotlib.pyplot.Axes
        the axe whose ticks and ticklabels we are going to modify.

    arr : array-like, numpy.ndarray,
        the array to create ticklabels.

    axis : 'x' or 'y'
        modify ticks and ticklabels of 'x' or 'y' axis.

    points : int or list of int or None, optional
        int : nbins of ticks
        list of int : list of ticks
        None : use the default ticks
        default : None

    fmt : string Formattor
        format to for string formatting ticklabels.

    rot : angle, [:math:`^\circ`]
        rotation angle of ticklabels.


    Returns
    -------
    None

    """

    if axis not in ('x', 'y'):
        raise ValueError(f"`axis` argument should either be 'x' or 'y'")

    #-- if integer format requied, convert array to np.int64
    if fmt[-1] == 'd':
        arr_ = arr.astype(np.int64)
    else:
        arr_ = arr.copy()

    #-- if points is integer, create equally spaced points
    if isinstance(points, int):
        temp = np.linspace(0, arr_.shape[0]-1, points)
        #points = ((temp[1:] + temp[:-1]) * 0.5).astype(np.int64)
        points = temp.astype(np.int64)
    #-- if points is None, use the default ticks
    elif points is None:
        if axis in ('x',):
            points = axe.get_xticks()[:-1].astype(np.int64)
        elif axis in ('y',):
            points = axe.get_yticks()[:-1].astype(np.int64)
        points = points[points>=0]

    #-- format ticklabels
    ticklabels = [("{:" + f"{fmt[1:]}" + "}").format(arr_[i]) for i in points]

    #-- set ticks and ticklabels
    if axis in ('x',):
        axe.set_xticks(points)
        axe.set_xticklabels(ticklabels, rotation=rot, fontsize=fontsize)
    elif axis in ('y',):
        axe.set_yticks(points)
        axe.set_yticklabels(ticklabels, rotation=rot, fontsize=fontsize)

    return None


def remove_tick_ticklabel_(*args, kind="xy"):
    r"""
    turn off x/y ticks and ticklables.

    4 |                |
    3 |                |
    2 |          --->  |
    1 |                |
      +------+         +-------+
      0   1   2
    """
    if 'x' in kind:
        for ax in args:
            ax.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                labelbottom=False,)
    if 'y' in kind:
        for ax in args:
            ax.tick_params(
                axis='y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                left=False,         # ticks along the left edge are off
                labelleft=False)

def remove_spline_(*args,pos=("left","right","top","bottom")):
    r"""
    turn off the spline of all axes(*args)
    """

    for p in pos:
        assert p in ("left","right","top","bottom"), "bad position argument"

        for ax in args:
            ax.spines[p].set_visible(False)


def axes_no_padding_(fig_kw={"figsize":(8,4),"dpi":100}, axe_kw={"ax1":[0,0,1,1]}):
    r""" """
    fig = plt.figure(figsize=fig_kw["figsize"], dpi=fig_kw["dpi"])

    axe_dict = {}
    for key, val in axe_kw.items():
        ax_ = fig.add_axes(val)
        axe_dict[key] = ax_
        remove_tick_ticklabel(ax_, kind="xy")
        remove_spline(ax_, pos=("left","right","top","bottom"))

    return fig, axe_dict

