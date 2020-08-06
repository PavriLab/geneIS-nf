#!/usr/bin/env python

import argparse as ap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.lines as mlines
from scipy.interpolate import interpn
from math import ceil
import logging
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# define function to compute the center of each bin and use these x/y coordinated to generate the plot accordingly

greymap = colors.LinearSegmentedColormap.from_list('greymap', ['Grey', 'Grey'], N = 256)
#Args = namedtuple('Args', ['xcol', 'ycol', 'xbins', 'ybins', 'xmin', 'xmax', 'ymin', 'ymax',
#                  'mastertable', 'figwidth', 'figheight', 'plotBackground', 'outputFile',
#                  'subsetColumn', 'colormaps', coarsen])
#args = Args(xcol = 'log2CH12_shLacZ', ycol = 'log2CH12_shMcm6', xbins = 300, ybins = 300, xmin = -1, xmax = 10,
#            ymin = -1, ymax = 10, mastertable = 'CH12_unstimulated_IS.master.tsv', figheight = 8, figwidth = 8,
#            plotBackground = False, subsetColumn = 'class', outputFile = None
#            colormaps = ['White,Blue', 'White,Purple', 'White,Yellow', 'White,Orange', 'White,Green'])

def coarseAverage(h, coarsen = 1):
    h = h.copy()

    if coarsen == 1:
        return h

    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            a = h[i: i + coarsen, j: j + coarsen]

            if not a.any():
                continue

            else:
                h[i: i + coarsen, j: j + coarsen][a != 0] = a[a != 0].mean()

    return h


def mapCenters(h, xbounds, ybounds):
    '''
    computed the centers of each 2d bin defined by xbounds and ybounds
    and returns x and y coordinated of these

    :param xbounds: x-axis boundaries of bins
    :param ybounds: y-axis boundaries of bins

    :return:        numpy.array holding x coordinates of bin centers, numpy.array holding y coordinates of bin centers
    '''
    # computing offsets
    yoffset = ybounds[1] - ybounds[0]
    xoffset = xbounds[1] - xbounds[0]

    tmpxcenters = xbounds[:-1] + xoffset
    tmpycenters = ybounds[:-1] + yoffset

    xcenters, ycenters = [], []

    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            if h[i, j] > 0:
                xcenters.append(tmpxcenters[i])
                ycenters.append(tmpycenters[j])

    return np.array(xcenters), np.array(ycenters)


def aggregateCoordinates(data, xcol, ycol, xbins, ybins, xrange, yrange, subsetcol = None):
    '''
    aggregates given points onto a 2D grid and returns the coordinated of the centers of
    each 2D bin holding at least one point

    :param data:        pandas.DataFrame containing x and y data
    :param xcol:        column of frame to use as x coordinate
    :param ycol:        column of frame to use as y coordinate
    :param xbins:       number of bins on x axis
    :param ybins:       number of bins on y axis
    :param xrange:      tuple giving range of x axis
    :param yrange:      tuple giving range of y axis
    :param subsetcol:   column to use to split data into subsets

    :return:            tuple of numpy.arrays holding x and y coordinates for each filled bin
                        or list of tuple of arrays of length len(data[subsetcol].unique())
    '''
    coords, hs, counts = [], [], []
    if not subsetcol:
        h, xe, ye = np.histogram2d(data[xcol], data[ycol], range=(xrange, yrange), bins=[xbins, ybins])
        coords.append(mapCenters(h, xe, ye))
        hs.append((h, xe, ye))

        counts.append(len(data))

    else:
        cls = sorted(data[subsetcol].unique())
        for c in cls:
            h, xe, ye = np.histogram2d(data[data[subsetcol] == c][xcol], data[data[subsetcol] == c][ycol],
                                       range=(xrange, yrange),
                                       bins=[xbins, ybins])

            coords.append(mapCenters(h, xe, ye))
            hs.append((h, xe, ye))

            counts.append(int(h.sum()))

    return coords, hs, counts


def aggregatePoints(data, xcol, ycol, xbins, ybins, xrange, yrange, subsetcol = None, density = False):
    '''
    aggregates given points onto a 2D grid

    :param data:        pandas.DataFrame containing x and y data
    :param xcol:        column of frame to use as x coordinate
    :param ycol:        column of frame to use as y coordinate
    :param xbins:       number of bins on x axis
    :param ybins:       number of bins on y axis
    :param xrange:      tuple giving range of x axis
    :param yrange:      tuple giving range of y axis
    :param subsetcol:   column to use to split data into subsets

    :return:            list of 2-dimensional arrays holding histogram counts for the given bins
    '''
    hs, counts = [], []
    if not subsetcol:
        h, xe, ye = np.histogram2d(data[xcol], data[ycol],
                                    range=(xrange, yrange),
                                    bins=[xbins, ybins],
                                    density = density)
        hs.append(h)
        counts.append(h.sum())

    else:
        cls = sorted(data[subsetcol].unique())
        for c in cls:
            h, xe, ye = np.histogram2d(data[data[subsetcol] == c][xcol], data[data[subsetcol] == c][ycol],
                                       range=(xrange, yrange),
                                       bins=[xbins, ybins],
                                       density=density)

            hs.append(h)
            counts.append(int(h.sum()))

    return hs, counts


def maskzeros(hs):
    '''
    masks all zeros in the given arrays

    :param hs:  list of 2-dimensional arrays holding data point counts

    :return:    list of 2-dimensional masked arrays holding data point counts
    '''

    maskhs = []
    for h in hs:
        maskhs.append(np.ma.masked_where(h == 0, h, copy = True))

    return maskhs


def plotmesh(hs, clrmaps, xlims, ylims, datalabels, background = None, ax = None, transpose = False,
             xlabel = None, ylabel = None, foldchange = None, grid = True, spacing = 2,
             gridcolor = 'White', legend = True, legendpos = 'best', text = None):
    '''
    visualizes the given numpy arrays. If background is given this is plotted first
    and hs arrays are plotted afterwards in the order given.

    :param hs:          list of 2d arrays containing information to plot
    :param clrmaps:     list of colormaps for the respective arrays in the hs in the same order as hs
    :param ylims:       limits of the y-axis
    :param xlims:       limits of the x-axis
    :param datalabels:  labels to use for each subset for the legend
    :param background:  2d array containing the background to plots
    :param ax:          optional matplotlib.Axes object to generate the plot in
    :param transpose:   if True, transposes the arrays before plotting
                        this is mainly due to flipped plotting of array by the
                        pcolormesh method
    :param xlabel:      label of the x-axis
    :param ylabel:      label of the y-axis
    :param foldchange:  value used to plot threshold indicators
    :param grid:        if True, also plots a grid
    :param spacing:     integerstep to determine the spacing between ticks
                        e.g. a spacing of 2 generates ticks at even positions starting from 0
    :param gridcolor:   color of the grid
    :param legend:      if True, inserts a legend into the plot
    :param legendpos:   specifies the position of the legend if legend is True
    :param text:        if given, is written into the lower left right corner of the plot

    :return:            fig, ax if ax is not given else just None, ax
    '''

    if not ax:
        fig, ax = plt.subplots()

    else:
        fig = None

    if type(background) == np.ndarray:
        ax.pcolormesh(background.T if transpose else background, cmap = greymap, zorder = 2)

    handles = []
    for h, cmap, label in zip(hs, clrmaps, datalabels):
        if legend:
            handles.append(mlines.Line2D([], [], color = cmap(1), marker = 's', label = label, ls = ''))

        ax.pcolormesh(h.T if transpose else h, cmap = cmap, zorder = 2)

    # computing tickpositions
    steps, lims = {}, {}
    for (min_, max_), setpos, setlabels, getlims, k in zip([xlims, ylims], [ax.set_xticks, ax.set_yticks],
                                                           [ax.set_xticklabels, ax.set_yticklabels],
                                                           [ax.get_xlim, ax.get_ylim], ['x', 'y']):
        axismin, axismax = getlims()
        lims[k + 'lim'] = (axismin, axismax)

        step = axismax/(max_ - min_)
        steps[k] = step
        integerscale = axismax - step * (max_ - min_ - int(max_ - min_))
        spacer = integerscale/int(max_ - min_)

        if ceil(min_)//2 != 0:
            startpos = (ceil(min_) - min_) * step + spacer
            startlabel = ceil(min_) + 1

        else:
            startpos = (ceil(min_) - min_) * step
            startlabel = ceil(min_)
        
        setpos(np.arange(startpos, axismax, spacer * spacing))
        setlabels(range(startlabel, int(max_), spacing))

        if grid:
            ax.grid(ls = '-', zorder = 1, color = gridcolor)

    if foldchange:
        ax.plot(lims['xlim'], lims['ylim'], ls = '--', color = 'black', lw = 1, zorder = 3)
        ax.plot(lims['xlim'], [lim + foldchange * steps['y'] for lim in lims['ylim']],
                ls='--', color='grey', lw=1, zorder = 3)
        ax.plot(lims['xlim'], [lim - foldchange * steps['y'] for lim in lims['ylim']],
                ls='--', color='grey', lw=1, zorder = 3)

    ax.set(**lims)

    if xlabel:
        ax.set_xlabel(xlabel)

    if ylabel:
        ax.set_ylabel(ylabel)

    if legend:
        ax.legend(handles = handles, loc = legendpos, frameon = False)

    if text:
        ax.text(lims['xlim'][1], lims['ylim'][0], text, va = 'bottom', ha = 'right')

    return fig, ax


def scatterplot(coords, hs, clrs, density, datalabels, bkgrdcoords = None, ax = None, xlabel = None,
                ylabel = None, ylims = None, xlims = None, foldchange = None, coarsen = 1, grid = True,
                legend = True, legendpos = 'best', text = None):
    '''
    visualizes the given numpy arrays. If background is given this is plotted first
    and hs arrays are plotted afterwards in the order given.
    If density is true the histogram countinformation is used to interpolate
    intensity values for the colormaps given via clrs

    :param coords:      list of tuples of x and y coordinates of filled bin centers
    :param hs:          list of 2d-histograms and their bin boundaries in x and y
    :param clrs:        list of colormaps or colors for the respective arrays in the coords in the same order as hs
    :param density:     list of boolean values indicating the use of histogram2d count
                        information of coords arrays for colormap intensities (note that if true is set
                        for a given choord array clrs must contain a colormap at the same position)
    :param datalabels:  labels to use for each subset for the legend
    :param background:  2d array containing the background to plots
    :param ax:          optional matplotlib.Axes object to generate the plot in
    :param xlabel:      label of the x-axis
    :param ylabel:      label of the y-axis
    :param ylims:       limits of the y-axis
    :param xlims:       limits of the x-axis
    :param foldchange:  used to plot threshold lines on both sides of the diagonal
    :param coarsen:     number of neighbouring pixels to consider for windowaveraging
    :param grid:        if True, adds a grid to the plot
    :param legend:      if True, inserts a legend into the plot
    :param legendpos:   specifies the position of the legend if legend is True
    :param text:        if given, is written into the lower left right corner of the plot


    :return:            fig, ax if ax is not given else just None, ax
    '''

    if not ax:
        fig, ax = plt.subplots()

    else:
        fig = None

    if grid:
        ax.grid(ls='-', zorder = 1)

    if bkgrdcoords:
        ax.scatter(bkgrdcoords[0], bkgrdcoords[1], c = 'grey', marker = '.', zorder = 2)

    for (x, y), (h, xe, ye), clr, d, label in zip(coords, hs, clrs, density, datalabels):
        if d:
            if coarsen > 1:
                h = coarseAverage(h, coarsen = coarsen)

            z = interpn((0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])), h, np.vstack([x, y]).T,
                        method='splinef2d',
                        bounds_error=False)

            ax.scatter(x, y, c = z, cmap = clr, marker = '.', zorder = 2, label = label)

        else:
            ax.scatter(x, y, c = clr, marker = '.', zorder = 2, label = label)

    if ylabel:
        ax.set_ylabel(ylabel)

    if xlabel:
        ax.set_xlabel(xlabel)

    for limtype, lims in zip(['xlim', 'ylim'], [xlims, ylims]):
        if lims:
            limdict = {limtype:lims}
            ax.set(**limdict)

    if foldchange:
        ax.plot(xlims, ylims, ls='--', color='black', lw=1, zorder = 3)
        ax.plot(xlims, [lim + foldchange for lim in ylims], ls='--', color='grey', lw=1, zorder = 3)
        ax.plot(xlims, [lim - foldchange for lim in ylims], ls='--', color='grey', lw=1, zorder = 3)

    if legend:
        ax.legend(loc = legendpos, frameon = False)

    if text:
        ax.text(xlims[1], ylims[0], text, va = 'bottom', ha = 'right')

    return fig, ax


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-mt', '--mastertable', required = True,
                    help = 'mastertable containing information to plot')
parser.add_argument('--xcol', required = True,
                    help = 'column of mastertable holding x coordinates')
parser.add_argument('--ycol', required = True,
                    help = 'column of mastertable holding y coordinates')
parser.add_argument('--labels', nargs = '+',
                    help = 'labels for the legend. Has to be in the same order as sorted categories of subsetcol')
parser.add_argument('--xbins', default = 300, type = int,
                    help = 'number of bins on x-axis for aggregation')
parser.add_argument('--ybins', default = 300, type = int,
                    help = 'number of bins on y-axis for aggregation')
parser.add_argument('--xmin', default = -1.0, type = float,
                    help = 'minimum value of x-axis')
parser.add_argument('--xmax', default = 10.0, type = float,
                    help = 'maximum value of x-axis')
parser.add_argument('--ymin', default = -1.0, type = float,
                    help = 'minimum value of y-axis')
parser.add_argument('--ymax', default = 10.0, type = float,
                    help = 'maximum value of y-axis')
parser.add_argument('--subsetColumn', default = None,
                    help = 'column in mastertable to use to determine subsets')
parser.add_argument('--xlabel', default = None,
                    help = 'x-axis label')
parser.add_argument('--ylabel', default = None,
                    help = 'y-axis label')
parser.add_argument('--coarsen', default = 3, type = int,
                    help = '''determined the neighbouring pixels to use for average computation
                              if density is True. set to 1 to turn it off''')
parser.add_argument('--usecounts', default = False, action = 'store_true',
                    help = '''if --legend is set, setting --usecounts displays counts of each
                              category in subsetColumn in the legend''')
parser.add_argument('--density', nargs = '+',
                    help = '''space separated list of boolean values in the same order as sorted categorical
                              values of subsetColumn. Indicates the use of histogram count data as colormap intensities
                              for plotting. (note that if True is set for a given position the corresponding position
                              in the colormaps has to be a valid colormap setting. Can either be T/True or F/False''')
parser.add_argument('--colormaps', nargs = '+',
                    help = '''space separated list of colormaps or colors to use. Can either be a list
                              of already existing colormaps or a list of comma-separated colors
                              that are then used to generate custom linear segmented colormaps with
                              N = 256 (e.g. White,Red Blue,Red will generate two colormaps one that
                              transitions from white to red and one from blue to red)
                              This has to be in the same order as the sorted categorical values of
                              subsetColumn''')
parser.add_argument('--foldchange', default = 0, type = float,
                    help = 'threshold value to determine cutoff lines in plot as log2 value')
parser.add_argument('--plotmethod', default = 'scatter', choices = ['scatter', 'mesh'],
                    help = '''determines which plotting scheme to use, scatter creates a scatterplot,
                              mesh uses the pcolormesh method (note that when using mesh all passed
                              colors have to be valid colormaps''')
parser.add_argument('--stylesheet', default = None,
                    help = 'stylesheet to use for plotting. See available at matplotlib')
parser.add_argument('--legend', default = False, action = 'store_true',
                    help = 'if set, adds a legend to the plot')
parser.add_argument('--legendposition', default = 'best',
                    choices = ['best', 'upper right', 'upper left', 'lower left', 'lower right', 'right',
                               'center left', 'center right', 'lower center', 'upper center', 'center'],
                    help = 'determines the position of the legend if --legend is set')
parser.add_argument('-o', '--outputFile', required = True,
                    help = 'file to save the plot in (must be a valid format supported by matplotlib)')
parser.add_argument('--figwidth', default = 8.0, type = float,
                    help = 'width of the generated figure')
parser.add_argument('--figheight', default = 8.0, type = float,
                    help = 'height of the generated figure')
parser.add_argument('--text', default = None,
                    help = 'text to write into the lower right corner of the plot')
parser.add_argument('--plotCounts', default = False, action = 'store_true',
                    help = '''if set plots counts of overlaps and total peaks into
                              lower right corner. --text overrides this flag''')
args = parser.parse_args()

# reading inputtable
table = pd.read_csv(args.mastertable, sep = '\t')

# reading colors and colormaps
clrs = []
for i, c in enumerate(args.colormaps):
    if len(c.split(',')) > 1:
        clrs.append(colors.LinearSegmentedColormap.from_list('custom' + str(i), c.split(','), N = 256))

    else:
        clrs.append(c)

# reading ranges
xrange = (args.xmin, args.xmax)
yrange = (args.ymin, args.ymax)

if args.stylesheet:
    plt.style.use(args.stylesheet)
    grid = False if args.plotmethod == 'scatter' else True
    gridcolor = 'White'

else:
    grid = True
    gridcolor = 'Grey'

fig, ax = plt.subplots()
if args.plotmethod == 'scatter':
    logging.info('computing 2d histograms and pixel centers')
    coords, hs, counts = aggregateCoordinates(table, args.xcol, args.ycol, args.xbins, args.ybins, xrange, yrange,
                                              subsetcol = args.subsetColumn)

    if args.usecounts:
        labels = [' '.join([l, '(' + str(c) + ')']) for l, c in zip(args.labels, counts)]

    else:
        labels = args.labels

    if args.text:
        text = args.text

    elif args.plotCounts:
        text = "{0}/{1}".format(counts[0], len(table))

    else:
        text = None

    logging.info('generating figure')
    tmp, ax = scatterplot(coords, hs, clrs, args.density, labels,
                          xlabel = args.xlabel,
                          ylabel = args.ylabel,
                          ylims = yrange,
                          xlims = xrange,
                          ax = ax,
                          foldchange = args.foldchange,
                          coarsen = args.coarsen,
                          grid = grid, 
                          text = text,
                          legend = args.legend,
                          legendpos = args.legendposition)

else:
    logging.info('aggregating points')
    hs, counts = aggregatePoints(table, args.xcol, args.ycol, args.xbins, args.ybins, xrange, yrange,
                                 subsetcol = args.subsetColumn)

    if args.usecounts:
        labels = [' '.join([l, '(' + str(c) + ')']) for l, c in zip(args.labels, counts)]

    else:
        labels = args.labels

    # if density is given for a particular subset
    # we compute the coarened values
    for i, h, d in zip(range(len(hs)), hs, args.density):
        if args.coarsen > 1:
            hs[i] = coarseAverage(h, args.coarsen)

    if args.text:
        text = args.text

    elif args.plotCounts:
        text = "{0}/{1}".format(counts[0], len(table))

    else:
        text = None

    logging.info('masking arrays')
    hs = maskzeros(hs)

    logging.info('generating figure')
    ax = plotmesh(hs, clrs, xrange, yrange, labels,
                  xlabel = args.xlabel,
                  ylabel = args.ylabel,
                  foldchange = args.foldchange,
                  grid = grid,
                  gridcolor = gridcolor,
                  ax = ax, transpose = True,
                  legend = args.legend,
                  legendpos = args.legendposition,
                  text = text)

fig.set_figwidth(args.figwidth)
fig.set_figheight(args.figheight)
fig.tight_layout()
fig.savefig(args.outputFile)
