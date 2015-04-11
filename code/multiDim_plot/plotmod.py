"""
Author : Karen Ng
This program takes the output pickle arrays of TSM.py and creates some plots
and statistics.  This is largely based on PlotTSM.py or PlotTSM_array.py by
Will Dawson.

The main function of this program is to generate report quality plots,
especially a covariance array plot

usage: if you use any part of this code for your publication, please cite Ng et
al. 2015 (a.k.a. http://arxiv.org/abs/1412.1826v1)
"""
from __future__ import division
import pylab
import numpy as np
import numpy
import pickle
import warnings
from astrostats import biweightLoc, bcpcl
from astropy.stats import biweight_location as C_BI
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
warnings.filterwarnings("ignore", category=DeprecationWarning)
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from astroML.plotting import hist as astroMLhist

from astroML import density_estimation as de
from statsmodels.api import nonparametric
from types import ModuleType

# -------- for masking and subsetting --------------------------------------


def loadcombo(prefix, index, suffix):
    """loads the data from pickle files
    Parameters:
    ===========
    prefix = string
        denotes the name of the path to the file
    index = number
        denotes the index to add to the file name
    filename consists of prefix + suffix+.pickle see below

    Returns
    =======
    numpy array
        data
    """
    array = []
    for i in index:
        filename = prefix + i + '_' + suffix + '.pickle'
        # read in the pickled array
        F = open(filename)
        tmp = pickle.load(F)
        F.close()
        array = numpy.append(array, tmp)

    #filename = prefix+suffix+'.pickle'
    # read in the pickled array
    #F = open(filename)
    #tmp = pickle.load(F)
    # F.close()
    #array = numpy.append(array,tmp)
    return array


def load_pickles_to_df(par, prefix, index, msun1e14=True,
                       msun_string=['m_1', 'm_2'], verbose=True,
                       save=False, path=None, key="df", mode="w",
                       complib="blosc"):
    """
    Parameters
    =========
    par = list of strings
        denotes the name of the parameters
    prefix = string
        prefix that denotes the path and file name excluding extension
        .pickle
    index = number
        how many pickle files there are
    msun1e14 = logical
        if we want units of mass to be in 1e14 M_sun
    msun_string = list of strings
        the strings are the keys to the df for changing the units
    verbose: logical
    save = logical, will save the data to a hdf5 file if true
    path = str, full path of the hdf5 including name to be saved
    key = str, corresponds to the key for the hdf5 file
    complib= str, compression format

    Returns
    ======
    Pandas Dataframe that have all nans removed

    """
    for i in range(len(par)):
        if verbose:
            print 'loading ' + par[i]
        d = loadcombo(prefix, index, par[i])
        if i == 0:
            data = pd.DataFrame(d, columns=[par[i]])
        else:
            data[par[i]] = pd.DataFrame(d)

    if verbose:
        print "dropping NA! Original data size = {0}".format(data.shape[0])
    data = data.dropna()

    if msun1e14 is True:
        for mass in msun_string:
            if verbose:
                print "converting entry " + mass + \
                    " to units of 1e14 m_sun"
            data[mass] = data[mass] / (1e14)

    if save is True:
        if path is None:
            path = input(
                "Please enter file path as a string inside quotes")

        print "New data size = {0}".format(data.shape[0])
        print "overwriting dataframe to ", path
        data.to_hdf(path, key=key, mode=mode, complib=complib)

    return data


def mask_bigger_than_age_of_universe(
        df, z, H0, Om0, T=None, TSM_0=None, TSM_1=None):
    """assumes a FlatLambdaCDM cosmology to calculate age of universe at
    particular redshift and returns a mask that masks time values
    that are bigger than the age of the universe or if they are undefined

    Parameters
    =========
    df = pandas dataframe - outputs from load_pickles_to_df
    z = float - redshift
    H0 = Hubble parameter
    Om0 = Relative density of matter
    T = string - denotes column name for T in the dataframe
    TSM_0 = string - denotes column name for TSM_0 in the dataframe
    TSM_1 = string - denotes column name for TSM_1 in the dataframe

    Returns
    ======
    combined mask = numpy array
    age = float - age of the universe at corresponding redshift in Gyrs
    """
    assert T is not None or TSM_0 is not None or TSM_1 is not None, \
        "no relevant col names for T, TSM_0 and TSM_1 specified"

    cosmology = FlatLambdaCDM(H0=H0, Om0=Om0)
    age_of_universe = cosmology.age(z)
    age = age_of_universe.value

    masks = {}
    if TSM_0 is not None:
        masks[TSM_0] = df[TSM_0] < age
    if TSM_1 is not None:
        masks[TSM_1] = df[TSM_1] < age
    if T is not None:
        masks[T] = np.isfinite(df[T])

    mask = np.ones(df.shape[0])
    mask = np.logical_and(masks[TSM_0], masks[T])
    print "# of masked rows = {0}".format(df.shape[0] - np.sum(mask))
    print "# of remaining rows = {0}".format(np.sum(mask))
    print "% of original data remaining = {0:.2f}".format(
        (np.sum(mask)) / df.shape[0] * 100)

    return mask, masks[TSM_1], age


def radio_dist_prior(d_3D, d_3Dmax=3.0, d_3Dmin=1.0):
    '''
    Stability: to be tested
    input:
    d_3D = numpy array to be masked, in unit of Mpc
    d_3Dmax = float, the upper limit to be masked out, in unit of Mpc
    d_3Dmin = float, the lower limit to be masked out, in unit of Mpc
    output:
    mask = numpy array that gives 1 if it is within the range
            0 if it is NOT within the specified range
    count = number of entries along the array that has value 1
    '''
    mask = np.logical_and(d_3D < d_3Dmax, d_3D >= d_3Dmin)
    count = np.sum(mask)
    return mask, count


def radio_polar_prior(alpha, alpha_min=0, alpha_max=35):
    '''smooth cut off of prior
    input:
    =====
    alpha = numpy array to be masked, in units of degrees
    alpha_min = float, if alpha is smaller than this value it's masked out
    alpha_max = float, if alpha is bigger than this value, it's masked out

    output:
    ======
    mask = numpy array that gives 1 if it is within the range and 0
            otherwise
    Stability: to be tested
    '''
    mask = np.logical_and(alpha < alpha_max, alpha > alpha_min)
    count = np.sum(mask)
    print "# of masked rows = {0}".format(alpha.size - count)
    print "# of remaining rows = {0}".format(count)
    print "% of original data remaining = {0:.2f}".format(
        count / alpha.size * 100.)
    return mask


def apply_radioprior(radiomask, dataarray):
    '''
    Checks if the length data array is the same as the prior mask
    if not, do nothing
    if lengths are the same, apply the prior
    starts examine if the length of mask is the same
    as the length of the array to be masked
    input:
    mask = numpy array with true or false as the values
    dataarray = numpy data array to be masked
    '''
    if len(radiomask) != len(dataarray):
        print 'length of mask and data array does not match!'
        print 'skipping the application of radio relic prior'
    else:
        # apply the mask
        temp = dataarray * radiomask

        counter = 0
        # removed the entries that are zero from the data array
        # to avoid the zero entries being binned
        for n in range(len(temp)):
            if temp[n] != 0.0:
                counter += 1
        # print 'number of non-zero entries after masking is', counter

        dataarray = numpy.zeros(counter)
        ncounter = 0
        for n in range(len(temp)):
            if temp[n] != 0.0:
                dataarray[ncounter] = temp[n]
                ncounter += 1
        # print 'number of non-zero entries for data array after masking is',
        # ncounter

    return dataarray


#-------helper functions-------------------------------------------------

def find_bin_ix(binedges, loc):
    """find the index in the numpy array binedges that corresponds to loc"""
    find_loc_i = binedges < loc
    return np.sum(find_loc_i)


def comb_zip(ls1, ls2):
    return [(lb1, lb2) for lb1 in ls1 for lb2 in ls2]


def round_to_n_sig_fig(x, n):
    if x > 1:
        return round(x, -int(np.log10(x)) + (n - 1))
    elif x < 1 and x > 0:
        return round(x, -int(np.log10(x)) + (n))


# --------modified version of  Will 's original functions----------------------
def histplot1d_pdf(x, prefix=None, prob=None, N_bins='knuth', histrange=None,
                   x_lim=None, y_lim=None, x_label=None, y_label=None,
                   legend=None, title=None, save=False, verbose=True,
                   plot=False):
    """plot data as normalized pdf
    plot the pdf of the histograms
    might want to return the pdf later on

    x = numpy array like object, can be dataframe columns
        data to be plotted on the x-axis
    prefix = string
        denotes the output file prefix
    prob = numpy array with same size as x
        denotes the weight to be put for correcting bias
    N_bins = integer
        denotes the number of bins
    histrange = a size 2 numpy array / list
        denotes the lower and upper range for making the histogram
    x_lim = a size 2 numpy array
    y_lim = a size 2 numpy array
    x_label = string
    y_label = string
    legend = string
    title = string
    """
    # compare bin width to knuth bin width
    if N_bins == 'knuth':
        binwidth, bins = de.knuth_bin_width(x, return_bins=True)
        knuth_N_bins = bins.size - 1
        N_bins = knuth_N_bins

    hist, binedges, tmp = \
        astroMLhist(x, bins=N_bins, histtype='step',
                    weights=prob, range=histrange, color='k', linewidth=2)

    # do not want to plot the graph without normalization
    # but the output of the binned array is needed for calculation
    # for location and confidence levels below
    pylab.close()

    #fig = pylab.figure()
    pylab.figure()
    # plot the normalized version (pdf) of the data
    # note that we use the astroML version of histogram
    # which infers the appropriate size of bin width
    pdf, binedges, holder = astroMLhist(x, bins=N_bins, histtype='step',
                                        weights=prob, range=histrange,
                                        color='k', linewidth=2,
                                        normed=True)
    pylab.close()
    # Calculate the location and %confidence intervals
    # Since my location and confidence calculations can't take weighted
    # data I need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned = \
                numpy.ones(hist[i]) * (binedges[i] + binedges[i + 1]) / 2
        else:
            x_temp = numpy.ones(hist[i]) * (binedges[i] + binedges[i + 1]) / 2
            x_binned = numpy.concatenate((x_binned, x_temp))

    loc = biweightLoc(x_binned)
    ll_68, ul_68 = bcpcl(loc, x_binned, 1)
    ll_95, ul_95 = bcpcl(loc, x_binned, 2)
    results = [loc, ll_68, ul_68, ll_95, ul_95]

    if verbose is True:
        print '{0}, {1:0.4f}, {2:0.4f},'.format(prefix, loc, ll_68) + \
            '{0:0.4f}, {1:0.4f}, {2:0.4f}'.format(ul_68, ll_95, ul_95)

    # Create location and confidence interval line plots
    # find the binedge that the location falls into
    # so that the line indicating the location only extends to top of
    # histogram
    loc_ix = find_bin_ix(binedges, loc)
    ll_68_ix = find_bin_ix(binedges, ll_68)
    ul_68_ix = find_bin_ix(binedges, ul_68)
    ll_95_ix = find_bin_ix(binedges, ll_95)
    ul_95_ix = find_bin_ix(binedges, ul_95)

    if save is False and plot is False:
        return results

    pylab.plot((loc, loc), (0, pdf[loc_ix - 1]), ls='--', lw=1, color="k")

    width = binedges[ll_68_ix + 1] - binedges[ll_68_ix]
    for i in range(ll_68_ix, ul_68_ix):
        pylab.bar(binedges[i], pdf[i], width, lw=0, color="b", alpha=.6)
    for i in range(ll_95_ix, ul_95_ix):
        pylab.bar(binedges[i], pdf[i], width, lw=0, color="b", alpha=.3)

    if x_label is not None:
        pylab.xlabel(x_label, fontsize=15)
    if y_label is not None:
        pylab.ylabel(y_label, fontsize=15)
    if x_lim is not None:
        pylab.xlim(x_lim)
    if y_lim is not None:
        pylab.ylim(y_lim)
    if legend is not None:
        pylab.legend()
    if title is not None:
        pylab.title(title, fontsize=16)

    # set font size
    # fontsize=14
    # ax = pylab.gca()
    # for tick in ax.xaxis.get_major_ticks():
        # tick.label1.set_fontsize(fontsize)

    if plot:
        pylab.show()

    if save is False:
        return results

    filename = prefix + '_histplot1D'
    pylab.savefig(filename + '.pdf', bbox_inches='tight')

    return results


def histplot1d(x, prefix=None, prob=None, norm=False, N_bins='knuth',
               histrange=None, x_lim=None, y_lim=None, x_label=None,
               y_label=None, legend=None, save=False, verbose=False,
               plot=False):
    """summarize the CI of the data
    plots data after binning and weighting the data appropriately
    """

    if N_bins == 'knuth':
        binwidth, bins = de.knuth_bin_width(x, return_bins=True)
        knuth_N_bins = bins.size - 1
        N_bins = knuth_N_bins
    # elif type(N_bins) is int:
    #    # compare bin width to knuth bin width
    #    print "specified bin width is {0}, Knuth bin size is {1}".format(
    #        N_bins, knuth_N_bins)

    #fig = pylab.figure()
    pylab.figure()
    hist, binedges, tmp = astroMLhist(
        x, bins=N_bins, histtype='step', weights=prob, range=histrange,
        color='k', linewidth=2, normed=norm)
    if plot is False:
        pylab.close()

    # Calculate the location and %confidence intervals
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations

    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned = numpy.ones(
                hist[i]) * (binedges[i] + binedges[i + 1]) / 2
        else:
            x_temp = numpy.ones(hist[i]) * (binedges[i] + binedges[i + 1]) / 2
            x_binned = numpy.concatenate((x_binned, x_temp))

    loc = biweightLoc(x_binned)
    ll_68, ul_68 = bcpcl(loc, x_binned, 1)
    ll_95, ul_95 = bcpcl(loc, x_binned, 2)
    results = [loc, ll_68, ul_68, ll_95, ul_95]

    if save is False and plot is False:
        return results

    # Create location and confidence interval line plots
    # find the binedge that the location falls into
    # so that the line indicating the location only extends to top of
    # histogram
    loc_ix = find_bin_ix(binedges, loc)
    ll_68_ix = find_bin_ix(binedges, ll_68)
    ul_68_ix = find_bin_ix(binedges, ul_68)
    ll_95_ix = find_bin_ix(binedges, ll_95)
    ul_95_ix = find_bin_ix(binedges, ul_95)

    pylab.plot((loc, loc), (0, hist[loc_ix - 1]), ls='--', lw=1, color="k")

    width = binedges[ll_68_ix + 1] - binedges[ll_68_ix]
    for i in range(ll_68_ix, ul_68_ix):
        pylab.bar(binedges[i], hist[i], width, lw=0, color="b", alpha=.6)
    for i in range(ll_95_ix, ul_95_ix):
        pylab.bar(binedges[i], hist[i], width, lw=0, color="b", alpha=.3)

    if x_label is not None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label is not None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim is not None:
        pylab.xlim(x_lim)
    if y_lim is not None:
        pylab.ylim(y_lim)
    if legend is not None:
        pylab.legend()
    # fontsize=14
    #ax = pylab.gca()
    # for tick in ax.xaxis.get_major_ticks():
        # tick.label1.set_fontsize(fontsize)

    if plot:
        pylab.show()

    if save is False:
        return results

    assert prefix is not None, "prefix cannot be None"

    filename = prefix + '_histplot1D'
    if save:
        print "saving file to {0}".format(filename)
        pylab.savefig(filename + '.png', dpi=300, bbox_inches='tight')

    if verbose:
        print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(
            prefix, loc, ll_68, ul_68, ll_95, ul_95)

    return results


def histplot1d_part(ax, x, prob=None, N_bins='knuth', histrange=None,
                    x_lim=None, y_lim=None):
    '''
    This take the additional value of an array axes. for use with subplots
    similar to histplot1d but for subplot purposes I believe
    '''
    # compare bin width to knuth bin width
    # if type(N_bins) is int:
    #    print "specified bin width is {0}, Knuth bin size is {1}".format(
    #        N_bins, knuth_N_bins)
    if N_bins == 'knuth':
        binwidth, bins = de.knuth_bin_width(x, return_bins=True)
        knuth_N_bins = bins.size - 1
        N_bins = knuth_N_bins

    hist, binedges, tmp = ax.hist(
        x, bins=N_bins, histtype='step', weights=prob, range=histrange,
        color='k', linewidth=1)

    # Calculate the location and %confidence intervals
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned = \
                numpy.ones(hist[i]) * (binedges[i] + binedges[i + 1]) / 2
        elif numpy.size(x_binned) == 0:
            x_binned = \
                numpy.ones(hist[i]) * (binedges[i] + binedges[i + 1]) / 2
        else:
            x_temp = \
                numpy.ones(hist[i]) * (binedges[i] + binedges[i + 1]) / 2
            x_binned = numpy.concatenate((x_binned, x_temp))
    loc = biweightLoc(x_binned)
    ll_68, ul_68 = bcpcl(loc, x_binned, 1)
    ll_95, ul_95 = bcpcl(loc, x_binned, 2)

    # Create location and confidence interval line plots
    # find the binedge that the location falls into
    # so that the line indicating the location only extends to top of
    # histogram
    loc_ix = find_bin_ix(binedges, loc)
    ll_68_ix = find_bin_ix(binedges, ll_68)
    ul_68_ix = find_bin_ix(binedges, ul_68)
    ll_95_ix = find_bin_ix(binedges, ll_95)
    ul_95_ix = find_bin_ix(binedges, ul_95)

    ax.plot((loc, loc), (0, hist[loc_ix - 1]), ls='--', lw=1, color="k")

    width = binedges[ll_68_ix + 1] - binedges[ll_68_ix]
    for i in range(ll_68_ix, ul_68_ix):
        ax.bar(binedges[i], hist[i], width, lw=0, color="b", alpha=.6)
    for i in range(ll_95_ix, ul_95_ix):
        ax.bar(binedges[i], hist[i], width, lw=0, color="b", alpha=.3)

    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)
    return loc, ll_68, ul_68, ll_95, ul_95


def histplot2d(x, y, prefix, prob=None, N_bins=100, histrange=None,
               x_lim=None, y_lim=None, x_label=None, y_label=None,
               legend=None, save=False):
    '''
    Input:
    plot 2d histogram of 2 data arrays
    x = [1D array of N floats]
    y = [1D array of N floats]
    prefix = [string] prefix of output file
    prob = [None] or [1D array of N floats] weights to apply to each (x,y) pair
    N_bins = [integer] the number of bins in the x and y directions
    histrange = [None] or [array of floats: (x_min,x_max,y_min,y_max)] the range
        over which to perform the 2D histogram and estimate the confidence
        intervals
    x_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    y_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    x_label = [None] or [string] the plot's x-axis label
    y_label = [None] or [string] the plot's y-axis label
    legend = [None] or [True] whether to display a legend or not
    '''
    # prevent masked array from choking up the 2d histogram function
    x = np.array(x)
    y = np.array(y)
    # Create the confidence interval plot
    if histrange is None:
        if prob is not None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, weights=prob)
        elif prob is None:
            H, xedges, yedges = numpy.histogram2d(x, y, bins=N_bins)
    else:
        if prob is not None:
            H, xedges, yedges = \
                numpy.histogram2d(x, y, bins=N_bins,
                                  range=[[histrange[0], histrange[1]],
                                         [histrange[2], histrange[3]]],
                                  weights=prob)
        elif prob is None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins,
                range=[[histrange[0], histrange[1]],
                       [histrange[2], histrange[3]]])
    H = numpy.transpose(H)
    # Flatten H
    h = numpy.reshape(H, (N_bins**2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum / h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum / h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]
    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x = xedges[:-1] + (xedges[1] - xedges[0]) / 2
    y = yedges[:-1] + (yedges[1] - yedges[0]) / 2
    X, Y = numpy.meshgrid(x, y)

    fig = pylab.figure()
    # Contours
    pylab.contour(X, Y, H, (h_2sigma, h_1sigma), linewidths=(2, 2))
    # imshow
    im = pylab.imshow(H, cmap=pylab.cm.gray_r, interpolation="gaussian")
    #pylab.pcolor(X, Y, H, cmap=pylab.cm.gray_r, interpolation="gaussian")

    if x_label is not None:
        pylab.xlabel(x_label, fontsize=14)
    if y_label is not None:
        pylab.ylabel(y_label, fontsize=14)
    if x_lim is not None:
        pylab.xlim(x_lim)
    if y_lim is not None:
        pylab.ylim(y_lim)
    if legend is not None:
        # Dummy lines for legend
        # 800000 is for Maroon color - 68% 1 sigma
        # 0000A0 is for blue color - 95% confidence 3 sigma
        pylab.plot((0, 1), (0, 1), c='#800000', linewidth=2, label=('68%'))
        pylab.plot((0, 1), (0, 1), c='#0000A0', linewidth=2, label=('95%'))
        pylab.legend(scatterpoints=1)
    fontsize = 15
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    if save:
        filename = prefix + '_histplot2d'
        pylab.savefig(filename, dpi=300, bbox_inches='tight')

    return fig


def histplot2d_part(ax, x, y, prob=None, N_bins=100, histrange=None,
                    x_lim=None, y_lim=None):
    '''
    similar to histplot2d
    This take the additional value of an array axes. for use with subplots
    Input:
    x = [1D array of N floats]
    y = [1D array of N floats]
    prefix = [string] prefix of output file
    prob = [None] or [1D array of N floats] weights to apply to each (x,y) pair
    N_bins = [integer] the number of bins in the x and y directions
    histrange = [None] or [array of floats: (x_min,x_max,y_min,y_max)] the range
        over which to perform the 2D histogram and estimate the confidence
        intervals
    x_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    y_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    x_label = [None] or [string] the plot's x-axis label
    y_label = [None] or [string] the plot's y-axis label
    legend = [None] or [True] whether to display a legend or not
    '''
    # prevent masked array from choking up the 2d histogram function
    x = np.array(x)
    y = np.array(y)

    # Create the confidence interval plot
    assert prob is not None, "there is no prob given for weighting"

    if histrange is None:
        if prob is not None:
            H, xedges, yedges = \
                numpy.histogram2d(x, y, bins=N_bins, weights=prob)
        elif prob is None:
            H, xedges, yedges = numpy.histogram2d(x, y, bins=N_bins)
    else:
        if prob is not None:
            H, xedges, yedges = \
                numpy.histogram2d(x, y, bins=N_bins,
                                  range=[[histrange[0], histrange[1]],
                                         [histrange[2], histrange[3]]],
                                  weights=prob)
        elif prob is None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]],
                                          [histrange[2], histrange[3]]])
    H = numpy.transpose(H)
    # Flatten H
    h = numpy.reshape(H, (N_bins ** 2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum / h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum / h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]

    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x = xedges[:-1] + (xedges[1] - xedges[0]) / 2
    y = yedges[:-1] + (yedges[1] - yedges[0]) / 2
    X, Y = numpy.meshgrid(x, y)

    # can use pcolor or imshow to show the shading instead
    ax.pcolormesh(X, Y, H, cmap=pylab.cm.gray_r, shading='gouraud')
    ax.contour(X, Y, H, (h_2sigma, h_1sigma), linewidths=(2, 2),
               colors=((158 / 255., 202 / 255., 225 / 255.),
                       (49 / 255., 130 / 255., 189 / 255.)))

    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)


def histplot2dTSC(x, y, prefix, prob=None, N_bins=100, histrange=None,
                  x_lim=None, y_lim=None, x_label=None, y_label=None,
                  legend=None):
    '''
    this is the one for generating the 2d plot in will 's paper?...
    Input:
    x = [1D array of N floats]
    y = [1D array of N floats]
    prefix = [string] prefix of output file
    prob = [None] or [1D array of N floats] weights to apply to each (x,y) pair
    N_bins = [integer] the number of bins in the x and y directions
    histrange = [None] or [array of floats: (x_min,x_max,y_min,y_max)] the range
        over which to perform the 2D histogram and estimate the confidence
        intervals
    x_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    y_lim = [None] or [array of floats: (x_min,x_max)] min and max of the range
        to plot
    x_label = [None] or [string] the plot's x-axis label
    y_label = [None] or [string] the plot's y-axis label
    legend = [None] or [True] whether to display a legend or not
    '''
    # Input calculated v and t parameters for other Dissociative Mergers
    v_bullet_analytic = 3400
    t_bullet_analytic = 0.218

    v_bullet_sf07 = 3400
    t_bullet_sf07 = 0.18

    v_macs = 2000
    t_macs = 0.255

    v_a520 = 2300
    t_a520 = 0.24

    v_pandora = 4045
    t_pandora = 0.162

    # Create the confidence interval plot
    if histrange is None:
        if prob is not None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, weights=prob)
        elif prob is None:
            H, xedges, yedges = numpy.histogram2d(x, y, bins=N_bins)
    else:
        if prob is not None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins,
                range=[[histrange[0], histrange[1]],
                       [histrange[2], histrange[3]]],
                weights=prob)
        elif prob is None:
            H, xedges, yedges = numpy.histogram2d(
                x, y, bins=N_bins, range=[[histrange[0], histrange[1]],
                                          [histrange[2], histrange[3]]])
    H = numpy.transpose(H)
    # Flatten H
    h = numpy.reshape(H, (N_bins ** 2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum / h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum / h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]
    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x = xedges[:-1] + (xedges[1] - xedges[0]) / 2
    y = yedges[:-1] + (yedges[1] - yedges[0]) / 2
    X, Y = numpy.meshgrid(x, y)

    fig = pylab.figure()
    # Countours
    CS = pylab.contour(X, Y, H, (h_2sigma, h_1sigma), linewidths=(2, 2))
    # imshow
    #im = pylab.imshow(H,cmap=pylab.cm.gray)
    pylab.pcolor(X, Y, H, cmap=pylab.cm.gray_r)

    # Data points for other dissociative mergers
    pylab.scatter(v_bullet_sf07, t_bullet_sf07, s=140,
                  c='k', marker='d', label="Bullet SF07")
    # pylab.scatter(v_macs,t_macs,s=140, c='0.4',markeredgecolor='0.4',
    #    marker='^',label='MACS J0025.4')
    # pylab.scatter(v_a520,t_a520,s=140,c='0.4',markeredgecolor='0.4',
    # marker='o',label='A520')
    # pylab.scatter(v_pandora,t_pandora,s=140,c='0.4',markeredgecolor='0.4',
    # marker='p',label='A2744')

    if x_label is not None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label is not None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim is not None:
        pylab.xlim(x_lim)
    if y_lim is not None:
        pylab.ylim(y_lim)
    if legend is not None:
        # Dummy lines for legend
        pylab.plot((0, 1), (0, 1), c='#800000', linewidth=2, label=('68%'))
        pylab.plot((0, 1), (0, 1), c='#0000A0', linewidth=2, label=('95%'))
        pylab.legend(scatterpoints=1)
    fontsize = 20
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    filename = prefix + '_histplot2dTSC.pdf'
    pylab.savefig(filename, dpi=300, bbox_inches='tight')

    return fig


def kde1d(x, plotObj, prob=None, kernel="gau", bw="scott", fft=True,
          gridsize=None, weights=None,
          adjust=1, cut=3, clip=(-np.inf, np.inf), xlabel=None, lim=None,
          labelsize=None, legendloc=None, **kwargs):
    """wrapper around the statsmodel.api.nonparametric.KDEUnivariate()
    and plot

    parameters
    ==========
    plotObj = matplotlib axis object, e.g. ax created by plt.subplots() or
        matplotlib.pyplot
    x = numpy array, the data that you try to visualize
    prob = numpy array, has the same as x
    kernel = string, what kernel to use, see KDEUnivariate() documentation
        for options
    bw = string or integer, denotes the binwidth
    fft = logical, if fast fourier transform should be used while
        computing the kde
    see KDEUnivaraiate documentation for the rest of the parameters

    label = string, denotes what would be put as xlabel of the plot
    lim = float tuple of length 2, denotes the xlim on the plot

    returns
    =======
    support = numpy array
        the corresponding x value of each element of the returned pdf
    pdf = numpy array, the corresponding pdf
    """
    kde = nonparametric.KDEUnivariate(x)
    kde.fit(kernel, bw, fft, gridsize=gridsize, cut=cut, clip=clip,
            weights=weights)
    support, density = kde.support, kde.density

    plotObj.plot(support, density, **kwargs)
    if isinstance(plotObj, ModuleType):
        plotObj.ylabel("PDF")
        if xlabel is not None:
            plotObj.xlabel(xlabel, fontsize=labelsize)
        if lim is not None:
            plotObj.xlim(lim)
        if legendloc is not None:
            plotObj.legend(loc=legendloc)
    else:
        if xlabel is not None:
            plotObj.set_xlabel(xlabel)
        if lim is not None:
            plotObj.set_xlim(lim)

    return support, density


def central_CI(support, density, level=68, lim=None):
    """returns the central credible intervals for a kde estimate from kde1d
    of a posterior
    parameters
    ==========
    support = numpy array
    density = numpy array that has the same length as the support
    level = float, indicates what percentile to include between
    lim = tuple of float of length 2

    returns
    ======
    low_ix = index of the lower limit, support[low_ix] gives the estimate
    up_ix = index of the upper limit, support[up_ix] gives the estimate

    warning: the index may be off by one depending on how you plot
    but kde is smooth enough it should not matter

    stability:
    ========
    work in progress
    """

    if lim is not None:
        lix = 0
        while(support[lix] < lim[0]):
            lix += 1
    else:
        lix = 0

    sig = (1 - level / 100.) / 2.
    total = np.sum(density[lix:])
    exclude_reg = total * sig

    low_exc = 0
    low_ix = 0
    while(low_exc < exclude_reg):
        low_exc += density[low_ix]
        low_ix += 1

    up_exc = 0
    up_ix = density.size - 1
    while(up_exc < exclude_reg):
        up_exc += density[up_ix]
        up_ix -= 1
    return low_ix, up_ix


def CI_loc_plot(x, plotObj, c='b', prob=None, kernel="gau", bw="scott",
                fft=True, gridsize=None, level=68, adjust=1, cut=3,
                clip=(-np.inf, np.inf), xlabel=None, lim=None,
                weights=None,
                labelsize=None, legendloc=None, **kwargs):
    support, den = kde1d(x, plotObj, prob=prob, kernel=kernel, bw=bw,
                         fft=fft, gridsize=gridsize, adjust=adjust, cut=cut,
                         clip=clip, xlabel=xlabel, lim=lim,
                         weights=weights,
                         labelsize=labelsize, legendloc=legendloc,
                         **kwargs)
    low68_ix, up68_ix = central_CI(support, den, level=68, lim=lim)
    low95_ix, up95_ix = central_CI(support, den, level=95, lim=lim)
    plt.fill_between(support[low68_ix: up68_ix],
                     den[low68_ix: up68_ix], alpha=0.5, color=c)
    plt.fill_between(support[low95_ix: up95_ix],
                     den[low95_ix: up95_ix], alpha=0.2, color=c)
    loc_ix = low68_ix
    loc = C_BI(x)
    while(support[loc_ix] < loc):
        loc_ix += 1
    ylim = plt.ylim()
    plt.axvline(loc, ymin=0,
                ymax=den[loc_ix] / ylim[1], ls='--', lw=3, c='k')

    return


def percentdiff(x, prefix, prob=None, N_bins=100, histrange=None, x_lim=None,
                y_lim=None, x_label=None, y_label=None, legend=None):
    """
    this function takes in parameter arrays
    bins the data then calculates the percentage difference
    actually forgot where I used this
    """
    #fig = pylab.figure()
    # find out size of array
    totalsize = len(x)
    print 'data size of each variable is ', totalsize
    # divide total number of data points into nparts-1
    nparts = 101
    d = (nparts - 1, 5)
    reduced_d = (nparts - 2, 5)
    data = numpy.zeros(d)
    x_perdiff = numpy.zeros(reduced_d)

    # iterate from 1 to nparts-1
    for n in range(1, nparts):
            # print i,"th iteration"
        # size of each of the n parts:
        partsize = totalsize * (n) / (nparts - 1)
        hist, binedges, tmp = \
            astroMLhist(x[:partsize], bins="knuth",
                        histtype='step', weights=prob[:partsize],
                        range=histrange, color='k', linewidth=2)

        # Calculate the location and confidence intervals
        # Since my location and confidence calculations can't take weighted data I
        # need to use the weighted histogram data in the calculations
        for i in numpy.arange(N_bins):
            if i == 0:
                x_binned = \
                    numpy.ones(hist[i]) * (binedges[i] + binedges[i + 1]) / 2
            else:
                x_temp = \
                    numpy.ones(hist[i]) * (binedges[i] + binedges[i + 1]) / 2
                x_binned = numpy.concatenate((x_binned, x_temp))
        # print 'len of x_binned is ',len(x_binned)
        loc = biweightLoc(x_binned)
        ll_68, ul_68 = bcpcl(loc, x_binned, 1)
        ll_95, ul_95 = bcpcl(loc, x_binned, 2)
        # this will store the data
        data[n - 1] = [loc, ll_68, ul_68, ll_95, ul_95]

    filename = prefix + '_histplot1D_percentdiff.png'

    pylab.axvline(loc, ls='--', linewidth=2, label='$C_{BI}$', color="k")
    pylab.axvline(ll_68, ls='-.', linewidth=2, color='#800000',
                  label='68% $IC_{B_{BI}}$')
    pylab.axvline(ul_68, ls='-.', linewidth=2, color='#800000')
    pylab.axvline(ll_95, ls=':', linewidth=2, color='#0000A0',
                  label='95% $IC_{B_{BI}}$')
    pylab.axvline(ul_95, ls=':', linewidth=2, color='#0000A0')

    pylab.savefig(filename, dpi=300, bbox_inches='tight')
    pylab.close()

    print '\n' + prefix + ' data is '
    print data
    print '      '

    for n in range(1, nparts - 1):
        x_perdiff[n - 1] = (data[nparts - 2] - data[n - 1]) * 2 * \
            100 / (data[nparts - 2] + data[n - 1])
    print prefix + ' per diff is '
    print x_perdiff
    print '      '

    # this will invert the array, now disabled
    #x_perdiff = x_perdiff[::-1]

    return x_perdiff


def boxplot(result, y, tick_label, key, size=16, save=False):
    """ plot my version of boxplot
    usage: apply to each row of a dataframe with the results from
    histplot1d

    Parameters
    ------------
    result = list of numbers of size 5
        denoting [location, CI68_lower, CI68_upper, CI95_lower, CI95_upper]
    y = vertical location
    tick_label = list of numbers
        to be used as y_ticks
    key = string
        denotes the name of the variable on the x axis
    size = integer, font size of labels

    example usage:
    -------------
    tick_label = np.arange(alpha_lower - alpha_interval,
                    alpha_upper + alpha_interval, alpha_interval)
    for key in par:
        for alpha in alpha0:
            # results is a dataframe with index = [[str(alpha)], [key]]
            boxplot(results.ix[str(alpha), key], alpha, tick_label, key)
            plt.savefig("CI_" + key + '.pdf', bbox_inches='tight')
        plt.close()

    """
    import matplotlib.pyplot as plt
    CI68 = np.arange(result[1], result[2], 0.01)
    CI95 = np.arange(result[3], result[4], 0.01)
    plt.plot(CI95, y * np.ones(CI95.size), 'g-', label='95% CI')
    plt.plot(CI68, y * np.ones(CI68.size), 'b-')
    plt.plot(CI68[0], y, 'b|', markeredgewidth=2, label='68% CI')
    plt.plot(CI68[-1], y, 'b|', markeredgewidth=2)
    plt.plot(result[0], y, 'r|', markeredgewidth=2, label='location')
    plt.yticks(tick_label)
    plt.ylabel(r'$\alpha_0$', size=size)
    plt.xlabel(key, size=size)
    plt.show()

    return


def N_by_M_plot_contour(data, Nvar_list, Mvar_list, space, axlims=None,
                        Nbins_2D=None, axlabels=None, xlabel_to_rot=None,
                        histran=None, figsize=6, fontsize=11, save=False,
                        path="./", prefix=None, suffix=".png"):
    """create a N by M matrix of 2D contour plots
    data = dataframe that contain the data of all the variables to be plots
    Nvar_list = list of strings - denotes the column header names
        that needs to be plotted on x-axes, col names correspond to xlabel
    Mvar_list = list of strings - denotes the column header names
        that needs to be plotted on y-axes, col names correspond to ylabel
    space = float, px of space that is added between subplots
    axlims = dictionary, keys are the strings in var_list,
        each value is a tuple of (low_lim, up_lim) to denote the limit
        of values to be plotted
    Nbins_2D = dictionary, keys are in format of tuples of
        (x_col_str, y_col_str) to denote which subplot you are referring to
    axlabels = dictionary, keys correspond to the variable names, values
        strings that correspond to the labels to be put on corresponding
        axes
    xlabel_to_rot = dictionary,
        key is the the key for the labels to be rotated,
        value is the degree to be rotated
    histran = dictionary,
        some keys has to be the ones for the plots, value are in
        form of (lowerhist_range, upperhist_range)
    figsize = integer, figuares are squared this refers to the side length
    fontsize = integer, denotes font size of the labels
    save = logical, denotes if plot should be saved or not
    prefix = string, prefix of the output plot file
    path = string, path of the output plot file
    suffix = string, file extension of the output plot file

    Stability: Not entirely tested, use at own risk

    Author: Karen Ng
    """
    from matplotlib.ticker import MaxNLocator

    N = len(Nvar_list)
    M = len(Mvar_list)

    print 'creating input-input figure with dimension ' + \
        '{0} rows by {1} cols'.format(M, N)

    # begin checking if inputs make sense
    assert N <= len(axlabels), "length of axlabels is wrong"
    assert M <= len(axlabels), "length of axlabels is wrong"

    compare_Nvar = np.sum([Nvar in data.columns for Nvar in Nvar_list])
    assert compare_Nvar == len(Nvar_list), "variable to be plotted not in df"

    compare_Mvar = np.sum([Mvar in data.columns for Mvar in Mvar_list])
    assert compare_Mvar == len(Mvar_list), "variable to be plotted not in df"

    keys = comb_zip(Nvar_list, Mvar_list)
    if Nbins_2D is not None:
        compare_keys = np.sum([key in Nbins_2D.keys() for key in keys])
        assert compare_keys == len(keys), "Nbins_2D key error"
    else:
        Nbins_2D = {key: 50 for key in keys}

    if axlims is None:
        axlims = {key: (None, None) for key in Nvar_list + Mvar_list}

    if save:
        assert prefix is not None, "prefix for output file cannot be none"

    # impossible for the matrix plot not to be squared in terms of dimensions
    # set each of the subplot to be squared with the figsize option
    f, axarr = plt.subplots(M, N, figsize=(figsize * N / M, figsize))
    f.subplots_adjust(wspace=space, hspace=space)

    # remove unwanted row axes tick labels
    plt.setp([a.get_xticklabels() for i in range(M - 1)
              for a in axarr[i, :]], visible=False)

    # remove unwanted column axes tick labels
    plt.setp([a.get_yticklabels() for i in range(1, N)
              for a in axarr[:, i]], visible=False)

    # rotate the xlabels appropriately
    if xlabel_to_rot is not None:
        match_ix = [Nvar_list.index(item) for item in Nvar_list]
        # ok to use for-loops for small number of iterations
        for ix in match_ix:
            labels = axarr[M - 1, ix].get_xticklabels()
            for label in labels:
                label.set_rotation(xlabel_to_rot[Nvar_list[ix]])

    # create axes labels
    if axlabels is not None:
        for j in range(M):
            axarr[j, 0].set_ylabel(axlabels[Mvar_list[j]], fontsize=fontsize)
        for i in range(N):
            axarr[M - 1, i].set_xlabel(axlabels[Nvar_list[i]],
                                       fontsize=fontsize)

    # fix the number of bins on the axis by nbins
    # avoid overlapping lowest and highest ticks mark with prune
    # option pruning both upper and lower bound
    for n in range(N):
        for m in range(M):
            ax2 = axarr[m, n]
            ax2.xaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))
            ax2.yaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))

    # start plotting the 2D contours
    for i in range(M):
        for j in range(N):
            # print "axarr[i, j] has indices {0}".format((i, j))
            # print "x axis label = {0}".format(Nvar_list[j])
            # print "y axis label = {0}".format(Mvar_list[i])
            histplot2d_part(axarr[i, j],
                            data[Nvar_list[j]],
                            data[Mvar_list[i]],
                            prob=data['prob'],
                            N_bins=Nbins_2D[(Nvar_list[j],
                                             Mvar_list[i])],
                            x_lim=axlims[Nvar_list[j]],
                            y_lim=axlims[Mvar_list[i]])

    if save:
        print "saving plot to {0}".format(path + prefix + suffix)
        plt.savefig(path + prefix + suffix, dpi=200, bbox_inches='tight')

    return


def N_by_N_lower_triangle_plot(data, space, var_list, axlims=None,
                               Nbins_2D=None, axlabels=None, N_bins=None,
                               xlabel_to_rot=None, histran=None, figsize=6,
                               fontsize=12, save=False, prefix=None,
                               suffix=".png", path="./"):
    """ create a N by N matrix of plots
    with the top plot of each row showing a density plot in 1D
    and the remaining plots being 2D contour plots
    df = dataframe that contain the data of all the variables to be plots
    space = float, px of space that is added between subplots
    var_list = list of strings - denotes the column header names
        that needs to be plotted
    axlims = dictionary, keys are the strings in var_list,
        each value is a tuple of (low_lim, up_lim) to denote the limit
        of values to be plotted
    Nbins_2D = dictionary, keys are in format of tuples of
        (x_col_str, y_col_str) to denote which subplot you are referring to
    axlabels = dictionary, keys correspond to the variable names
    xlabel_to_rot = dictionary,
        key is the the key for the labels to be rotated,
        value is the degree to be rotated
    histran = dictionary,
        some keys has to be the ones for the plots, value are in
        form of (lowerhist_range, upperhist_range)
    figsize = integer, figuares are squared this refers to the side length
    fontsize = integer, denotes font size of the labels
    save = logical, denotes if plot should be saved or not
    prefix = string, prefix of the output plot file
    path = string, path of the output plot file
    suffix = string, file extension of the output plot file

    Stability: Not entirely tested, use at own risk
    """
    from matplotlib.ticker import MaxNLocator

    # begin checking if inputs make sense
    N = len(var_list)
    assert N <= len(axlabels), "length of axlabels is wrong"
    assert N >= 2, "lower triangular contour plots require more than 2\
        variables in the data"

    for var in var_list:
        assert var in data.columns, "variable to be plotted not in df"

    if axlabels is None:
        axlabels = {key: key for key in var_list}

    if xlabel_to_rot is None:
        xlabel_to_rot = {key: 0 for key in var_list}

    if histran is None:
        histran = {key: None for key in var_list}

    if axlims is None:
        axlims = {key: (None, None) for key in var_list}

    if Nbins_2D is None:
        keys = comb_zip(var_list, var_list)
        Nbins_2D = {key: 50 for key in keys}

    if N_bins is None:
        N_bins = {key: 'knuth' for key in var_list}

    if save:
        assert prefix is not None, "prefix for output file cannot be none"

    # impossible for the matrix plot not to be squared in terms of dimensions
    # set each of the subplot to be squared with the figsize option
    f, axarr = pylab.subplots(N, N, figsize=(figsize, figsize))
    f.subplots_adjust(wspace=space, hspace=space)

    # remove unwanted plots on the upper right
    plt.setp([a.get_axes() for i in range(N - 1)
              for a in axarr[i, i + 1:]], visible=False)

    # remove unwanted row axes tick labels
    plt.setp([a.get_xticklabels() for i in range(N - 1)
              for a in axarr[i, :]], visible=False)

    # remove unwanted column axes tick labels
    plt.setp([axarr[0, 0].get_yticklabels()], visible=False)
    plt.setp([a.get_yticklabels() for i in range(N - 1)
              for a in axarr[i + 1, 1:]], visible=False)

    # create axes labels
    if axlabels is not None:
        for j in range(1, N):
            axarr[j, 0].set_ylabel(axlabels[var_list[j]], fontsize=fontsize)
        for i in range(N):
            axarr[N - 1, i].set_xlabel(axlabels[var_list[i]],
                                       fontsize=fontsize)

    for n in range(N):
        # avoid overlapping lowest and highest ticks mark
        # print "setting x and y tick freq for {0}".format((n, n))
        ax2 = axarr[n, n]
        ax2.xaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))

    # print "setting x and y tick freq for {0}".format((i, j))
    for i in range(N):
        for j in range(N):  # range(i)
            ax2 = axarr[i, j]
            ax2.yaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))
            ax2.xaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))

    # rotate the xlabels appropriately
    if xlabel_to_rot is not None:
        match_ix = [var_list.index(item) for item in var_list]
        # ok to use for-loops for small number of iterations
        for ix in match_ix:
            labels = axarr[N - 1, ix].get_xticklabels()
            for label in labels:
                label.set_rotation(xlabel_to_rot[var_list[ix]])

    # start plotting the diagonal
    for i in range(N):
        print "N_bins = {0}".format(N_bins[var_list[i]])
        histplot1d_part(axarr[i, i], np.array(data[var_list[i]]),
                        np.array(data['prob']),
                        N_bins=N_bins[var_list[i]],
                        histrange=histran[var_list[i]],
                        x_lim=axlims[var_list[i]])

    # start plotting the lower triangle when row no > col no
    for i in range(N):
        for j in range(i):
            histplot2d_part(axarr[i, j], np.array(data[var_list[j]]),
                            np.array(data[var_list[i]]),
                            prob=np.array(data['prob']),
                            N_bins=Nbins_2D[(var_list[j], var_list[i])],
                            x_lim=axlims[var_list[j]],
                            y_lim=axlims[var_list[i]])

    if save:
        print "saving plot to {0}".format(path + prefix + suffix)
        plt.savefig(path + prefix + suffix, dpi=200, bbox_inches='tight')

    return


def N_by_M_lower_tri_plot(Ndata, Mdata, space, Nvar_list, Mvar_list,
                          axlims=None, Nbins_2D=None, axlabels=None,
                          xlabel_to_rot=None, histran=None, figsize=6,
                          fontsize=12, save=False, prefix=None, path="./",
                          suffix=".png"):
    """ create a N by N lower triangular matrix of 2D contour plots
    This is a sensitivity analysis of how the dynamical simulation affects
    the output parameters

    Parameters
    ==========
    Ndata = dataframe that contain the data of all the variables to be plots
    Mdata = dictionary that contain data to be plotted
    space = float, px of space that is added between subplots
    var_list = list of strings - denotes the column header names
        that needs to be plotted
    axlims = dictionary, keys are the strings in var_list,
        each value is a tuple of (low_lim, up_lim) to denote the limit
        of values to be plotted
    Nbins_2D = dictionary, keys are in format of tuples of
        (x_col_str, y_col_str) to denote which subplot you are referring to
    axlabels = dictionary, keys correspond to the variable names
    xlabel_to_rot = dictionary,
        key is the the key for the labels to be rotated,
        value is the degree to be rotated
    histran = dictionary,
        some keys has to be the ones for the plots, value are in
        form of (lowerhist_range, upperhist_range)
    figsize = integer, figuares are squared this refers to the side length
    fontsize = integer, denotes font size of the labels
    save = logical, denotes if plot should be saved or not
    prefix = string, prefix of the output plot file
    path = string, path of the output plot file
    suffix = string, file extension of the output plot file

    Stability: Not entirely tested, use at own risk
    """
    from matplotlib.ticker import MaxNLocator

    # begin checking if inputs make sense
    N = len(Nvar_list)
    M = len(Mvar_list)
    if axlabels is not None:
        assert N <= len(axlabels), "length of axlabels is wrong"
        assert M <= len(axlabels), "length of axlabels is wrong"
    assert N >= 2, "lower triangular contour plots require more than 2\
        variables in the data"
    assert M == N, "lower triangular contour plots require more than \
        x dimensions == y dimensions"

    for var in Nvar_list:
        assert var in Ndata.columns, "N variable to be plotted not in df"
    for var in Mvar_list:
        assert var in Mdata.columns, "M variable to be plotted not in df"

    if axlabels is None:
        axlabels = {key: key for key in Nvar_list + Mvar_list}

    if xlabel_to_rot is None:
        xlabel_to_rot = {key: 0 for key in Nvar_list}

    if histran is None:
        histran = {key: None for key in Nvar_list}

    if axlims is None:
        axlims = {key: (None, None) for key in Nvar_list + Mvar_list}

    if Nbins_2D is None:
        keys = comb_zip(Nvar_list, Mvar_list)
        Nbins_2D = {key: 50 for key in keys}

    if save:
        assert prefix is not None, "prefix for output file cannot be none"

    # impossible for the matrix plot not to be squared in terms of dimensions
    # set each of the subplot to be squared with the figsize option
    f, axarr = pylab.subplots(N, N, figsize=(figsize, figsize))
    f.subplots_adjust(wspace=space, hspace=space)

    # remove unwanted plots on the upper right
    plt.setp([a.get_axes() for i in range(N - 1)
              for a in axarr[i, i + 1:]], visible=False)

    # remove unwanted row axes tick labels
    plt.setp([a.get_xticklabels() for i in range(N - 1)
              for a in axarr[i, :]], visible=False)

    # remove unwanted column axes tick labels
    #plt.setp([axarr[0, 0].get_yticklabels()], visible=False)
    plt.setp([a.get_yticklabels() for i in range(N - 1)
              for a in axarr[i + 1, 1:]], visible=False)

    # create axes labels
    # if axlabels is not None:
    for j in range(M):
        axarr[j, 0].set_ylabel(axlabels[Mvar_list[j]], fontsize=fontsize)
    for i in range(N):
        axarr[N - 1, i].set_xlabel(axlabels[Nvar_list[i]],
                                   fontsize=fontsize)

    # avoid overlapping lowest and highest ticks mark
    # for n in range(N):
    #    ## avoid overlapping lowest and highest ticks mark
    #    #ax1 = axarr[n, 0]
    #    #ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))

    #    #print "setting x tick freq for {0}".format((n, n))
    #    ax2 = axarr[n, n]
    #    ax2.xaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
    #    ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))

    for i in range(N):
        for j in range(N):
            # print "setting x tick freq for {0}".format((i, j))
            ax2 = axarr[i, j]
            ax2.xaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
            ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))

    # rotate the xlabels appropriately
    if xlabel_to_rot is not None:
        match_ix = [Nvar_list.index(item) for item in xlabel_to_rot.keys()]
        # ok to use for-loops for small number of iterations
        for ix in match_ix:
            labels = axarr[N - 1, ix].get_xticklabels()
            for label in labels:
                label.set_rotation(xlabel_to_rot[Nvar_list[ix]])

    # start plotting the diagonal
    # have to handle how prob is applied for Mdata !!!!!!
    # note that the xlim and ylims are supposed to be the same
    for i in range(N):
        histplot2d_part(axarr[i, i],
                        Ndata[Nvar_list[i]],
                        Mdata[Mvar_list[i]],
                        prob=Ndata['prob'],
                        N_bins=Nbins_2D[(Nvar_list[i], Mvar_list[i])],
                        x_lim=axlims[Mvar_list[i]],
                        y_lim=axlims[Mvar_list[i]])

    # start plotting the lower triangle when row no > col no
    # have to handle how prob is applied for Mdata !!!!!!
    for i in range(N):
        for j in range(i):
            histplot2d_part(axarr[i, j],
                            Ndata[Nvar_list[j]],
                            Mdata[Mvar_list[i]],
                            prob=Ndata['prob'],
                            N_bins=Nbins_2D[(Nvar_list[j],
                                             Mvar_list[i])],
                            x_lim=axlims[Mvar_list[j]],
                            y_lim=axlims[Mvar_list[i]])

    if save:
        print "saving plot to {0}".format(path + prefix + suffix)
        plt.savefig(path + prefix + suffix, dpi=200, bbox_inches='tight')

    return


def summarize_CI_as_df(data, par, N_bins="knuth", columns=None,
                       histran=None, verbose=False):
    """summarize CI as a dataframe / other suitable output formats

    Parameters
    ==========
    data = pandas dataframe with corresponding data variables, 'prob' has
        to be one of the variable in the dataframe
    par = list of strings, corresponds to the names of the variables to be
        summarized
    N_bins = integers, number of bins to be used for computing histograms
        by default use "knuth" rule, or accept a dictionary that has keys
        to be par, the values should be the bin widths
    columns = list of 5 strings denoting the location, 68 and 95 upper and
        lower credible intervals
        ["loc", "ci68low", "ci68up", "ci95low", "ci95up"]

    Returns
    ======
    dataframe

    Stability: works - if the res are nonsense for some variables,
    try adjusting the histran / leaving it as None
    """

    assert "prob" in data.columns, "prob is not one of the variables " + \
        " in the dataframe"

    if columns is None:
        columns = ["loc", "CI68low", "CI68up", "CI95low", "CI95up"]

    if histran is None:
        histran = dict.fromkeys(data.columns)

    # sometimes Knuth bins are problematic
    if isinstance(N_bins, str) or isinstance(N_bins, int):
        N_bins = {key: N_bins for key in par}

    # use list comprehension to loop through the results
    res = [
        histplot1d(
            np.array(
                data[key]),
            prob=np.array(
                data['prob']),
            N_bins=N_bins[key],
            histrange=histran[key],
            save=False,
            verbose=verbose,
            plot=False) for key in par]
    res = pd.DataFrame(res, columns=columns, index=par)

    return res


def save_CI_df_to_h5(res, par, prefix=None, mode="w", complib="blosc",
                     path="./", verbose=True):
    """
    res = pandas dataframe to be saved to file should contain
    par = list of strings, corresponds to the names of the variables to be
        summarized
    path = string, denotes path to for outputting the file
    verbose = logical, denotes if messages should be displayed or not
    mode = string, "w" or "a", for write only or append mode
    complib = "string", read df.to_hdf for options, default = "blosc"
    """
    assert prefix is not None, "output file prefix cannot be None"
    if verbose:
        "writing file to {0}".format(path + prefix + ".h5")
    res.to_hdf(path + prefix + ".h5", key="df", mode=mode, complib="blosc")

    return


def texline(res, m_res, columns, key, dp, sf):
    """more formatted table for the El Gordo paper
    res = panda df results with appropriate row labels
    m_res = panda df masked df results with appropriate row labels
    textfile = string path to the output texfile
    columns = column names for the res and m_res df
    key = string - row name for the key in the dataframe
    dp = string - format option for the precision
    sf = integer - significant figure
    """
    line = "{0:.{5}}&{1:.{5}}-{2:.{5}}&{3:.{5}}-{4:.{5}}&&".format(
        round_to_n_sig_fig(res.ix[key][columns[0]], sf),
        round_to_n_sig_fig(res.ix[key][columns[1]], sf),
        round_to_n_sig_fig(res.ix[key][columns[2]], sf),
        round_to_n_sig_fig(res.ix[key][columns[3]], sf),
        round_to_n_sig_fig(res.ix[key][columns[4]], sf), dp) + \
        "{0:.{5}}&{1:.{5}}-{2:.{5}}&{3:.{5}}-{4:.{5}}\\\\".format(
            round_to_n_sig_fig(m_res.ix[key][columns[0]], sf),
            round_to_n_sig_fig(m_res.ix[key][columns[1]], sf),
            round_to_n_sig_fig(m_res.ix[key][columns[2]], sf),
            round_to_n_sig_fig(m_res.ix[key][columns[3]], sf),
            round_to_n_sig_fig(m_res.ix[key][columns[4]], sf), dp) + "\n"
    return line


def save_CI_df_to_tex_table(res, m_res, texfile, columns=None, verbose=True):
    """more formatted table for the El Gordo paper
    res = panda df results with appropriate row labels
    m_res = panda df masked df results with appropriate row labels
    textfile = string path to the output texfile
    columns = column names for the res and m_res df
    """
    if columns is None:
        columns = ["loc", "CI68low", "CI68up", "CI95low", "CI95up"]
    F = open(texfile, "w+")
    if verbose:
        print "writing file to ", texfile
    F.write("\\begin{table*} \n")
    F.write("\\begin{minipage}{170mm} \n")
    F.write("\caption{Table of the output PDF properties of the model " +
            "variables and output variables from Monte Carlo simulation\n")
    F.write("\label{tab:outputs}}\n")
    F.write("\\begin{tabularx}{\\textwidth}" +
            "{@{\extracolsep{\\fill}}lccccccccc@{}}\n")
    F.write("\hline\n")
    F.write("\hline\n")
    F.write(
        "&&&&Default priors & & & & Default + polarization priors  \\\\ \n")
    F.write("\cmidrule{4-6} \cmidrule{8-10} \n")
    F.write("Variables & Units && Location & 68$\%$ CI $^{\dagger}$ &" +
            "95$\%$ CI && Location & 68$\%$ CI  & 95$\%$ CI \\\\ \n")
    F.write("\hline \n")
    key = "alpha"
    sf = 2
    dp = "0f"
    F.write("$\\alpha$ &(degree)&&" +
            texline(res, m_res, columns, key, dp, sf))
    dp = "2g"
    key = "d_proj"
    F.write("$d_{proj}$ &Mpc&&" +
            texline(res, m_res, columns, key, dp, sf))
    key = "d_max"
    F.write("$d_{max}$ &Mpc&&" +
            texline(res, m_res, columns, key, dp, sf))
    key = "d_3d"
    F.write("$d_{3D}$ &Mpc&&" +
            texline(res, m_res, columns, key, dp, sf))
    dp = "2g"
    key = "TSM_0"
    F.write("$TSC_0$&Gyr&&" +
            texline(res, m_res, columns, key, dp, sf))
    key = "TSM_1"
    F.write("$TSC_1$&Gyr&&" +
            texline(res, m_res, columns, key, dp, sf))
    key = "T"
    F.write("$T$&Gyr&&" +
            texline(res, m_res, columns, key, dp, sf))
    key = "v_3d_obs"
    dp = "0f"
    F.write("$v_{3D}(t_{obs})$ & \kilo \meter~\second$^{-1}$ &&" +
            texline(res, m_res, columns, key, dp, sf))
    key = "v_rad_obs"
    dp = "0f"
    F.write("$v_{rad}(t_{obs})$ & \kilo \meter~\second$^{-1}$ &&" +
            texline(res, m_res, columns, key, dp, sf))
    key = "v_3d_col"
    F.write("$v_{3D}(t_{col})$ & \kilo \meter~\second$^{-1}$ &&" +
            texline(res, m_res, columns, key, dp, sf))
    F.write("\\bottomrule" + "\n")
    F.write("\end{tabularx}\\\\" + "\n")
    F.write("\\footnotesize{$\dagger$ CI stands for credible interval}\\\\" +
            "\n")
    F.write("\end{minipage}" + "\n")
    F.write("\end{table*}" + "\n")
    F.close()
    return


def save_CI_df_to_tex(res, par, columns=None, prefix=None,
                      labels=None, units=None, verbose=True,
                      svar=None, path="./"):
    """
    res = pandas dataframe to be saved to file should contain
    par = dictionary, keys correspond to the names of the variables to be
        summarized, values correspond to the sign fig, 0 < value < 1
    columns = list of 5 strings denoting the location, 68 and 95 upper and
        lower credible intervals
        ["loc", "ci68low", "ci68up", "ci95low", "ci95up"]
    labels = dictionaries, keys are the same as data columns, values are
        strings that correspond to the labels
    units = dictionaries, keys are the same as data columns, values are
        strings that correspond to the units
    verbose = logical, denotes if messages should be displayed or not
    svar = dictionary, keys denote list of variables that would require
        different display of sign fig.,  values are integers to denote the
        number of sign. fig and these values > 1
    path = string, output path
    """
    if columns is None:
        columns = ["loc", "CI68low", "CI68up", "CI95low", "CI95up"]

    assert prefix is not None, "output file prefix cannot be None"
    assert labels is not None, "labels cannot be missing"
    assert units is not None, "units cannot be missing"

    if svar is not None:
        svar_check = np.sum([key in res.index for key in svar.keys()])
        assert svar_check == len(svar), "some svar keys are not in df"

    # prepare all the strings first
    # {0} is label name {1} is unit
    # {2} is loc {3} is CI68low {4} is CI68up {5} is CI95low
    # {6} is CI95up
    # if value < 1, then have to use {0:.{s}f}
    # or else value > 1 then have to use {0:.{s}g} where s is the integer
    # indicating the sign fig
    STR = ["${0}$&{1}&{2:.{s}f}&{3:.{s}f}-".format(
           labels[key], units[key],
           round_to_n_sig_fig(res.ix[key][columns[0]], par[key]),
           round_to_n_sig_fig(res.ix[key][columns[1]], par[key]),
           s=par[key]) +
           "{0:.{s}f}&{1:.{s}f}-{2:.{s}f}".format(
               round_to_n_sig_fig(res.ix[key][columns[2]], par[key]),
               round_to_n_sig_fig(res.ix[key][columns[3]], par[key]),
               round_to_n_sig_fig(res.ix[key][columns[4]], par[key]),
               s=par[key]) +
           "\\\\ \n" if res.ix[key][columns[0]] < 1 else
           "${0}$&{1}&{2:.{s}g}&{3:.{s}g}-".format(
               labels[key], units[key], res.ix[key][columns[0]],
               res.ix[key][columns[1]], s=par[key]) +
           "{0:.{s}g}&{1:.{s}g}-{2:.{s}g}".format(
               res.ix[key][columns[2]], res.ix[key][columns[3]],
               res.ix[key][columns[4]], s=par[key]) +
           "\\\\ \n" for key in par.keys()]
    if svar is not None:
        STR1 = ["${0}$&{1}&{2:.0f}&{3:.0f}-{4:.0f}&{5:.0f}-{6:.0f}".format(
                labels[key], units[key],
                round_to_n_sig_fig(res.ix[key][columns[0]], svar[key]),
                round_to_n_sig_fig(res.ix[key][columns[1]], svar[key]),
                round_to_n_sig_fig(res.ix[key][columns[2]], svar[key]),
                round_to_n_sig_fig(res.ix[key][columns[3]], svar[key]),
                round_to_n_sig_fig(res.ix[key][columns[4]], svar[key])) +
                "\\\\ \n" for key in svar.keys()]
        STR = STR + STR1
    if verbose:
        "writing file to {0}".format(path + prefix + ".tex")
    f = open(path + prefix + ".tex", "w")
    f.writelines(STR)
    f.close()
    return


# ------ newly added functions for comparing before and after adding prior--

def resample_to_size(x, n):
    """resample the array x to size n """
    x = np.array(x)
    ix = np.random.randint(0, x.size, n)
    return np.array([x[index] for index in ix])


def plot_perdiff(perdiff, labels, title):
    """
    Plot the percentage differences between different % of data within
    Will's   2 million iterations
    """
    fig = plt.figure()
    for i in range(5):
        plt.plot(range(1, len(perdiff[:, 1]) + 1),
                 perdiff[:, i], label=labels[i])
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc=4, bbox_to_anchor=(1.375, .5))
    minorLocator = MultipleLocator(2)

    plt.title(title)
    plt.grid()
    plt.grid(True, which='minor')
    plt.xlabel('% of iterations')
    plt.ylabel('Percent difference')
    ax.xaxis.set_minor_locator(minorLocator)
    plt.savefig(title + '.png', dpi=300, bbox_inches='tight')
    return fig


def prior_diff(data1, data2, prefix, prob1=None, prob2=None, N_bins=100,
               histrange=None, x_lim=None, y_lim=None, x_label=None,
               y_label=None, legend=None):
    """
    Plot the histograms between the same set of data with
    and without applying prior on the same plot
    """

    fig = pylab.figure()

    # bin data 1 first
    hist2, binedges2, tmp2 = \
        pylab.hist(data2, bins=N_bins, histtype='step',
                   weights=prob2, range=histrange, color='#ff0000', linewidth=2)
    hist1, binedges1, tmp1 = \
        pylab.hist(data1, bins=N_bins, histtype='step',
                   weights=prob1, range=histrange, color='#0000ff', linewidth=2)

    # Calculate the location and %confidence intervals for data 1
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned1 = numpy.ones(
                hist1[i]) * (binedges1[i] + binedges1[i + 1]) / 2
        else:
            x_temp1 = numpy.ones(
                hist1[i]) * (binedges1[i] + binedges1[i + 1]) / 2
            x_binned1 = numpy.concatenate((x_binned1, x_temp1))
    loc1 = biweightLoc(x_binned1)
    ll_68_1, ul_68_1 = bcpcl(loc1, x_binned1, 1)
    ll_95_1, ul_95_1 = bcpcl(loc1, x_binned1, 2)

    # Create location and confidence interval line plots
    # this should totally be replaced
    pylab.plot((loc1, loc1), (pylab.ylim()[0], pylab.ylim()[1]),
               '--', linewidth=2, color='#6495ed', label='$C_{BI}$')
    pylab.plot((ll_68_1, ll_68_1), (pylab.ylim()[0], pylab.ylim()[1]),
               '-.', linewidth=2, color='#7fffd4', label='68% $IC_{B_{BI}}$')
    pylab.plot(
        (ul_68_1,
         ul_68_1),
        (pylab.ylim()[0],
         pylab.ylim()[1]),
        '-.',
        linewidth=2,
        color='#7fffd4')
    pylab.plot((ll_95_1, ll_95_1), (pylab.ylim()[0], pylab.ylim()[1]),
               ':', linewidth=2, color='#87ceeb', label='95% $IC_{B_{BI}}$')
    pylab.plot(
        (ul_95_1,
         ul_95_1),
        (pylab.ylim()[0],
         pylab.ylim()[1]),
        ':',
        linewidth=2,
        color='#87ceeb')

    # Calculate the location and %confidence intervals for data 2
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned2 = numpy.ones(
                hist2[i]) * (binedges2[i] + binedges2[i + 1]) / 2
        else:
            x_temp2 = numpy.ones(
                hist2[i]) * (binedges2[i] + binedges2[i + 1]) / 2
            x_binned2 = numpy.concatenate((x_binned2, x_temp2))
    loc2 = biweightLoc(x_binned2)
    ll_68_2, ul_68_2 = bcpcl(loc2, x_binned2, 1)
    ll_95_2, ul_95_2 = bcpcl(loc2, x_binned2, 2)

    # Create location and confidence interval line plots
    pylab.plot((loc2, loc2), (pylab.ylim()[0], pylab.ylim()[1]),
               '--', linewidth=2, color='#ff4500', label='$C_{BI}$')
    pylab.plot((ll_68_2, ll_68_2), (pylab.ylim()[0], pylab.ylim()[1]),
               '-.', linewidth=2, color='#ff8c00', label='68% $IC_{B_{BI}}$')
    pylab.plot(
        (ul_68_2,
         ul_68_2),
        (pylab.ylim()[0],
         pylab.ylim()[1]),
        '-.',
        linewidth=2,
        color='#ff8c00')
    pylab.plot((ll_95_2, ll_95_2), (pylab.ylim()[0], pylab.ylim()[1]),
               ':', linewidth=2, color='#ffa500', label='95% $IC_{B_{BI}}$')
    pylab.plot(
        (ul_95_2,
         ul_95_2),
        (pylab.ylim()[0],
         pylab.ylim()[1]),
        ':',
        linewidth=2,
        color='#ffa500')

    # create labels for the plots
    if x_label is not None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label is not None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim is not None:
        pylab.xlim(x_lim)
    if y_lim is not None:
        pylab.ylim(y_lim)
    if legend is not None:
        pylab.legend()
    # fontsize=14
    #ax = pylab.gca()
    # for tick in ax.xaxis.get_major_ticks():
        # tick.label1.set_fontsize(fontsize)

    filename = prefix + '_prior_diff'
    pylab.savefig(title + '.png', dpi=300, bbox_inches='tight')

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(
        prefix, loc1, ll_68_1, ul_68_1, ll_95_1, ul_95_1)

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(
        prefix, loc2, ll_68_2, ul_68_2, ll_95_2, ul_95_2)

    return loc2, ll_68_2, ul_68_2, ll_95_2, ul_95_2


def prior_diff_pdf(data1, data2, prefix, prob1=None, prob2=None, N_bins=100,
                   histrange=None, x_lim=None, y_lim=None, x_label=None,
                   y_label=None, legend=None):
    """
    Plot the pdf between the same set of data with
    and without applying prior on the same plot
    """

    # bin data 2 first
    hist2, binedges2, tmp2 = \
        pylab.hist(data2, bins=N_bins, histtype='step', weights=prob2,
                   range=histrange, color='#ff0000', linewidth=2)
    # bin data 1
    hist1, binedges1, tmp1 = \
        pylab.hist(data1, bins=N_bins, histtype='step', weights=prob1,
                   range=histrange, color='#0000ff', linewidth=2)
    pylab.close()

    fig = pylab.figure()
    # plot the pdf for the data
    pdf2, histbin2, tmp_pdf2 = \
        pylab.hist(data2, bins=N_bins, normed=True,
                   histtype='step', weights=prob2, range=histrange,
                   color='#ff0000', linewidth=2)

    pdf1, histbin1, tmp_pdf1 = pylab.hist(data1, bins=N_bins, normed=True,
                                          histtype='step', weights=prob1,
                                          range=histrange, color='#0000ff',
                                          linewidth=2)

    # Calculate the location and %confidence intervals for data 1
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned1 = numpy.ones(histbin1[i]) *\
                (binedges1[i] + binedges1[i + 1]) / 2
        else:
            x_temp1 = numpy.ones(
                hist1[i]) * (binedges1[i] + binedges1[i + 1]) / 2
            x_binned1 = numpy.concatenate((x_binned1, x_temp1))
    loc1 = biweightLoc(x_binned1)
    ll_68_1, ul_68_1 = bcpcl(loc1, x_binned1, 1)
    ll_95_1, ul_95_1 = bcpcl(loc1, x_binned1, 2)

    # adjust the max ylim so it does not look weird
    ylim_max = pdf2.max() * 1.2

    # Create location and confidence interval line plots
    pylab.plot((loc1, loc1), (pylab.ylim()[0], ylim_max), '--', linewidth=2,
               color='#6495ed', label='$C_{BI}$')
    pylab.plot((ll_68_1, ll_68_1), (pylab.ylim()[0], ylim_max), '-.',
               linewidth=2, color='#7fffd4', label='68% $IC_{B_{BI}}$')
    pylab.plot((ul_68_1, ul_68_1), (pylab.ylim()[0], ylim_max), '-.',
               linewidth=2, color='#7fffd4')
    pylab.plot((ll_95_1, ll_95_1), (pylab.ylim()[0], ylim_max), ':',
               linewidth=2, color='#87ceeb', label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95_1, ul_95_1), (pylab.ylim()[0], ylim_max), ':',
               linewidth=2, color='#87ceeb')

    # Calculate the location and %confidence intervals for data 2
    # Since my location and confidence calculations can't take weighted data I
    # need to use the weighted histogram data in the calculations
    for i in numpy.arange(N_bins):
        if i == 0:
            x_binned2 = numpy.ones(
                hist2[i]) * (binedges2[i] + binedges2[i + 1]) / 2
        else:
            x_temp2 = numpy.ones(
                hist2[i]) * (binedges2[i] + binedges2[i + 1]) / 2
            x_binned2 = numpy.concatenate((x_binned2, x_temp2))
    loc2 = biweightLoc(x_binned2)
    ll_68_2, ul_68_2 = bcpcl(loc2, x_binned2, 1)
    ll_95_2, ul_95_2 = bcpcl(loc2, x_binned2, 2)

    # Create location and confidence interval line plots
    # if y_lim == None:
    # else:
    #    ylim_max = y_lim[1]
    #    ylim_min = y_lim[0]
    # if x_lim == None:
    #    xlim_max = pylab.xlim()[1]
    #    xlim_min = pylab.xlim()[0]
    # else:
    #    xlim_max = x_lim[1]
    #    xlim_min = x_lim[0]
    pylab.plot((loc2, loc2), (pylab.ylim()[0], ylim_max), '--', linewidth=2,
               color='#ff4500', label='$C_{BI}$')
    pylab.plot((ll_68_2, ll_68_2), (pylab.ylim()[0], ylim_max), '-.',
               linewidth=2, color='#ff8c00', label='68% $IC_{B_{BI}}$')
    pylab.plot((ul_68_2, ul_68_2), (pylab.ylim()[0], ylim_max), '-.',
               linewidth=2, color='#ff8c00')
    pylab.plot((ll_95_2, ll_95_2), (pylab.ylim()[0], ylim_max), ':',
               linewidth=2, color='#ffa500', label='95% $IC_{B_{BI}}$')
    pylab.plot((ul_95_2, ul_95_2), (pylab.ylim()[0], ylim_max), ':',
               linewidth=2, color='#ffa500')

    # create labels for the plots
    if x_label is not None:
        pylab.xlabel(x_label, fontsize=20)
    if y_label is not None:
        pylab.ylabel(y_label, fontsize=20)
    if x_lim is not None:
        pylab.xlim(x_lim)
    if y_lim is not None:
        pylab.ylim(y_lim)
    if legend is not None:
        pylab.legend()
    # fontsize=14
    #ax = pylab.gca()
    # for tick in ax.xaxis.get_major_ticks():
        # tick.label1.set_fontsize(fontsize)
    pylab.ylim(0, pdf2.max() * 1.2)

    filename = prefix + '_prior_diff'
    pylab.savefig(filename, dpi=300, bbox_inches='tight')

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(
        prefix, loc1, ll_68_1, ul_68_1, ll_95_1, ul_95_1)

    print '{0}, {1:0.4f}, {2:0.4f}, {3:0.4f}, {4:0.4f}, {5:0.4f}'.format(
        prefix, loc2, ll_68_2, ul_68_2, ll_95_2, ul_95_2)

    return loc2, ll_68_2, ul_68_2, ll_95_2, ul_95_2


def histplot2d_2contour(x1, y1, x, y, prefix, prob1=None, prob=None,
                        N_bins=100, histrange=None, x_lim=None, y_lim=None,
                        x_label=None, y_label=None, legend=None, ls1=None,
                        alpha1=1.0, tick_fontsize=14, label_fontsize=15):
    """
    My function for plotting 2 sets of 2d contour
    to create the confidence interval plot
    see plt.hist2d for detailed documentation

    x1 = numpy array, what's plotted on the x axis of the greyed out
        contour, or the values from default priors
    y1 = numpy array, what's plotted on the y axis of the greyed out
        contour, or the values from default priors
    x = numpy array, what's plotted on the x axis of the greyed out
        contour, or the values from additional priors
    y = numpy array, what's plotted on the y axis of the greyed out
        contour, or the values from additional priors
    ls1 = linestyle of contour plot 1 for default arguments
    prefix =
    histrange : array_like shape(2, 2), optional, default: None
         The leftmost and rightmost edges of the bins along each dimension
              (if not specified explicitly in the bins parameters):
                  [[xmin, xmax], [ymin, ymax]].
              All values outside of this range will be
                   considered outliers and not tallied in the histogram.
    """

    if histrange is not None and histrange.shape != (2, 2):
        histrange = \
            [[histrange[0], histrange[1]], [histrange[2], histrange[3]]]

    H1, xedges1, yedges1 = \
        numpy.histogram2d(x1, y1, bins=N_bins, weights=prob1, range=histrange)

    H1 = numpy.transpose(H1)
    # Flatten H
    h = numpy.reshape(H1, (N_bins ** 2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum / h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum / h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]

    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x1 = xedges1[:-1] + (xedges1[1] - xedges1[0]) / 2
    y1 = yedges1[:-1] + (yedges1[1] - yedges1[0]) / 2
    X1, Y1 = numpy.meshgrid(x1, y1)

    fig = pylab.figure()
    # Countours
    # imshow
    #im = pylab.imshow(H1,cmap=pylab.cm.gray)
    # pylab.pcolor(X1,Y1,H1,cmap=pylab.cm.white)

    if x_label is not None:
        pylab.xlabel(x_label, fontsize=label_fontsize)
    if y_label is not None:
        pylab.ylabel(y_label, fontsize=label_fontsize)
    if x_lim is not None:
        pylab.xlim(x_lim)
    if y_lim is not None:
        pylab.ylim(y_lim)
    # if legend is not None:
    # Dummy lines for legend
        # 800000 is for light gray - 68% 1 sigma
        # 0000A0 is for whitesmoke - 95% confidence 3 sigma
    # pylab.plot((0,1),(0,1),c='#f5fffa',linewidth=2,label=('68%'))
    # pylab.plot((0,1),(0,1),c='#d3d3d3',linewidth=2,label=('95%'))
        # pylab.legend(scatterpoints=1)
    fontsize = 20
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)

    # first contour
    pylab.contour(X1, Y1, H1, (h_2sigma, h_1sigma),
                  linewidths=(2, 2), colors=('#a4a4a4', '#6e6e6e'),
                  ls=ls1, alpha=alpha1)
    # second contour
    # Create the confidence interval plot for the second sets of contour
    H, xedges, yedges = \
        numpy.histogram2d(x, y, bins=N_bins, range=histrange, weights=prob)

    H = numpy.transpose(H)
    # Flatten H
    h = numpy.reshape(H, (N_bins ** 2))
    # Sort h from smallest to largest
    index = numpy.argsort(h)
    h = h[index]
    h_sum = numpy.sum(h)
    # Find the 2 and 1 sigma levels of the MC hist
    for j in numpy.arange(numpy.size(h)):
        if j == 0:
            runsum = h[j]
        else:
            runsum += h[j]
        if runsum / h_sum <= 0.05:
            # then store the value of N at the 2sigma level
            h_2sigma = h[j]
        if runsum / h_sum <= 0.32:
            # then store the value of N at the 1sigma level
            h_1sigma = h[j]
    # Create the contour plot using the 2Dhist info
    # define pixel values to be at the center of the bins
    x = xedges[:-1] + (xedges[1] - xedges[0]) / 2
    y = yedges[:-1] + (yedges[1] - yedges[0]) / 2
    X, Y = numpy.meshgrid(x, y)

    # Coutours
    pylab.contour(X, Y, H, (h_2sigma, h_1sigma), linewidths=(2, 2))
    # imshow
    #im = pylab.imshow(H,cmap=pylab.cm.Blues)
    pylab.pcolor(X, Y, H, cmap=pylab.cm.gray_r)

    if x_label is not None:
        pylab.xlabel(x_label, fontsize=label_fontsize)
    if y_label is not None:
        pylab.ylabel(y_label, fontsize=label_fontsize)
    if x_lim is not None:
        pylab.xlim(x_lim)
    if y_lim is not None:
        pylab.ylim(y_lim)

    if legend is not None:
        # Dummy lines for legend
        # 800000 is for Maroon color - 68% 1 sigma
        # 0000A0 is for blue color - 95% confidence 3 sigma
        pylab.plot((0, 1), (0, 1), c='#87cefa', linewidth=2, label=('68%'))
        pylab.plot((0, 1), (0, 1), c='#6495ed', linewidth=2, label=('95%'))
        pylab.legend(scatterpoints=1)
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)

    filename = prefix + '2contour2d.png'
    pylab.savefig(filename, dpi=300, bbox_inches='tight')

    return fig


def plot_perdiff(perdiff, labels, title):
    """
    Plot the percentage differences between different % of data within
    Will's   2 million iterations
    """
    fig = plt.figure()
    for i in range(5):
        plt.plot(range(1, len(perdiff[:, 1]) + 1),
                 perdiff[:, i], label=labels[i])
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc=4, bbox_to_anchor=(1.375, .5))
    minorLocator = MultipleLocator(2)

    plt.title(title)
    plt.grid()
    plt.grid(True, which='minor')
    plt.xlabel('% of iterations')
    plt.ylabel('Percent difference')
    ax.xaxis.set_minor_locator(minorLocator)
    plt.savefig(title + '.png', dpi=300, bbox_inches='tight')
    return fig


def hist_scenario(this_ax, data1, data2, data1_mean, data2_mean, beta, i,
                  r95low, r95up, weight=None, cluster_label="NW",
                  beta_label_loc='left'):
    """plot in multipanel probability of different scenarios based on
    different assumptions, see the TSC_CM_different_prior-no approx.ipynb
    for usage

    parameters
    =========
    this_ax = matplotlib ax object
    data1 = np.array, pdf of the projected outgoing relic location
    data2 = np.array, pdf of the projected returning relic location

    returns
    =======
    prob2 = returning probability within r95low and r95up
    prob1 = outgoing probability within r95low and r95up

    note
    ====
    right now I just try to pick bin width to be as small as possible
    and add up all the bins for computing the probability. This should be
    improved if the data quality ever gets better.
    """
    n1, bins1, patches1 = \
        this_ax.hist(data1, histtype='step', bins=500, label=r'$s-outgoing$',
                     normed=True, weights=weight)
    n2, bins2, patches2 = \
        this_ax.hist(data2, histtype='step', label=r'$s-returning$',
                     bins=1000, normed=True, weights=weight)
    this_ax.axvspan(r95low, r95up, color='red',
                    label=r'true $s_{proj}$', alpha=0.3)
    # this_ax.axvline(NW_subclstr_sep.value, color='orange',
    #                label=r'obs. NW sep')
    this_ax.set_xlim(-0.1, 1.7)
    this_ax.set_ylabel('PDF')
    xmin, xmax = this_ax.get_xlim()
    ymin, ymax = this_ax.get_ylim()

    if beta_label_loc == "left":
        beta_label_xloc = xmin + .05
        beta_label_yloc = ymax - 0.2
    else:
        beta_label_xloc = xmax - .3
        beta_label_yloc = ymax - 0.2

    this_ax.text(beta_label_xloc, beta_label_yloc,
                 # r"Assume $<v_{relic}> / " +
                 # r"v_{{3D, {0}}}(t_{{col}})".format(cluster_label) +
                 r"$\beta =$ {0}".format(beta), va='top', size=15)

    # to compute the Bayes factor
    low_nbin1 = 0
    up_nbin1 = 0

    # sanity check to make sure that the red area is within the
    # bounds of the histograms
    assert bins1[len(bins1) - 1] > r95low, \
        "beta = {0} all hist1 are smaller than s_proj".format(beta)
    assert bins1[0] < r95up, \
        "beta = {0} all hist1 are bigger than s_proj".format(beta)

    # find bins that correspond to the red area of the curve
    while(bins1[low_nbin1] < r95low):
        low_nbin1 += 1

    if bins1[len(bins1) - 1] > r95up:
        while(bins1[up_nbin1] < r95up):
            up_nbin1 += 1
    else:
        up_nbin1 = len(bins1) - 1
        print "beta = {0} s_proj is bigger than some hist1".format(beta)

    assert low_nbin1 < up_nbin1, "something wrong with bins "

    low_nbin2 = 0
    up_nbin2 = 0

    # sanity check to make sure that the red area is within the
    # bounds of the histograms
    assert bins2[len(bins2) - 1] > r95low, \
        "beta ={0} all hist2 are smaller than s_proj".format(beta)
    assert bins2[0] < r95up, \
        "beta = {0} all hist2 are bigger than s_proj".format(beta)

    while(bins2[low_nbin2] < r95low):
        low_nbin2 += 1

    if bins2[len(bins2) - 1] > r95up:
        while(bins2[up_nbin2] < r95up):
            up_nbin2 += 1
    else:
        up_nbin2 = len(bins2) - 1
        print "beta = {0} s_proj is bigger than some hist2".format(beta)

    assert low_nbin2 < up_nbin2, \
        "beta = {0} something wrong with bins".format(beta)

    # we have fixed bin width for each set of histogram
    # compute probability roughly by adding the area of the bins
    bw1 = bins1[1] - bins1[0]
    bw2 = bins2[1] - bins2[0]
    prob1 = np.sum(n1[low_nbin1: up_nbin1]) * bw1
    prob2 = np.sum(n2[low_nbin2: up_nbin2]) * bw2
    # print "Beta = {0}, P_1 / P_0 = {1:1.4f} / {2:1.4f}".format(
    #    beta, prob2, prob1) + \
    #    " Bayes factor = {0:1.2f}".format(prob2 / prob1)
    # this_ax.text(xmax - 0.8, ymax - 1.75,
    #             "Bayes factor = "
    #             " {0:.2f}".format(prob2 / prob1), va='top', size=14)

    return prob2, prob1


# def radio_dist_prior(d_3D, mask, d_3Dmax = 3.0, d_3Dmin = 1.0):
#    '''
#    Create the prior from the radio relic constraints
#    we know that if the two subclusters are within 0.5 Mpc to 1.5 Mpc
#    then the detection of a radio relic is possible
#    d_3d is in Mpc
#    return a mask that can be applied to other data arrays
#    also return how many non-zero entries
#    to be modified such that this takes in the range for the uniform prior
#    input:
#    d_3D = numpy array to be masked
#    radiomask = numpy array that is contains either 1 or 0
#    d_3Dmax = float, the upper limit to be masked out
#    d_3Dmin = float, the lower limit to be masked out
#    '''
#    count = 0
#    for r in range(len(d_3D)):
#        if (d_3D[r]>d_3Dmax)  or (d_3D[r]< d_3Dmin):
#            radiomask[r] = 0
#            count += 1
#    count = len(radiomask) - count
#    return radiomask, count
