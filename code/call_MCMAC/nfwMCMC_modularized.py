from __future__ import division
import numpy
import numpy as np
import pylab
import cosmo
import tools
#import sys
import pickle
import time
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.special import erf
#from astropy import wcs

# Will 's modules
from profiles import nfw_Sigma
from profiles import nfw_Sigmabar
from profiles import nfwparam

# use cPickle to store debugging info may change to use HDF5 instead
import cPickle

"""
nfwMCMC is a program that simultaniously fits multiple NFW profiles to the
Weak Lensing data of in input catalog of galaxy shapes (e1, e2) and redshifts
(z).  This is very similar to lensMCMC_v0.0.py, except that instead of fitting
for both rho_s and r_s we are only fitting for r_s, essentially.
Techically this program fits for M_200 of each halo, which is realted to
rho_s and r_s through M_200/c emperical scaling relations (e.g. Duffy et
al. 2008) and general properties of the NFW profile (see e.g. Springer's
text).  More details are
contained in my OneNote notebook Programs>Miscelanious Programs>nfwMCMC.
"""

##########################################################################
# Functions
##########################################################################


def calc_kappa(del_c, r_s, r, z_halo, invSigmacr,
               h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0):
    '''
    del_c = characteristic overdensity of the CDM halo
    r_s = scale radius of the halo (Mpc)
    r = radius of interest (Mpc)
    z_halo = halo redshift
    h_scale = hubble scale H = h*100 km/s/Mpc
    Om = matter energy density
    Ol = dark energy density
    Or = radiation energy density

    Read Umetsu 2010 Table 1 and p. 21 for details
    '''
    Sigma = nfw_Sigma(del_c, r_s, r, z_halo, h_scale, Om, Ol, Or)
    return Sigma * invSigmacr, Sigma


def calc_real_gamma(del_c, r_s, r, z_halo, invSigmacr, Sigma,
                    h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0,
                    debug=False, prefix=None):
    """compute gamma in vectorized way
    """
    return (nfw_Sigmabar(del_c, r_s, r, z_halo, h_scale, Om, Ol, Or)
            - Sigma) * invSigmacr


def correct_reduced_shear(del_c, r_s, r, z_halo, invSigmacr,
                          h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0, debug=False,
                          prefix=None):
    """ apply additional corrections for El Gordo
    kappa = np array - of shape ngal * ncluster

    Stability : to be debugged and tested
    """
    kappa, Sigma = \
        calc_kappa(del_c, r_s, r, z_halo, invSigmacr, h_scale, Om, Ol, Or)
    real_gamma = calc_real_gamma(del_c, r_s, r, z_halo, invSigmacr, Sigma)

    # the following should produce a np array with
    # g.shape =  (n_gal, ncluster)
    g = real_gamma / (1 - np.sum(kappa))

    ## was not adding components of shear, i.e. the denominator correctly
    ## for the two subclusters
    # g = gamma / (1 - kappa)
    if numpy.sum(numpy.isnan(g)) != 0:
        raise ValueError('shear: There is a problem, ' +
                         ' a nan gamma value was calculated')

    # correcting for how reduced shear g' is related to true shear g
    # see Seitz & Schneider 1997
    # should only enable this when doing El Gordo analysis
    # apply corrections
    g = g * (1. + 0.79 * kappa)

    # g always has to be smaller than one in magnitude,
    # flip numerater with denominator
    thisMask = numpy.abs(g) > 1.0
    g[thisMask] = 1. / g[thisMask]

    if debug is True:
        debugFile = "debug_shear"
        if prefix is not None:
            debugFile += prefix
        print "dumping local variables to {0}.pkl".format(debugFile)
        cPickle.dump(locals(), open(debugFile + ".pkl", "w"))

    return g


def shear(del_c, r_s, r, z_halo, invSigmacr,
          h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0, debug=False, prefix=None):
    '''
    del_c = characteristic overdensity of the CDM halo
    r_s = scale radius of the halo (Mpc)
    r = radius of interest (Mpc)
    z_halo = halo redshift
    h_scale = hubble scale H = h*100 km/s/Mpc
    Om = matter energy density
    Ol = dark energy density
    Or = radiation energy density

    Read Umetsu 2010 Table 1 and p. 21 for details
    '''
    Sigma = nfw_Sigma(del_c, r_s, r, z_halo, h_scale, Om, Ol, Or)
    kappa = Sigma * invSigmacr
    # this is really an approximate expression for g+
    # see eqn near the end of second paragraph on p.21 of Umetsu
    g = (nfw_Sigmabar(del_c, r_s, r, z_halo, h_scale, Om, Ol, Or) - Sigma) * \
        invSigmacr / (1 - kappa)
    if numpy.sum(numpy.isnan(g)) != 0:
        raise ValueError('shear: There is a problem, ' +
                         ' a nan gamma value was calculated')

    # correcting for how reduced shear g' is related to true shear g
    # see Seitz & Schneider 1997
    # should only enable this when doing El Gordo analysis
    # apply corrections
    g = g * (1. + 0.79 * kappa)

    # g always has to be smaller than one in magnitude,
    # flip numerater with denominator
    thisMask = numpy.abs(g) > 1.0
    g[thisMask] = 1. / g[thisMask]

    if debug is True:
        debugFile = "debug_shear"
        if prefix is not None:
            debugFile += prefix
        print "dumping local variables to {0}.pkl".format(debugFile)
        cPickle.dump(locals(), open(debugFile + ".pkl", "w"))

    return g


def logLcalc(
        M_200, e1, e2, var_e1, var_e2, del_r, invSigmacr, cos2phi, sin2phi,
        z_halo, mask, h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0, debug=False,
        prefix=None):
    '''
    Calculates the likelihood of the NFW models given the ellipticity data.
    Same inputs as mcmcengine.
    M_200 = [list of length N_halo; units:1e14 M_sun] mass of each of the
    halos
    '''
    N_halo = numpy.shape(del_r)[1]
    N_gal = numpy.shape(del_r)[0]
    gamma = numpy.zeros((N_gal, N_halo))
    # Calculate the NFW parameters corresponding to M_200
    del_c, r_s = nfwparam(M_200, z_halo, h_scale=h_scale, Om=Om, Ol=Ol, Or=Or)

    for h in numpy.arange(N_halo):
        # loop through each halo calculating the expected absolute shear
        # due to that individual halo
        # *# Could place a mapping based multiprocessing function here
        # calculate the expected ellipticity components of each galaxy
        # OLD SHEAR FUNCTION
        gamma[mask[:, h], h] = \
            shear(del_c[h], r_s[h], del_r[mask[:, h], h], z_halo[h],
                  invSigmacr[mask[:, h], h], h_scale, Om, Ol, Or,
                  debug=debug, prefix=prefix)

       # gamma[mask[:, h], h] = \
       #     correct_reduced_shear(del_c[h], r_s[h], del_r[mask[:, h], h],
       #     z_halo[h], invSigmacr[mask[:, h], h], h_scale, Om, Ol, Or,
       #     debug=debug, prefix=prefix)
    mask_e = numpy.sum(mask, axis=1) != 0
    # calculate the expected ellipticities of each galaxy due to all halos
    e1_exp = numpy.sum(-gamma[mask_e, :] * cos2phi[mask_e, :], axis=1)
    e2_exp = numpy.sum(-gamma[mask_e, :] * sin2phi[mask_e, :], axis=1)

    # calculate the chi^2
    chi2 = numpy.sum((e1[mask_e] - e1_exp) ** 2 / var_e1[mask_e]) +\
        numpy.sum((e2[mask_e] - e2_exp) ** 2 / var_e2[mask_e])
    # chi2 = numpy.sum((e1-e1_exp)**2/var_e1)+numpy.sum((e2-e2_exp)**2/var_e2)
    # calculate the likelihood of the starting position
    # *# Note that this liklihood calculations is not really
    # correct there should be
    # *# another term added since there are errors in the expected
    # ellipticities
    logL = -chi2 / 2.

    if numpy.isnan(logL):
        raise ValueError('logLcalc problem: a nan likelihood was calculated')

    if debug is True:
        if prefix is None:
            debugFile = "debug_logLcal"
        else:
            debugFile = "debug_logLcal" + prefix
        print "dumping local variables to {0}.pkl".format(debugFile)
        cPickle.dump(locals(), open(debugFile + ".pkl", "w"))

    return logL


def mcmcengine(
        e1, e2, var_e1, var_e2, del_r, invSigmacr, cos2phi, sin2phi, z_halo,
        mask, M_200_start, parambounds, cpropdist, N_check, chainset, N_burn,
        del_a, del_d, halos, centroid_step,
        cat, key, coord, r_bounds, cd,
        h_scale=0.7, Om=0.3, Ol=0.7, Or=0, verbose=True, debug=False,
        prefix=None, float_centroid=False):
    '''
    Input:
    e1 = [array, size(N_galaxy)] measured e1 ellipticity component of each
        galaxies
    e2 = [array, size(N_galaxy)] measured e2 ellipticity component of each
        galaxies
    var_e1 = [array, size(N_galaxy,N_halo)] essentially weight to apply to
        deviation of measured e1 from <e1>
    var_e2 = [array, size(N_galaxy,N_halo)] essentially weight to apply to
        deviation of measured e1 from <e2>
    del_r = [array, size(N_galaxy,N_halo)] radial distance of each galaxy
            from the respective halo
    invSigmacr = [array, size(N_galaxy,N_halo)] (Sigma_critical)^(-1) of each
        galaxy with respect to each halo
    cos2phi = [array, size(N_galaxy,N_halo)]
    sin2phi = [array, size(N_galaxy,N_halo)]
    z_halo = [array, size(N_halo)] redshift of each halo
    mask = [array, size(N_galaxy, N_halo)] mask of only valid galaxies
    M_200_start = [array, size(N_halo)] M_200 starting values for each halo
    parambounds = [array, size(2*N_halo)]
    cpropdist = [float] [sigma_M_200] standard deviations of the
        normal proposal distributions in unit of 1e14 Msun
    N_check = [integer]
        Number of successful steps to take in this section of
        chain
    N_burn = integer
        number of chain before burn in cut off
    Output:
    M_200_array = [array N_check long] and contains all the successful steps
    chainset = integer, that represents the N-th chain set that we are
                running the MCMC, we will tweak the acceptance rate for the
                0th chain set only
    '''

    ### helper function for debugging
    if debug is True:
        if prefix is None:
            debugFile = "debug_mcmcengine"
        else:
            debugFile = "debug_mcmcengine" + prefix
        print "dumping local variables to {0}.pkl".format(debugFile)
        cPickle.dump(locals(), open(debugFile + ".pkl", "w"))

    N_iter = 0
    N_halo = numpy.size(z_halo)
    N_gal = numpy.size(e1)

    # check how many draws of the halos get stuck
    stuckTimes = np.zeros(N_halo)

    M_200_old = M_200_start
    centroid_old = halos

    # record the mass for making trace plots
    trace_mass_NW = []
    trace_mass_SE = []

    # record acceptance rate
    accept = .5

    # record global number of successful steps and failed steps
    success = 0
    fail = 0

    # record the number of failed attempts to jump to a new value
    consec_failtimes = 0

    # create the blank arrays for recording the masses
    M_200 = numpy.zeros((N_check, N_halo))

    # create the blank arrays for recording the centroids
    # the shape is according to the halo specification
    # wastes a lot of space for the redshift entry but oh well
    centroids = numpy.zeros((N_check, N_halo * 3))

    # Calculate the likelihood of the starting (old) model
    logL_old = logLcalc(M_200_start, e1, e2, var_e1, var_e2, del_r, invSigmacr,
                        cos2phi, sin2phi, z_halo, mask, h_scale, Om, Ol,
                        Or, debug=debug, prefix=prefix)

    while N_iter < N_check:
        # determine the next random step
        M_200_new = numpy.random.randn(N_halo) * cpropdist + M_200_old.copy()

        # check if the new step is outside the parameter prior distribution
        # if so then redraw a new random step
        # parallelizable if slow
        for h in range(N_halo):
            while(outsideBound(h, M_200_new, parambounds)):
                M_200_new[h], stuckTimes[h] = \
                    draw_new_mass_for_halo_h(
                        h, M_200_new, stuckTimes, M_200_old, cpropdist, N_iter,
                        parambounds, verbose)

        if float_centroid:
            # initialize an array to store the drawn step
            # we would be drawing a pair of centroid location each time
            centroid_new = \
                get_new_centroids(N_halo, centroid_old, centroid_step)

            # calculate the new del_r, del_a, sigma etc. due to updating the
            # centroid
            invSigmacr, del_a, del_d = \
                calc_stuff(N_halo, N_gal, del_a, del_d, cat, h_scale, Om, Ol,
                           centroid_new, key, coord, invSigmacr, verbose)

            # mask the galaxies according to new centroids
            mask = \
                mask_galaxies(cat, N_gal, N_halo, del_r,
                              r_bounds, halos, key, coord)

            cos2phi, sin2phi = compute_rotation(cd, del_a, del_d)

        # calculate the likelihood of the new model
        logL_new = \
            logLcalc(M_200_new, e1, e2, var_e1, var_e2, del_r,
                     invSigmacr, cos2phi, sin2phi, z_halo, mask, h_scale,
                     Om, Ol, Or)

        # this is the Metropolis algorithm
        # since we are using symmetric proposal functions
        if logL_new >= logL_old:
            logL_old = logL_new
            M_200[N_iter, :] = M_200_new.copy()
            M_200_old = M_200_new.copy()
            N_iter += 1
            # consec_failtimes is Will 's l, reset here since we succeeded
            consec_failtimes = 0
            success += 1
            if N_iter % 100 == 0 and N_iter != 0:
                accept = success / (success + fail)
                print 'step {0} / {1} completed'.format(N_iter, N_check) +\
                    ' ,percentage of acceptance is ' + \
                    '{0:.2f}%'.format(accept * 100)
        else:
            logProb = numpy.log(numpy.random.rand())
            if logProb < (logL_new - logL_old):
                logL_old = logL_new
                M_200[N_iter, :] = M_200_new.copy()
                M_200_old = M_200_new.copy()
                N_iter += 1
                success += 1
                if N_iter % 100 == 0 and N_iter != 0:
                    accept = success / (success + fail)
                    print 'step {0} / {1} completed'.format(N_iter, N_check) +\
                        ' ,percentage of acceptance is ' + \
                        '{0:.2f}%'.format(accept * 100)
            else:
                consec_failtimes += 1
                fail += 1
                if consec_failtimes > int(N_check * 0.8):
                    print_little_success(N_iter, consec_failtimes, logProb,
                                         logL_new, logL_old)

        # record the mass for trace plots
        #trace_mass_NW.append(M_200_old[0])
        #trace_mass_SE.append(M_200_old[1])

        # tweak the propdist to have a desired acceptance rate
        # before the burn in
        #if N_iter % 100 == 0 and N_iter != 0:
        #    cpropdist = \
        #        adjust_accept_rate(N_iter, cpropdist,
        #                           accept, N_burn, N_check, chainset)

    print 'mcmcengine: finished {0} chainlinks'.format(N_check)

    return M_200, trace_mass_NW, trace_mass_SE, accept * 100, cpropdist


def bcpcl(T, T_p, N_sigma):
    '''
    Calculates the bias corrected percent confidence limits.
    -- Suppose that we have observed data (y1, y2, ..., yn) and use it to
    estimate a population parameter Q (e.g. Q could be the true mean of the
    entire population).
    -- T is a statistic that estimates Q. For example T could be an estimate
    of the true mean by calculating the mean of  (y1, y2, ..., yn).
    -- Suppose that we create m bootstrap samples (y_p_1j, y_p_2j,...,j_p_nj)
    from observed sample  (y1, y2, ..., yn), where j is the jth bootstrap
    sample.
    -- Then T_p_j is the jth bootstrap observation of T.
    For example this could be the mean of (y_p_1j, y_p_2j, ...,j_p_nj).

    T = [float] e.g. biweight Location for (y1, y2, ..., yn)
    T_p = [vector array] biwieght Locations for the bootstrap samples
    N_sigma = the number of sigma to report the confidence limits for
        e.g. for 95% confidence limits N_sigma=2
    Return (lower, upper) confidence limits

    -Note that this was taken from CAT.py and is not specifically intended for
    application to MCMC analysis but rather bootstrap resampling analysis but
    I think that the general concept of its application to MCMC datasets should
    be valid.
    '''
    # Percentile confidence interval is defined as 100%(1-a), thus for 1sigma
    # a=0.32
    a = 1 - erf(N_sigma / numpy.sqrt(2))
    # order the bootstrap sample values smallest to largest
    index = numpy.argsort(T_p)
    T_p = T_p[index]
    # Number of bootstrap samples
    m = numpy.size(T_p)
    # Calculate the bias correction term
    mask = T_p < T
    z_0 = norm.ppf(numpy.sum(mask) / m)
    # Calculate the a1 and a2 values
    a1 = norm.cdf(2 * z_0 + norm.ppf(a / 2))
    a2 = norm.cdf(2 * z_0 + norm.ppf(1 - a / 2))
    # Calculate the lower and upper indicies of lower and upper confidence
    # intervals
    id_L = numpy.int(m * a1) - 1
    id_U = numpy.int(m * a2)
    # Find the lower an upper confidence values
    T_L = T_p[id_L]
    T_U = T_p[id_U]

    return T_L, T_U


def chi2calc(e1, e2, var_e1, var_e2, del_r, invSigmacr, cos2phi, sin2phi,
             z_halo, mask, M_200, h_scale=0.7, Om=0.3, Ol=0.7, Or=0,
             verbose=True):
    '''
    This function calculates and reports the chi^2 and reduced chi^2 of a
    specific population of halos with given M_200.  It follows the same basic
    calculations as logLcalc.
    '''
    N_halo = numpy.shape(del_r)[1]
    N_gal = numpy.shape(del_r)[0]
    gamma = numpy.zeros((N_gal, N_halo))
    # Calculate the NFW parameters corresponding to M_200
    # the units used here is entirely in terms of M_sun NOT e14 M_sun
    del_c, r_s = nfwparam(M_200, z_halo, h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0)
    for h in numpy.arange(N_halo):
        # loop through each halo calculating the expected absolute shear
        # due to that individual halo
        # *# Could place a mapping based multiprocessing function here
        # calculate the expected ellipticity components of each galaxy
        gamma[mask[:, h], h] = shear(del_c[h], r_s[h], del_r[mask[:, h], h],
                                     z_halo[h], invSigmacr[mask[:, h], h],
                                     h_scale, Om, Ol, Or)
    mask_e = numpy.sum(mask, axis=1) != 0
    # calculate the expected ellipticities of each galaxy due to all halos
    e1_exp = numpy.sum(-gamma[mask_e, :] * cos2phi[mask_e, :], axis=1)
    e2_exp = numpy.sum(-gamma[mask_e, :] * sin2phi[mask_e, :], axis=1)

    # calculate the chi^2
    chi2 = numpy.sum((e1[mask_e] - e1_exp) ** 2 / var_e1[mask_e]) +\
        numpy.sum((e2[mask_e] - e2_exp) ** 2 / var_e2[mask_e])
    # calculate the reduced chi^2
    N_data = numpy.sum(mask_e)
    nu = N_data - N_halo
    redchi2 = chi2 / nu
    return N_data, chi2, redchi2


def fitM200(halos, catalog, coord, ellip_meas, ellip_intstd, r_bounds,
            parambounds, N_maxlinks, N_chain, N_burn, N_check, propdist,
            prefix, centroid_step=None,
            float_centroid=False,
            seedbounds=None, bins=20, halo_names=None,
            cd=((-1, 0), (0, 1)),
            h_scale=0.7, Om=0.3, Ol=0.7, Or=0.0, verbose=True, debug=False):
    ##########################################################################
    # PROGRAM
    ##########################################################################
    # Read in the catalog and header information
    print '============================================================'
    print 'WARNING: Check nfwMCMC.shear(): special correction done for' +\
        ' El Gordo'
    print 'disable if you are doing analysis for another cluster'
    print '============================================================'
    cat = tools.readcatalog(catalog)
    key = tools.readheader(catalog)
    N_gal = numpy.shape(cat)[0]  # number of galaxies

    propdist = np.ones(N_chain) * propdist

    # Determine the number of halo's to model
    N_halo = numpy.size(halos) // 3

    # Perform some basic user input checks
    input_checks(N_check, N_burn, halos, halo_names, N_chain, parambounds,
                 seedbounds, N_halo)

    if float_centroid is True:
        assert centroid_step is not None, "centroid_step is not specified"

    # Don't allow seedbound to exceed parambounds
    for i in range(N_chain):
        check_seedbound_against_parambounds(parambounds, seedbounds, i)

    ###
    # Calculate the properties of the galaxies with respect to each halo
    ###
    print ('fitM200: calculating properties of the galaxies with ' +
           'respect to each halo')

    invSigmacr, del_a, del_d  = \
        calc_stuff(N_halo, N_gal, cat, h_scale, Om, Ol,
                   halos, key, coord, verbose)

    ##### See wikipedia
    ##http://upload.wikimedia.org/math/d/a/0/
    ##da0c0f7924ada281d56e483d2ecd5427.png
    ### the following made use of small angle approximate of
    ### tan (del_r) \approx del_r
    # note del_r is in units of Mpc
    del_r = numpy.sqrt(del_a ** 2 + del_d ** 2)

    mask = \
        mask_galaxies(cat, N_gal, N_halo, del_r, r_bounds, halos, key, coord)

    # Find the unit relation between del_x, del_y and del_a, del_d
    cos2phi, sin2phi = compute_rotation(cd, del_a, del_d)


    """
    Calculate the combined variance for each galaxy's ellipticity components
    *# Karen: this is also what James Jee uses for his code
    *# Will: I am not sure if the following is correct since the last part
    has to do with
    *# and error on the expected ellipticity not on the measeured like
    the firstinSigmacr
    *# two parts. See Andi Mahdavi's Davis talk.
    """
    var_e1 = \
        cat[:, key[ellip_meas[2]]] ** 2. + ellip_intstd ** 2.
        # +var_invSigmacr

    var_e2 = \
        cat[:, key[ellip_meas[3]]] ** 2. + ellip_intstd ** 2.
        # +var_invSigmacr


    ###
    # Initiate the chains
    ###
    """
    Each of the N_chains will be run in parallel so that we can perform
    convergence checks every so often
    Create the blank parameter arrays
    """
    M_200 = numpy.zeros((N_maxlinks, N_halo, N_chain))

    ####del_c = numpy.zeros((N_maxlinks, N_halo, N_chain))
    ####r_s = numpy.zeros((N_maxlinks, N_halo, N_chain))
    # Extract each halos' redshift information and create a redshift array
    z_halo = numpy.zeros(N_halo)
    for h in numpy.arange(N_halo):
        z_halo[h] = halos[h * 3 + 2]

    if debug is True:
        debugFile = "debug_fitM200_" + prefix
        print "dumping local variables to {0}.pkl".format(debugFile)
        cPickle.dump(locals(), open(debugFile + ".pkl", "w"))

    # Create the chains
    chainlink = 0
    while chainlink < numpy.floor(N_maxlinks / N_check):
        t_start = time.time()
        if verbose:
            print 'fitM200: creating {0} '.format(N_chain * N_check) + \
                'chainlinks, set {0} '.format(chainlink)
        if chainlink == 0:
            # Draw a starting location for each chain
            if seedbounds is None:
                # then use the parameter bounds to limit the initial random
                # step
                print 'using parambounds as starting mass'
                M_200[chainlink, :, :] = \
                    numpy.random.rand(N_halo, N_chain) *\
                    (parambounds[1] - parambounds[0]) + parambounds[0]
            else:
                # check that the seed bound input has the correct number of
                # inputs
                if numpy.size(seedbounds) / N_chain != 2:
                    raise ValueError(
                        'fitM200: incorrect number of parameters ' +
                        'entered for seedbounds, exiting')

                for chain in numpy.arange(N_chain):
                    M_200[chainlink, :, chain] = \
                        numpy.random.rand(N_halo) * \
                        (seedbounds[chain * 2 + 1] - seedbounds[chain * 2]) +\
                        seedbounds[chain * 2]
            # if chainlink!=0 the last step will be used as the starting point
            # for the next
            # section of chain

        # Run through each chain N_check at a time then compare for convergence
        # *# This would be an ideal location to take advantage of
        # multiprocessors
        for chain in numpy.arange(N_chain):
            if verbose:
                print 'fitM200: creating {0} '.format(N_check) + \
                    'chainlinks of chain {0}'.format(chain)
            # KN: WHY ????
            if chainlink == 0:
                M_200_start = M_200[chainlink * N_check, :, chain]
            else:
                M_200_start = M_200[chainlink * N_check - 1, :, chain]
            # print 'M_200_start is ', M_200_start
            # Call mcmcengine to determine the next N_check steps 1 chain
            # at a time
            # M_200[chainlink*N_check:(chainlink+1)*N_check,:,c]=\
            #                mcmcengine(cat[:,key[ellip_meas[0]]],
            #                cat[:,key[ellip_meas[1]]], var_e1, var_e2,
            #                del_r, invSigmacr, cos2phi, sin2phi, z_halo,
            #                mask, M_200_start, parambounds,
            #                propdist, N_check, h_scale, Om, Ol, Or,verbose)

            # this function is too complicated :/
            M_200[chainlink * N_check: (chainlink + 1) * N_check, :,
                  chain], trace_NW, trace_SE, accept, propdist[chainlink] =\
                mcmcengine(cat[:, key[ellip_meas[0]]],
                           cat[:, key[ellip_meas[1]]],
                           var_e1, var_e2,
                           del_r, invSigmacr, cos2phi, sin2phi, z_halo,
                           mask, M_200_start, parambounds,
                           propdist[i], N_check, chain, N_burn, del_a,
                           del_d, halos, centroid_step,
                           cat, key, coord, r_bounds, cd,
                           h_scale, Om, Ol, Or, verbose, debug, prefix)

            # want to plot the trace plots only after convergence
            # testing the trace plot list
            # have to think about how to do this well
            #trace_plot(c, accept, trace_NW, trace_SE, N_burn, prefix)

        ###
        # Test of Convergence of Chains
        ###
        # test for convergence using Gelmans and Ruben test ref. Cos. on Beach
        # Verde lecture 2, see also Convergence st Chains onenote note
        # calculate the mean of each chain
        if verbose:
            print 'fitM200: testing if chains have converged after ' + \
                '{0} chainlinks'.format((chainlink + 1) * N_check)
        N_links = (chainlink + 1) * N_check - N_burn

        # calculate the mean of each chain
        mean_M_200_chain = \
            numpy.sum(M_200[N_burn: (chainlink + 1) * N_check, :, :],
                      axis=0) / N_links
        # calculate the mean of the distribution
        mean_M_200_dist = numpy.sum(mean_M_200_chain, axis=1) / N_chain
        # Calculate the vairance between chains
        B_M_200 = \
            numpy.sum((mean_M_200_chain -
                       numpy.reshape(mean_M_200_dist, (N_halo, 1))) ** 2,
                      axis=1) / (N_chain - 1)

        # Calculate the variance within chains
        W_M_200 = \
            numpy.sum(
                numpy.sum((M_200[N_burn: (chainlink + 1) * N_check, :, :] -
                           mean_M_200_chain) ** 2, axis=0), axis=1) /\
            (N_chain * (N_links - 1))

        # Calculate the estimated parameter variance
        V_M_200 = \
            (N_links - 1) / N_links * W_M_200 + B_M_200 * (1 + 1 / N_chain)

        # Calculate the Gelman-Rubin statisitc
        R_M_200 = V_M_200 / W_M_200
        # if all the R values are < 1.03 then chains have converged and can
        # exit
        if numpy.sum(R_M_200 > 1.03) == 0:
            # then chains have converged and the mcmc analysis can terminate
            chainlink += 1
            print 'fitM200: chains have converged after' + \
                ' {0} links'.format(chainlink * N_check)
            break
        else:
            # continue with MCMC analysis
            if verbose:
                print 'fitM200: Convergence test results:'
                print 'R_M_200 = ' + \
                    '{0} [Halo_0,... ,halo_(N_halo-1)]'.format(R_M_200)
                print 'We require that R < 1.03 for all halos'
            chainlink += 1
            print 'fitM200: chains have not converged after' + \
                ' {0} links, continuing'.format(chainlink * N_check)
        del_t = (time.time() - t_start) / 60
        print 'fitM200: It took {0:0.0f} minutes'.format(del_t) + \
            'to complete {0} chains of {1} links'.format(N_chain, N_check)

    # Calculate the chi^2 and reduced chi^2 of the best fit model (i.e. mean of
    # the margonalized M_200 chains for each halo)
    N_data, chi2, redchi2 = \
        chi2calc(cat[:, key[ellip_meas[0]]],
                 cat[:, key[ellip_meas[1]]], var_e1, var_e2,
                 del_r, invSigmacr, cos2phi, sin2phi, z_halo,
                 mask, mean_M_200_dist, h_scale, Om, Ol, Or, verbose)

    # trim the unused end values from the del_c and r_s arrays
    M_200 = M_200[:chainlink * N_check, :, :]
    # pickle array
    filename = prefix + '_M_200.pickle'
    F = open(filename, 'w')
    pickle.dump(M_200, F)
    F.close()

    ##########################################################################
    # PRESENT RESULTS
    ##########################################################################
    if verbose:
        'fitM200: finished MCMC analysis, creating output and plots'

    # create the blank confidence limit arrays
    M_200_ul = numpy.zeros(N_halo)
    M_200_ll = numpy.zeros(N_halo)

    # create M_200 array that has been collapsed along the chains dimension
    rows = numpy.shape(M_200)[0] - N_burn
    M_200_flat = numpy.zeros((rows * N_chain, N_halo))
    for h in numpy.arange(N_halo):
        M_200_flat[:, h] = \
            numpy.reshape(M_200[N_burn:, h, :], (rows * N_chain,))
        # calculate the bias corrected percent confidence limits
        M_200_ll[h], M_200_ul[h] = \
            bcpcl(mean_M_200_dist[h], M_200_flat[:, h], 1)

    # print mean and estimated variance of parameters
    print "Results in units of 1e14 solar mass"
    for h in numpy.arange(N_halo):
        if halo_names is not None:
            print 'Halo {0}:'.format(halo_names[h])
        else:
            print 'Halo {0}:'.format(h)
        print 'M_200 = {0:1.2e} '.format(mean_M_200_dist[h]) + \
            '+/- {0:1.2e} \t'.format(numpy.sqrt(V_M_200[h])) + \
            'Gelmans & Rubin Confidence Limits'
        print 'M_200 Mean and Bias Corrected 1 sigma Confidence Limits:'
        print 'Mean = {0:1.2e} '.format(mean_M_200_dist[h]) + \
            'UCL = {0:1.2e} LCL ={1:1.2e}\n\n'.format(M_200_ul[h],
                                                      M_200_ll[h])
    # report the chi^2 of the mean model
    print 'Comparing the mean model with {0}'.format(N_data) + \
        ' source galaxies results in:'
    print 'chi^2 = {0:0.0f}'.format(chi2)
    print 'reduced chi^2 = {0:0.2f}'.format(redchi2)

    outputname = prefix + '_results.txt'
    print 'Results are stored as ', outputname
    R = open(outputname, 'w')
    save_num_output(R, prefix, N_halo, halo_names, mean_M_200_dist,
                    V_M_200, redchi2, M_200_ul, M_200_ll, N_data, chi2)
    R.close()

    # plot some of the results
    #----------------------------------------------------------------------
    color = (
        'b', 'g', 'r', 'c', 'm', 'y', 'k', '0.9', '0.8', '0.7', '0.6',
        '0.5', '0.4', '0.3', '0.2', '0.1')
    if N_halo > 1:
        # Then plot M_200_halo(i) vs M_200_halo(i+1)
        #from math import factorial
        #N_plots = factorial(N_halo) / 2
        for i in range(N_halo):
            for j in range(i + 1, N_halo):
                fig = pylab.figure()
                for c in numpy.arange(N_chain):
                    pylab.loglog(M_200[:N_burn, i, c], M_200[:N_burn, j, c],
                                 c=color[c], alpha=0.3)
                    pylab.loglog(M_200[N_burn:, i, c], M_200[N_burn:, j, c],
                                 c=color[c], ls='-',
                                 label='Chain {0}'.format(c))
                    # pylab.plot(M_200[:N_burn,i,c],M_200[:N_burn,j,c],
                    #            c=color[c],alpha=0.3)
                    #pylab.plot(M_200[N_burn:,i,c],M_200[N_burn:,j,c],
                    #           c=color[c],ls='-',label='Chain {0}'.format(c))
                pylab.legend(loc=0)
                if halo_names is not None:
                    pylab.title(
                        'MCMC Chains for Halo' +
                        ' {0} & {1}'.format(halo_names[i], halo_names[j]))
                    pylab.xlabel('Halo ' + halo_names[i] +
                                 ' $M_{200}$ $(10^{14} M_{\odot})$')
                    pylab.ylabel(
                        'Halo ' + halo_names[j] +
                        ' $M_{200}$ $(10^{14} M_{\odot})$')
                    figname = \
                        prefix + \
                        '_M200_halo{0}_halo{1}'.format(halo_names[i],
                                                       halo_names[j])
                else:
                    pylab.title('MCMC Chains for Halo {0} & {1}'.format(i, j))
                    pylab.xlabel(
                        'Halo ' + str(i) + ' $M_{200}$ $(10^{14} M_{\odot})$')
                    pylab.ylabel(
                        'Halo ' + str(j) + ' $M_{200}$ $(10^{14} M_{\odot})$')
                    figname = prefix + '_M200_halo{0}_halo{1}'.format(i, j)
                pylab.savefig(figname + '.pdf', bbox_inches='tight')

    # plot M_200 distribution
    for h in numpy.arange(N_halo):
        fig = pylab.figure()
        for c in numpy.arange(N_chain):
            bin_array = \
                numpy.logspace(numpy.log10(numpy.min(M_200[N_burn:, h, c])),
                               numpy.log10(numpy.max(M_200[N_burn:, h, c])),
                               bins)
            #pylab.hist(M_200[N_burn:,h,c],bins=bins,ec=color[c],fill=False,
            #label='Chain {0}'.format(c))
            pylab.hist(
                M_200[N_burn:, h, c], bins=bin_array, ec=color[c],
                fill=False, label='Chain {0}'.format(c))
        #pylab.hist(M_200_flat[:, h], bins=bins, ec='k', lw=3, fill=False,
        #label='All Chains')
        bin_array = numpy.logspace(
            numpy.log10(numpy.min(M_200_flat[:, h])),
            numpy.log10(numpy.max(M_200_flat[:, h])), bins)
        pylab.hist(M_200_flat[:, h], bins=bin_array, ec='k', lw=3,
                   fill=False, label='All Chains')
        pylab.legend(loc=0)
        if halo_names is not None:
            pylab.title('Halo {0}'.format(halo_names[h]))
            figname = prefix + '_M_200_halo{0}'.format(halo_names[h])
        else:
            pylab.title('Halo {0}'.format(h))
            figname = prefix + '_M_200_halo{0}'.format(h)
        pylab.xlabel('$M_{200} 10^{14} (M_{\odot})$')
        pylab.ylabel('Number')
        pylab.xscale('log')
        pylab.savefig(figname + '.pdf', bbox_inches='tight')

    # Make contour plots of r_s vs del_c for each halo
    if N_halo > 1:
        # Then plot M_200_halo(i) vs M_200_halo(i+1)
        #from math import factorial
        #N_plots = factorial(N_halo) / 2
        for i in range(N_halo):
            for j in range(i + 1, N_halo):
                fig = pylab.figure()
                # pylab.hexbin(M_200_flat[:,i],M_200_flat[:,j],gridsize=bins)
                pylab.hexbin(M_200_flat[:, i],
                             M_200_flat[:, j],
                             gridsize=bins,
                             xscale='log',
                             yscale='log')
                cb = pylab.colorbar()
                cb.set_label('$N_{steps}$')
                if halo_names is not None:
                    pylab.title(
                        'MCMC Step Density for Halo ' +
                        '{0} & {1}'.format(halo_names[i], halo_names[j]))
                    pylab.xlabel(
                        'Halo ' + halo_names[i] +
                        ' $M_{200}$ $  (10^{14} M_{\odot})$')
                    pylab.ylabel(
                        'Halo ' + halo_names[j] +
                        ' $M_{200}$ $ (10^{14} M_{\odot})$')
                    figname = prefix + \
                        '_hist_M200_halo{0}_halo{1}'.format(
                            halo_names[i],
                            halo_names[j])
                else:
                    pylab.title(
                        'MCMC Step Density for Halo {0} & {1}'.format(i, j))
                    pylab.xlabel(
                        'Halo ' + str(i) + ' $M_{200}$ $ (10^{14} M_{\odot})$')
                    pylab.ylabel(
                        'Halo ' + str(j) + ' $M_{200}$ $ (10^{14} M_{\odot})$')
                    figname = prefix + \
                        '_hist_M200_halo{0}_halo{1}'.format(i, j)
                pylab.savefig(figname + '.pdf', bbox_inches='tight')
    pylab.show()
    return


def trace_plot(c, accept, trace_NW, trace_SE, N_burn, prefix):
    plt.title("trace plot for chain {0} w/ ".format(c) +
              " {0:.0f}% acceptance".format(accept))
    plt.plot(trace_NW, label='NW')
    plt.plot(trace_SE, label='SE')
    #trace_NW.sort()
    plt.axvspan(0, N_burn, color='grey', label='burn in', alpha=0.3)
    #plt.vlines(N_burn, trace_NW[0], trace_NW[len(trace_NW) - 1])
    plt.legend(loc='best')
    plt.savefig(
        prefix + '_trace_{0}.png'.format(c), bbox_inches='tight')
    plt.close()


def save_num_output(F, prefix, N_halo, halo_names, mean_M_200_dist,
                    V_M_200, redchi2, M_200_ul, M_200_ll, N_data, chi2):
    F.write('These results were generated by nfwMCMC.fitM200()\n' +
            'M_200 units are in 1e14 *h^{-1} solar mass.\n\n')
    for h in numpy.arange(N_halo):
        if halo_names is not None:
            F.write('Halo {0}:\n'.format(halo_names[h]))
        else:
            F.write('Halo {0}:\n'.format(h))
        F.write(
            'M_200 = {0:1.2e} +/- '.format(mean_M_200_dist[h]) +
            '{0:1.2e} \t'.format(numpy.sqrt(V_M_200[h])) +
            ' Gelmans & Rubin Confidence Limits\n')
        F.write('M_200 Mean and Bias Corrected 1 sigma Confidence Limits:\n')
        F.write(
            'Mean = {0:1.2e} '.format(mean_M_200_dist[h]) +
            'UCL = {0:1.2e} '.format(M_200_ul[h]) +
            'LCL = {0:1.2e}\n\n'.format(M_200_ll[h]))
    F.write(
        'Comparing the mean model with {0} '.format(N_data) +
        'source galaxies results in:\n')
    F.write('chi^2 = {0:0.0f}\n'.format(chi2))
    F.write('reduced chi^2 = {0:0.2f}'.format(redchi2))
    F.close()


def input_checks(N_check, N_burn, halos, halo_names, N_chain, parambounds,
                 seedbounds, N_halo):
    if N_check <= N_burn:
        raise ValueError('fitM200: N_check must be greater than N_burn,' +
                         'exiting')

    if numpy.size(halos) % 3. != 0:
        raise ValueError('fitM200: an invalid "halos" input was ' +
                         'specified, exiting')
    if halo_names is not None:
        if numpy.size(halo_names) != N_halo:
            raise ValueError('the list size of halo_names must be equal ' +
                             ' to the number of halos specified in ' +
                             ' halos, exiting')


def print_little_success(i, k, x, logL_new, logL_old):
    print 'mcmcengine @ step {0}:\n over'.format(i) +\
        ' {0} N_iterations without a successful step'.format(k)
    print 'log probability of making bad step = {0}'.format(x)
    print 'logL_new - logL_old = ' + \
        '{0}'.format(logL_new - logL_old)


def print_stuck(i, M_200_new, parambounds):
    print 'nfwMCMC.mcmcengine(): N_iteration ', i
    print 'has been stuck for 10000 draws\n'
    print 'M_200_new = {0} \n '.format(M_200_new)
    print 'parambounds are {0}'.format(parambounds)


def adjust_accept_rate(N_iter, cpropdist, accept, N_burn, N_check, chainset):
    # check for crazy values
    if cpropdist == np.inf:
        raise ValueError("cpropdist = inf!")
    if cpropdist == np.nan:
        raise ValueError("cpropdist = nan!")

    # only adjust acceptance rate for first set of chains
    # only allow propdist to be adjusted before burn-in
    if chainset == 0 and N_iter <= N_burn:
        if accept > 0.6 and cpropdist < 1e1:
            cpropdist = increase_propdist(cpropdist)
            print 'increasing propdist to {0}'.format(cpropdist)
        if accept < 0.2 and cpropdist > 1e-2:
            cpropdist = decrease_propdist(cpropdist)
            print 'decreasing propdist to {0}'.format(cpropdist)

    return cpropdist


def increase_propdist(cpropdist):
    if cpropdist > 1.0:
        return cpropdist ** 1.3
    elif 1.0 - cpropdist < 1e-2:
        return 1.01
    else:
        return cpropdist ** 0.7


def decrease_propdist(cpropdist):
    if cpropdist > 1.0:
        return cpropdist ** 0.7
    elif cpropdist - 1.0 < 1e-2:
        return 0.99
    else:
        return cpropdist ** 1.3


def outsideBound(index, M_200_new, parambounds):
    """
    check if M_200_new(index) is outside parambound
    """
    return np.logical_or(M_200_new[index] < parambounds[0],
                         M_200_new[index] > parambounds[1])

# ------------- functions to be tested below this line  ---------------------

def check_seedbound_against_parambounds(parambounds, seedbounds, i):
    if parambounds[0] > seedbounds[2 * i]:
        raise ValueError(
            'seedbound of chain ' +
            '{0} : {1} <  upper '.format(i, seedbounds[2 * i]) +
            'parabounds : {0}'.format(parambounds[0]))
    if parambounds[1] < seedbounds[2 * i + 1]:
        raise ValueError(
            'seedbound of chain ' +
            '{0} = {1} >  upper '.format(i, seedbounds[2 * i + 1]) +
            'parabounds = {0}'.format(parambounds[1]))
    if seedbounds[2 * i] >= seedbounds[2 * i + 1]:
        raise ValueError('lower seedbound >= upper seedbound for ' +
                         'chain {0}'.format(i))
    return


def calc_stuff(N_halo, N_gal, cat, h_scale, Om, Ol, halos,
               key, coord, verbose=True):
    """
    calculate del_a and del_d and invSigmacr
    Args :
    del_a = np.array - unit : Mpc
    del_d = np.array - unit : Mpc
    invSigmacr = np.array
    """
    # initialize a bunch of arrays for use
    #var_invSigmacr = numpy.zeros((N_gal,N_halo))
    invSigmacr = numpy.zeros((N_gal, N_halo))
    del_a = numpy.zeros((N_gal, N_halo))
    del_d = numpy.zeros((N_gal, N_halo))

    for h in numpy.arange(N_halo):
        if verbose:
            print 'fitM200: evaluating galaxies with respect to halo' + \
                '{0}'.format(h)
        for g in numpy.arange(N_gal):
            # Calculate the angular separations
            del_a[g, h], del_d[g, h] = \
                tools.angcomp(cat[g, key[coord[0]]], cat[g, key[coord[1]]],
                              halos[h * 3], halos[h * 3 + 1])

            # *# Note that I could speed up the code by masking galaxies
            # outside
            # *# the apodizing radii and not calculating their invSigmacr
            # Calulate inverse sigma critical
            if halos[h * 3 + 2] < cat[g, key[coord[2]]]:
                invSigmacr[g, h] = \
                    1 / cosmo.lensgeo(
                        halos[h * 3 + 2], cat[g, key[coord[2]]],
                        h_scale, Om, Ol)['sigcr']

            # Calculate the standard error on sigma critcal
            ###################################################
            # Need to add this section
            ###################################################
        # convert deg anglular separation to Mpc separation
        del_a[:, h] *= 60 * \
            cosmo.ProjectedLength(halos[h * 3 + 2], h_scale, Om, Ol)
        del_d[:, h] *= 60 * \
            cosmo.ProjectedLength(halos[h * 3 + 2], h_scale, Om, Ol)
    return invSigmacr, del_a, del_d


def mask_galaxies(cat, N_gal, N_halo, del_r, r_bounds, halos, key, coord):
    """
    Mask galaxies that are outside the apodizing radii or foreground of
    the halo
    Stability - untested - should give 2386 galaxies for James El Gordo
    data
    """
    mask_r_inner = numpy.zeros((N_gal, N_halo)) != 0
    mask_r_outer = numpy.zeros((N_gal, N_halo)) != 0
    mask_z = numpy.zeros((N_gal, N_halo)) != 0

    for h in numpy.arange(N_halo):
        mask_r_inner[:, h] = del_r[:, h] > r_bounds[0]

        # mask out galaxies that are not within any outer radii of a halo
        mask_r_outer[:, h] = del_r[:, h] < r_bounds[1]

        # needs background galaxies, not foreground
        mask_z[:, h] = cat[:, key[coord[2]]] > halos[h * 3 + 2]


    # mask out all the galaxies that are within any of the inner radii of
    # either halo
    mask_r_inner = numpy.sum(mask_r_inner, axis=1) == N_halo

    #mask = mask_r*mask_z
    mask = \
        numpy.reshape(mask_r_inner, (numpy.size(mask_r_inner), 1)) * \
        mask_r_outer * mask_z

    for h in numpy.arange(N_halo):
        print(
        'fitM200: There are {0} '.format(numpy.sum(mask[:, h])) +
        'background galaxies for Halo {0}'.format(h))

    return mask


def compute_new_centroid(centroid_offset, centroid_old):
    """to be tested
    centroid_offset = np.array with shape = (N_halo, 2)
    centroid_old = assumes the same structure as the "halos" parameter
    which looks like [ra1, dec1, z1, ..., ra_n, dec_n, z_n]
    """
    new_centroid = [[centroid_offset[i, 0] + centroid_old[3 * i],
                     centroid_offset[i, 1] + centroid_old[3 * i + 1],
                     centroid_old[3 * i + 2]]
                     for i in range(centroid_offset.shape[0])]

    # this straightens the nested list of lists as one list
    new_centroid = sum(new_centroid, [])

    return new_centroid


def ElGo_centroid_outside_bd_or_if_centroids_exchanged(
        centroid_offset, centroid_step, centroid_old, centroid_new,
        verbose=True):
    """ check if we drew values that tells us that the centroid is outside
    a believable range - also stops centroids from being exchanged

    This is customized for El Gordo with NW and SE cluster
    You will need to modify this function if you are doing analysis for
    another cluster

    Return
    bad_step = number, > 0 for a bad step
    """
    # loop through the del_a and del_d to make sure
    # np.sqrt(del_a ** 2 + del_d ** 2) < 15 * 2D centroid step
    outside_bound = \
        [np.sqrt(centroid_offset[i, 0] ** 2 + centroid_offset[i, 1] ** 2) >
         15 * np.sqrt(2 * centroid_step ** 2)
         for i in range(centroid_offset.shape[0])]

    if verbose:
        print "custom centroid location check for El Gordo "

    # note that the inequality signs have to be changed according to
    # the relative location of the subclusters
    # RA_NW should be < RA_SE RA increases in the east direction
    # DEC_NW should be > DEC_SE
    centroid_position_exchanged = \
        centroid_new[0] > centroid_new[3 * 1] and\
        centroid_new[1] < centroid_new[3 * 1 + 1]

    # print "outside bound = {0}".format(outside_bound)
    # print "centroid_position_exchanged = {0}".format(
    #     centroid_position_exchanged)
    bad_step = np.sum(outside_bound) + centroid_position_exchanged

    return bad_step


def compute_rotation(cd, del_a, del_d):
    """
    Calculate the angle of the galaxy CCW from +x axis of the halo
    use trig identities while preserving magnitude
    sinphi = del_y/del_r
    cos2phi = 1- 2sinphi**2
    sin2phi = 2sinphi*cosphi
    """

    cd_inv = numpy.linalg.inv(numpy.array(cd))
    cd_mag = numpy.sqrt(numpy.sum(cd_inv[0, :] ** 2))
    del_x = cd_inv[0, 0] / cd_mag * del_a + cd_inv[0, 1] / cd_mag * del_d
    del_y = cd_inv[1, 0] / cd_mag * del_a + cd_inv[1, 1] / cd_mag * del_d

    del_r = np.sqrt(del_x ** 2 + del_y ** 2)
    cos2phi = 1. - 2. * (del_y / del_r) ** 2
    sin2phi = 2. * del_x * del_y / (del_r ** 2)
    return cos2phi, sin2phi


def draw_new_mass_for_halo_h(h, M_200_new, stuckTimes, M_200_old,
                             cpropdist, N_iter, parambounds, verbose=False):
    M_200_new[h] = \
        np.random.randn() * cpropdist + M_200_old[h].copy()
    stuckTimes[h] += 1

    # when not it's not the first iteration of any loops
    # and it's been stuck for multiples of 10000 times
    # print warning
    if stuckTimes[h] != 0 and \
            stuckTimes[h] % 10000 == 0 and\
            verbose and N_iter != 0:
        print_stuck(N_iter, M_200_new[h], parambounds)

    return M_200_new[h], stuckTimes[h]


def get_new_centroids(N_halo, centroid_old, centroid_step, verbose=False):
    """
    N_halo = int
    centroid_old = assumes the same structure as the "halos" parameter
    which looks like [ra1, dec1, z1, ..., ra_n, dec_n, z_n]
    centroid_step = float, in arcseconds
    """
    centroid_offset = \
        np.array([numpy.random.randn(2) * centroid_step
                  for k in range(N_halo)])
    centroid_new = \
        compute_new_centroid(centroid_offset, centroid_old)

    while(ElGo_centroid_outside_bd_or_if_centroids_exchanged(
        centroid_offset, centroid_step, centroid_old,
            centroid_new, verbose)):

        centroid_offset = \
            np.array([numpy.random.randn(2) * centroid_step
                     for k in range(N_halo)])
        centroid_new = \
            compute_new_centroid(centroid_offset, centroid_old)
    return centroid_new


# def nfwparam(M_200,z,h_scale=0.7,Om=0.3,Ol=0.7,Or=0.0):
#    '''
#    Inputs:
#    M_200 = [array of floats; units= M_sun]
#    Outputs:
#    del_c, r_s = characteristic overdensity of the CDM halo, scale radius of
#    the halo (Mpc)
#    '''
#    M_200 *= 1e14
# conversions
# kginMsun = 1.98892*10**30 #kg in a solar mass
# minMpc = 3.08568025*10**22 # m in a Megaparsec
#
#    rho_cr = cosmo.rhoCrit(z,h_scale,Om,Ol,Or)/kginMsun*minMpc**3
# calculate the r_200 radius
#    r_200 = (M_200*3/(4*numpy.pi*200*rho_cr))**(1/3.)
# calculate the concentration parameter based on Duffy et al. 2008
#    c = 5.71/(1+z)**0.47*(M_200*h_scale/2e12)**(-0.084)
#    del_c = 200/3.*c**3/(numpy.log(1+c)-c/(1+c))
#    r_s = r_200/c
#    return del_c, r_s


fitM200.__doc__ = \
    '''
    Simultaniously fits multiple NFW profiles to the Weak Lensing data of in
    input catalog of galaxy shapes (e1, e2) and redshifts (z).  This is very
    similar to lensMCMC_v0.0.py, except that instead of fitting for both rho_s
    and r_s we are only fitting for r_s, essentially. Techically this program
    fits for M_200 of each halo, which is realted to rho_s and r_s through
    M_200/c emperical scaling relations (e.g. Duffy et al. 2008) and general
    properties of the NFW profile (see e.g. Springer's text).  More details are
    contained in my OneNote notebook Programs>Miscelanious Programs>nfwMCMC.

    Input:
    halos = [(list of floats length in multiples of 3); units:(degrees,
            degrees,none)] a list of floats providing the (ra, dec, redshift)
            of each halo. There is no limit to the number of halos, e.g.:
            (ra_0,dec_0,z_0,ra_1,dec_1,z_1,...ra_n,dec_n,z_n)
    catalog = [string] name of the ttype'd catalog containing all the galaxy
              data
    coord = [(string, string, string)] list of ttype names for the galaxy (ra,
            dec, redshift) columns in catalog. Note that the ra and dec values
            of the catalog should be in units of Degrees.
    ellip_meas = [(string,string,string,string)] list of ttype names for the
                 galaxy (e1, e2, e1_error, e2_error). Ellipticities should be
                 defined with respect to the standard x,y coordinate system and
                 it is assumed that ra and dec axes are orientated such that
                 +x=-RA and +y=+Dec unless a cd matrix is input
    ellip_intstd = [float] standard deviation of intrinsic ellipticities
    r_bounds = [(float,float); units:(Mpc,Mpc)] the inner and outer radial
               bounds of galaxies to be used in the lensing analysis of each
               halo. No galaxies within the inner bound of any halo will be
               used, while no galaxy outside of all the halo outer bounds will
               be used.
    parambounds = [(float,float); units:(1e14 h^{-1} M_sun,
                   1e14 h^{-1} M_sun)] defines the bounds
                  (M_200_min, M_200_max) of the uniform prior to apply to each
                  M_halo parameter
    N_maxlinks = [integer] program will stop after this number of chain link
                 attempts even if the chains have not converged
    N_chain = [integer] Number of chains to run in parallel and to check
              convergence against
    N_burn = [integer] Number of initial chain links to discard
    N_check = [integer] Check convergence at multiples of this interval of
              successful links. Note that it must be greater than N_burn
    propdist = [float; units:M_sun] Program assumes a normal proposal
               distribution with sigma_M_200 = propdist and mean=0
    seedbounds = defined bounds of uniform distribution from which to randomly
                 draw the starting point of each chain, defaults to None in
                 which case the parambounds will be used but allows each halos
                 starting point bounds to be specified. Note that if one halos
                 seed bounds are specified then all must be
                 (M_200_min_0,M_200_max_0,...
                 ,M_200_min_(N_halo-1),M_200_max_(N_halo-1))
    bins = [integer] number of bins to use when presenting the histograms of
           step distributions. Defaults to 20
    halo_names = [(list of strings size N_halo)] This is an optional input that
                 is used to format the output figures and results.
    cd = [((float,float),(float,float))] This is the CD matrix typically found
         in WCS header definitions which provides the scale and orientation of
         the RA,Dec coordinate system with the x,y pixel coordinate system such
         that [[x];[y]] = [[CD_11, CD_12];[CD_21, CD_22]] [[RA]; [Dec]]. For
         the purposes of this program the pixel scale is not important as we
         are only after the respective orientation of the two coordinates.
    h_scale = 0.7 Hubble parameter scale H = h*100km/s/Mpc
    Om = 0.3 Matter enegry density
    Ol = 0.7 Dark Energy energy density
    Or = 0.0 Radiation energy density
    verbose = True
    centroid_step = this is the one dimensional step size that you want the
        centroid to move, e.g. true 2D stepsize is sqrt(2 * centroid_step ** 2 )

    Output:
    '''
