'''
Author : Karen Yin-Yee Ng, based on work by Will Dawson
Comment:
This is the code that gets the fits file from James for the ElGordo weak
lensing mass analysis and massage it so that it can be digested by Will 's
MCMAC.py code parallel version!

Usage: examine EVERY VARIABLE assignment in this file !
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import cPickle
import pandas as pd
from astropy.coordinates import ICRS
import astropy.coordinates.angle_utilities as ang_util
from astropy.cosmology import FlatLambdaCDM
import h5py
import profiles
from MCMAC import vrad
from scipy.stats import norm
import os

# home brewed modules
#from MCMAC import MCengine
#from nfwMCMC_modularized import mask_galaxies, calc_stuff, \
#    compute_new_centroid,\
#    ElGo_centroid_outside_bd_or_if_centroids_exchanged,\
from nfwMCMC_modularized import get_new_centroids

# inputs for the TSM module
# N_sample is the number of iterations of the MC runs
compress_fmt = "gzip"  # default = None
debug = True

fitsfile = "./mcmc.fits"
fits_wcs = "./header.fits"

out_prefix = "PtOneMnWBSG"
print "making directory with name " + out_prefix
os.system("mkdir ./" + out_prefix)
out_prefix = "./" + out_prefix + "/" + out_prefix
h5_outfile = out_prefix + "_inputs.h5"

###
# Cluster Mass distribution arrays
###
# have to look at where the burn-in is and chop it off
# chopping off the first 250 steps of each chain
# the cluster has different upper and lower error limits
# call function from profile.py
# m_sub = profiles.nfwM200(cSE,A200,B200,C200,zSE) #this should be SE
# m_main = profiles.nfwM200(cNW,A200,B200,C200,zNW)  #this should be NW

dat = pd.read_table("./el_gordo_chain.txt", sep="\s*", skiprows=3)
dat["m_NW"] = \
    profiles.nfwM200(dat["cNW"], A200=5.71, B200=-0.084, C200=-0.47, z=0.87)
dat["m_SE"] = \
    profiles.nfwM200(dat["cSE"], A200=5.71, B200=-0.084, C200=-0.47, z=0.87)
m_main = dat["m_NW"]
m_sub = dat["m_SE"]


###
# Cluster Redshift distribution arrays
###
# use M13 Table 1 number
#zNW = (0.86849, 0.00020)
#zSE = (0.87175, 0.00024)

# my bootstrapped estimates with 51 vs 36 members excluding z > 0.886
#zNW = (0.86850, 0.00019378)  # cPickle.load(open("NW_specz.pickle"))
#zSE = (0.87190, 0.00023488)  # cPickle.load(open("SE_specz.pickle"))

zNW = cPickle.load(open("NW_bootstrapRedshift.pickle"))
zSE = cPickle.load(open("SE_bootstrapRedshift.pickle"))
# main cluster is named to be the NW cluster
z_main = zNW
# sub cluster is named to be the SE cluster
z_sub = zSE


###
# find D_proj
###
# draw the same number of samples as the # of mass values
no = dat["m_NW"].shape[0]
# these values r from El Gordo WL paper published values
NW_cent = ICRS('01h02m50.601s -49d15m04.48s')
SE_cent = ICRS('01h02m56.312s -49d16m23.15s')

# assume the location of the centroids to be the published value
# draw values from 4 normal distribution
centroid_old = [[NW_cent.ra.deg, NW_cent.dec.deg, 0.87,
                SE_cent.ra.deg, SE_cent.dec.deg, 0.87]] * no

# hopefully that the arrays shape will be broadcast
N_halos = [2] * no
centroid_steps = [1. / np.sqrt(2) / 60. / 60.] * no
cent_new = map(get_new_centroids, N_halos, centroid_old, centroid_steps)
cent_new = np.array(cent_new)
[NW_ra, NW_dec, NW_z, SE_ra, SE_dec, SE_z] = cent_new.transpose()

# inputs of angular separation has to be in units of radians
ang_sep = \
    ang_util.angular_separation(NW_ra / 180. * np.pi,
                                NW_dec / 180. * np.pi,
                                SE_ra / 180. * np.pi,
                                SE_dec / 180. * np.pi)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
D_proj = cosmo.angular_diameter_distance(0.87) * ang_sep


###
# Debug input
###

def plot_dproj():
    plt.title('D-proj')
    plt.xlabel('Mpc')
    plt.hist(D_proj, bins=100, normed=True, histtype='step')
    plt.savefig(out_prefix + '_D_proj.png')
    plt.close()


def plot_hist_m_main():
    plt.title('m_main')
    plt.xlabel('$10^{15} M_{sun}$')
    plt.hist(m_main / 1.0e15, bins=100, normed=True, histtype='step')
    plt.savefig(out_prefix + '_m_main.png')
    plt.close()


def plot_hist_m_sub():
    plt.title('m_sub')
    plt.xlabel('$10^{15} M_{sun}$')
    plt.hist(m_sub / 1.0e15, bins=100, normed=True, histtype='step')
    plt.savefig(out_prefix + '_m_sub.png')
    plt.close()


if debug is True:
    plot_dproj()
    plot_hist_m_main()
    plot_hist_m_sub()

    if type(z_main) is tuple or type(z_main) is list:
        assert len(z_main) == 2 and len(z_sub) == 2,\
            "check z_main and z_sub type use numpy array if they are not\
            Gaussian"
        z_m = norm.rvs(loc=z_main[0], scale=z_main[1], size=N_sample)
        z_s = norm.rvs(loc=z_sub[0], scale=z_sub[1], size=N_sample)
        v_rad_diff = vrad(z_m, z_s)

        plt.title('z_main')
        plt.xlabel('$redshift$')
        plt.hist(z_m, bins=100, normed=True, histtype='step')
        plt.savefig(out_prefix + '_z_main.png')
        plt.close()

        plt.title('z_sub')
        plt.xlabel('$redshift$')
        plt.hist(z_s, bins=100, normed=True, histtype='step')
        plt.savefig(out_prefix + '_z_sub.png')
        plt.close()
    else:
        v_rad_diff = vrad(z_main, z_sub)
        plt.title('z_main')
        plt.xlabel('$redshift$')
        plt.hist(z_main, bins=100, normed=True, histtype='step')
        plt.savefig(out_prefix + '_z_main.png')
        plt.close()

        plt.title('z_sub')
        plt.xlabel('$redshift$')
        plt.hist(z_sub, bins=100, normed=True, histtype='step')
        plt.savefig(out_prefix + '_z_sub.png')
        plt.close()
    plt.title('v_rad diff')
    plt.hist(v_rad_diff, bins=100, normed=True, histtype='step')
    plt.savefig(out_prefix + '_v_rad_diff.png')
    plt.close()

# try outputting the arrays in compressed form

seed = np.random.get_state()

print "outputting dataset to {0}".format(h5_outfile)
f = h5py.File(h5_outfile, "w")
f.create_dataset("m_main", data=m_main, compression=compress_fmt)
f.create_dataset("m_sub", data=m_sub, compression=compress_fmt)
f.create_dataset("z_main", data=z_main, compression=compress_fmt)
f.create_dataset("z_sub", data=z_sub, compression=compress_fmt)
f.create_dataset("D_proj", data=D_proj, compression=compress_fmt)
f.create_dataset("NW_ra", data=NW_ra, compression=compress_fmt)
f.create_dataset("NW_dec", data=NW_dec, compression=compress_fmt)
f.create_dataset("SE_ra", data=SE_ra, compression=compress_fmt)
f.create_dataset("SE_dec", data=SE_dec, compression=compress_fmt)
f.create_dataset("seed[1]", data=seed[1], compression=compress_fmt)

f["m_main"].attrs["info"] = "May 2014 chains from James"
f["z_main"].attrs["info"] = "July 2014 bootstrap estimates assuming \
    with Will 's CAT zVdisp code"
f["D_proj"].attrs["info"] = \
    "generated by two gaussians with around published centroids\
    with one arcsecs as standard dev"
f["seed[1]"].attrs["info"] = \
    "use >>> np.random.set_state(('MT19937', seed[1], 624, 0, 0.0)) to \
    recover the random state"
f.close()
