"""
Read in h5 files containing inputs of MCMAC and output a table for paper
"""

from __future__ import division
import sys
import os
import h5py
sys.path.append(os.path.abspath("../multiDim_plot"))
from plotmod import round_to_n_sig_fig
from astropy.stats import biweight_midvariance as S_BI
from astropy.stats import biweight_location as C_BI

# import homebrewed module
import pandas as pd

if len(sys.argv) < 2:
    sys.exit("ERROR " +
             "you have to specify which hdf5 file you want as argument")
else:
    h5_file = sys.argv[1]

f = h5py.File(h5_file, "r")
m_main = f["m_main"][...] / (1e14)
m_sub = f["m_sub"][...] / (1e14)
z_main = f["z_main"][...]
z_sub = f["z_sub"][...]
d_proj = f["D_proj"][...]
f.close()


dat = pd.read_table("../call_MCMAC/el_gordo_chain.txt", sep="\s*", skiprows=3)
c_main = dat["cNW"]
c_sub = dat["cSE"]
# convert the mass into concentration
# estimate the mean and 1 sigma confidence level of the inputs

# we want to be able to output these filling the gaps to a file
# this would be incorporated as part of ElGordo1.tex
paper1 = "/Users/karenyng/Documents/Research/papers/ElGo_paper/texFiles/"
F = open(paper1 + "input_table.tex", "w")
F.write("\setcounter{table}{0}\n")
F.write("\\begin{table}\n")
F.write(
    "\caption{Properties of the sampling PDFs of the Monte Carlo " +
    "simulation} \n")
F.write("\\begin{center} \n")
F.write("\\begin{tabular}{@{}lcccc}\n")
F.write("\hline \hline Data & Units & Location & Scale & Ref\\\\ \hline\n")
F.write(
    "$M_{200c_{\mathrm{NW}}}$ & $10^{14} h_{70}^{-1}$ M$_{\odot}$" +
    " &{0:.1f}&{1:.1f}& \citetalias{{Jee13}}\\\\ \n".format(
        round_to_n_sig_fig(C_BI(m_main), 2), round(S_BI(m_main), 2)))
F.write("c$_{{\mathrm{{NW}}}}$ &  & {0:.2f}& {1:.2f}& ".format(
    round_to_n_sig_fig(C_BI(c_main), 2), round(S_BI(c_main), 2)) +
    "\citetalias{{Jee13}} \\\\ \n")
F.write("$M_{200c_{\mathrm{SE}}}$ & $10^{14} h_{70}^{-1}$ M$_{\odot}$" +
        " &{0:.1f}&{1:.1f} & \citetalias{{Jee13}}\\\\ \n".format(
            round_to_n_sig_fig(C_BI(m_sub), 2), round(S_BI(m_sub), 2)))
F.write("$c_{{\mathrm{{SE}}}}$ &  & {0:.2f} & {1:.2f}& ".format(
    round_to_n_sig_fig(C_BI(c_sub), 2), round(S_BI(c_sub), 2)) +
        "\citetalias{Jee13}\\\\ \n")
F.write("$z_{{\mathrm{{NW}}}}$ &  & {0:.5f}".format(
        round_to_n_sig_fig(C_BI(z_main), 5)) +
        " & {0:.5f}& \citetalias{{M12}}, \citetalias{{Sifon13}}\\\\ \n".format(
        round_to_n_sig_fig(S_BI(z_main), 5)))
F.write("$z_{{\mathrm{{SE}}}}$ &  & {0:.5f} & ".format(C_BI(z_sub)) +
        "{0:.5f}& \citetalias{{M12}}, \citetalias{{Sifon13}}\\\\ \n".format(
        round_to_n_sig_fig(S_BI(z_sub), 5)))
F.write("d$_{{\mathrm{{proj}}}}$ & Mpc & {0:.2g} &".format(
    round_to_n_sig_fig(C_BI(d_proj), 2)) +
        "{0:.2g} & \citetalias{{Jee13}} \\\\ \n".format(
            round_to_n_sig_fig(S_BI(d_proj), 2)))
F.write("\hline \n")
F.write("\end{tabular} \n")
F.write("\end{center} \n")
F.write("\label{tab:inputs} \n")
F.write("\end{table} \n")
F.close()
