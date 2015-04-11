# Files
This contains most of the important codes for the analysis for the El Gordo
paper

## Dependencies 
See python_packages.txt for all the dependencies.
For setting up a virtual environment for computing use

```
$ source elgoenv/bin/activate
$ pip install -r python_packages.txt 
```

Or 
```
$ pip install -I $package_name==$package_version 
```

# To view specific IPython notebook without installation, type:
`http://nbviewer.ipython.org/url/LINK_TO_RAW_IPYTHON_NOTEBOOK_IN_THIS_REPO` in
your browser and you should modified the `LINK_TO_RAW_IPYTHON_NOTEBOOK_IN_THIS_REPO`

## folders 
* Position - Fig. 1 - in `plot_Fig1.ipynb` 
* sensitivity_analysis 
	* contains sensitivity analysis  
	* the calculation of the probability of different merger scenario 
    * in `different scenario probability - polar prior E.ipynb`
    * in `different scenario probability - polar prior NW.ipynb`
  
* multiDim_plot 
	* contains the main plots in appendix B of the paper 
* Polarization 
	* contains calculation and derivation of polarization weights 
