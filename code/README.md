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

## folders 
* Fig. 1 - in `Position/plot_Fig1.ipynb` 
* fits_image - contains fits files for configuration plots 
* sensitivity_analysis 
	* contains sensitivity analysis 
	* the calculation of the probability of different merger scenario 
* multiDim_plot 
	* contains the main plots in appendix B of the paper 
* Polarization 
	* contains calculation and derivation of polarization weights 
