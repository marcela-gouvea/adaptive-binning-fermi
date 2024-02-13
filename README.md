# Adaptive Binning for Fermi-LAT Analysis

This repository presents an Adaptive-Binning method to be used for Fermi-LAT data analysis using the fermipy lib

Requisites:

Numpy - https://numpy.org/

Matplotlib - https://matplotlib.org/

Fermipy - https://fermipy.readthedocs.io/en/latest/

Astropy - https://www.astropy.org/

To use our codes, you have to:
1- Create a "files.txt" file which contains the path to the photon files
2- Create a configuration file for the analysis and a copy (one named "config.yaml" and "config1.yaml")

First, run the "adaptive.py" file and then run the "fix-spectral.py" file. 
The final lightcurve generated can be found at "fixed-bins.txt", in which each column represent, respectively: time, time error, flux, flux error, binsize
The temporal evolution of the spectral index can be found at "spectral-index.txt", in which each column represent, respectively: alpha, alpha error, beta, beta error

If you have any questions regarding our codes or the output results, please contact Marcela Gouvêa at marcela.gouvearrb@gmail.com
