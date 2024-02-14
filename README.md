# Adaptive Binning for Fermi-LAT Analysis

This repository presents an Adaptive-Binning method to be used for Fermi-LAT data analysis using the fermipy lib

## Requisites:

Numpy - https://numpy.org/

Matplotlib - https://matplotlib.org/

Fermipy - https://fermipy.readthedocs.io/en/latest/

Astropy - https://www.astropy.org/

## Usage

To use our codes, you should know:
- First, download our codes and create a folder for the analysis;
- Before running the code, you must create a "files.txt" file which contains the whole path to the photon files. Also, create a configuration file for the analysis and a copy (one named "config.yaml" and "config1.yaml"). You must put the whole path to the spacecraft file in the configuration;
- Run the "adaptive.py" file and then run the "fix-spectral.py" file; 
- The final lightcurve generated can be found at "fixed-bins.txt", in which each column represent, respectively: time, time error, flux, flux error, binsize;
- The temporal evolution of the spectral index can be found at "spectral-index.txt", in which each column represent, respectively: alpha, alpha error, beta, beta error.
- If you want to plot the lightcurve and the spectral index temporal evolution considering a simple power law, you must run "plot.py".

## Questions and Contact

If you have any questions regarding our codes or the output results, please contact Marcela Gouvêa at marcela.gouvearrb@gmail.com
