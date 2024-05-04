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
- Run the "adaptive.py" file. You can choose the minium uncertainty for each time bin by changin the 'unc' value in the beggining of the code; 
- The final lightcurve generated along with associated spectral indexes can be found at "total-lightcurve.txt", in which each column represent, respectively: time, time error, flux, flux error, binsize, alpha, alpha_error, beta, beta_error;

## Questions and Contact

If you have any questions regarding our codes or the output results, please contact Marcela Gouvêa at marcela.gouvearrb@gmail.com
