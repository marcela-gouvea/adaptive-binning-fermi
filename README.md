# Adaptive Binning for Fermi-LAT Analysis

This repository presents an Adaptive-Binning method to be used for Fermi-LAT data analysis using the fermipy lib. The method allows the creation of more significant light curves. Each bin size is adjusted by requiring a minimum significance ```TS >= 𝜹₀``` in which TS is the test statistics and 𝜹₀ is what we call acceptance criteria. The minimun bin size is established by a variable we call ```unc``` 

Deafults: 

```𝜹₀ = 25```

```unc = 1```

Remember to change these according to your purposes

## Requisites:

Numpy - https://numpy.org/

Fermipy - https://fermipy.readthedocs.io/en/latest/

Astropy - https://www.astropy.org/

You can use the same environment used for Fermipy analysis

## Usage

To use our codes, you should know:
- First, download our codes and create a folder for the analysis;
- Change ```𝜹₀``` and ```unc``` according to your criteria;
- Before running the code, you must create a "files.txt" file which contains the whole path to the photon files. Also, create a configuration file for the analysis and a copy (one named "config.yaml" and "config1.yaml"). You must put the whole path to the spacecraft file in the configuration;
- Run the "adaptive.py" file. You can choose the minium uncertainty for each time bin by changin the 'unc' value in the beggining of the code; 
- The final lightcurve generated along with associated spectral indexes can be found at "total-lightcurve.txt", in which each column represent, respectively: **time (MJD), time error (MJD), flux (MeV / cm² s), flux error (MeV / cm² s), binsize, alpha, alpha_error, beta, beta_error;**

## Questions and Contact

If you have any questions regarding our codes or the output results, please contact Marcela Gouvêa at marcela.gouvearrb@gmail.com
