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
- Before running the code, you must create a "files.txt" file which contains the whole path to the photon files. Also, create a configuration file for the analysis and a copy (one named "config.yaml" and "config1.yaml").
- Your config file must be similar to the following exemple in order for the code to work. You must put the whole path to the spacecraft file in the configuration. Also, make sure that the ```tmin```, ```tmax``` and ```target``` are on lines 11, 12 and 14, respectively. We are still working to improve the config reading. The code does not support double selection.
  
 ```
data:
  evfile: files.txt
  scfile: L23102221342953B6354C62_SC00.fits
binning:
  roiwidth: 15.0
  binsz: 0.1
  binsperdec: 3
selection:
  emin: 100
  emax: 1500000
  zmax: 90
  evtype: 3
  tmin: 694224005
  tmax: 718588805
  filter: null
  target: 3C454.3
gtlike:
  edisp: true
  irfs: P8R3_SOURCE_V3
  edisp_disable:
  - isodiff
  - galdiff
model:
  src_roiwidth: 15.0
  galdiff: gll_iem_v07.fits
  isodiff: iso_P8R3_SOURCE_V3_v1.txt
  catalogs:
  - 4FGL
```

- Run the "adaptive.py" file. You can choose the minium uncertainty for each time bin by changin the 'unc' value in the beggining of the code; 
- The final lightcurve generated along with associated spectral indexes can be found at "total-lightcurve.txt", in which each column represent, respectively: **time (MJD), time error (MJD), flux (MeV / cm² s), flux error (MeV / cm² s), binsize, alpha, alpha_error, beta, beta_error;**

## Questions and Contact

If you have any questions regarding our codes or the output results, please contact Marcela Gouvêa at marcela.gouvearrb@gmail.com
