import os
import numpy as np
import subprocess
from astropy.table import Table
from fermipy.gtanalysis import GTAnalysis

# Constants
MJDREFF = 51910
MJDREFFI = 7.428703703703703e-4

def met_to_mjd(time):
    """Convert mission elapsed time to mean julian date."""
    return MJDREFF + MJDREFFI + np.array(time) / 86400.

def mjd_to_met(time):
    """Convert mean julian date to mission elapsed time."""
    return (np.array(time) * 86400.) - (MJDREFF * 86400) - (MJDREFFI * 86400)

def read_config(filename):
    """Read configuration from file."""
    with open(filename, 'r') as file:
        lines = file.readlines()
    t_min = int(lines[12].split()[1])
    t_max = int(lines[13].split()[1])
    source = ' '.join(lines[15].split()[1:3])
    source2 = '_'.join(lines[15].split()[1:3])
    return t_min, t_max, source, source2

def update_config(tmin, tmax):
    """Update config.yaml with new time range."""
    subprocess.run(f'sed -i "13s/.*/  tmin: {tmin}/" config.yaml', shell=True)
    subprocess.run(f'sed -i "14s/.*/  tmax: {tmax}/" config.yaml', shell=True)

def perform_analysis(gta, source, period, tmin, tmax):
    """Perform Fermi-LAT analysis."""
    gta.setup()
    gta.optimize()
    gta.free_sources(distance=10.0, pars='norm')
    gta.free_source('galdiff', pars='norm')
    gta.free_source('isodiff', pars='norm')
    gta.fit()
    sed = gta.sed(source, outfile='sed.fits')
    gta.write_roi(f"{source}_{period}_fit_{tmin}_{tmax}")
    return sed

def save_results(sed, results, tmin, tmax, source, period, j, count, home):
    """Save analysis results."""
    np.savetxt('sed.txt', np.c_[sed['dnde'], sed['e2dnde'], sed['e2dnde_err']], 
               delimiter=';', header='dnde;e2dnde;e2dnde_err')
    np.savetxt('spectral_index.txt', np.c_[sed['param_values'][1], sed['param_errors'][1],
                                           sed['param_values'][2], sed['param_errors'][2]], 
               delimiter=';', header='alpha;alpha_err;beta;beta_err')
    np.savetxt(f'lightcurve_{j}.txt', np.c_[(np.mean([met_to_mjd(tmin), met_to_mjd(tmax)])),
                                            count/2, results['eflux'][0], results['eflux_err'][0], count], 
               delimiter=';', header='time;time_err;flux;flux_error;binsize')
    os.chdir(home)
    with open('total_lightcurve.txt', 'a') as file:
        file.write(f'\n{np.mean([met_to_mjd(tmin), met_to_mjd(tmax)])};{count/2};'
                   f'{results["eflux"][0]};{results["eflux_err"][0]};{count};'
                   f'{sed["param_values"][1]};{sed["param_errors"][1]};'
                   f'{sed["param_values"][2]};{sed["param_errors"][2]}')

def main():
    unc = 1  # uncertainty - days
    Delta0 = 25  # ts threshold
    period = '2023'
    
    t_min, t_max, source, source2 = read_config('config1.yaml')
    home = os.getcwd()

    j = 0
    count = 1

    while (t_min + (j * unc * 86400)) <= t_max:
        tmin = t_min + unc * j * 86400
        tmax = t_min + unc * (j + count) * 86400
        
        update_config(tmin, tmax)
        
        try:
            gta = GTAnalysis('config.yaml', logging={'verbosity': 3}, 
                             fileio={'outdir': f"{tmin}_{tmax}"})
            sed = perform_analysis(gta, source, period, tmin, tmax)
            
            os.chdir(f"{tmin}_{tmax}")
            results = Table.read(f"{source}_{period}_fit_{tmin}_{tmax}.fits")
            
            if results['ts'][0] >= Delta0:
                save_results(sed, results, tmin, tmax, source, period, j, count, home)
                j += count
                count = 1
            else:
                print(f'\nBin size: {count}\n')
                count += 1
        
        except Exception as e:
            print(f'Error in analysis: {e}')
            with open('error_bins.txt', 'a') as file:
                file.write(f'\n{met_to_mjd(tmin)};{met_to_mjd(tmax)};{count}')
            count += 1
        
        print(f'\nLoading {((tmin - t_min) / (t_max - t_min)) * 100:.2f}% >> time {tmin}')
        os.chdir(home)
        
        if tmax >= t_max:
            break

if __name__ == "__main__":
    main()
