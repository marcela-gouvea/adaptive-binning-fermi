from fermipy.gtanalysis import GTAnalysis
import numpy as np
import subprocess
from matplotlib import pyplot as plt
from astropy.table import Table
import os
import shutil
import faulthandler

faulthandler.enable()

def met_to_mjd(time):
    #Convert mission elapsed time to mean julian date
    MJDREFF = 51910
    MJDREFFI = 7.428703703703703e-4
    
    if isinstance(time, (int,float)):
        # Convert from elapsed seconds to elapsed days
        elapsed_days = time / 86400.
        
        # Add the elapsed days to the refrence epoch
        MJD = MJDREFF + MJDREFFI + elapsed_days

    else:
        MJD = []
        for i in range(0,len(time)):
            # Convert from elapsed seconds to elapsed days
            elapsed_days = time[i] / 86400.
            MJD.append(MJDREFF + MJDREFFI + elapsed_days)

    # Return the MJD in TT format
    return MJD

def mjd_to_met(time):
    #Convert mean julian date to mission elapsed time
    MJDREFF = 51910 * 86400
    MJDREFFI = 7.428703703703703e-4 * 86400
    
    if isinstance(time, (int,float)):
        # Convert from elapsed days to elapsed seconds
        elapsed_seconds = time * 86400.
        
        # Add the elapsed seconds to the refrence epoch
        MET = elapsed_seconds - MJDREFF - MJDREFFI

    else:
        MET = []
        for i in range(0,len(time)):
            # Convert from elapsed days to elapsed seconds
            elapsed_seconds = time[i] * 86400.
            MET.append(elapsed_seconds - MJDREFF - MJDREFFI)

    # Return the MJD in TT format
    return MET


unc = 1 #uncertainty - days
Delta0 = 25 # ts 
directory = []


fname = 'config.yaml'

file = open('config1.yaml','r')
lista = file.readlines()
t_min = int(lista[12].split()[1])
t_max = int(lista[13].split()[1])

if len(lista[15].split()) == 3:
    source  = str(lista[15].split()[1])+' '+str(lista[15].split()[2])
    source2 = str(lista[15].split()[1])+'_'+str(lista[15].split()[2])
else:
    source  = str(lista[15].split()[1])
    source2 = str(lista[15].split()[1])

file.close()

period = '2023' # organize folders


home = os.getcwd() # get our current working directory and call it home

#n = (t_max-t_min)/(60*60*24*unc) #number of bins

j = 0 # iterator which defines tmin
count = 1 # iterator which defines tmax

while (t_min+(j*unc*86400)) <= t_max:

    try:    
        tmin = t_min + unc*j*86400
        tmax = t_min + unc*(j+count)*86400
        
        subprocess.run('sed -i \'13s/.*/  tmin: {}/\' config.yaml'.format(tmin), shell=True)
        subprocess.run('sed -i \'14s/.*/  tmax: {}/\' config.yaml'.format(tmax), shell=True)
        
        print('                                                          START ANALYSIS                                                                     \n')
        
        gta = GTAnalysis('config.yaml', logging={'verbosity': 3}, fileio={'outdir': "{}_{}".format(tmin,tmax)}) #define our gta object, but with the tiperiod from our list
        gta.setup() #photon selection, good time intervals, livetime cube, binning etc
        
        print('\n                                                          PASS GTA.SETUP()                                                                     \n')
        
        gta.optimize() #initial optimise
        gta.free_sources(distance=10.0,pars='norm') #frees the point sources
        gta.free_source('galdiff', pars='norm') #frees the galactic diffuse
        gta.free_source('isodiff', pars='norm') #frees the isotropic diffuse
        gta.fit() #full likelihood fit
        gta.sed('{}'.format(source))#, make_plots="True") #do an SED, we'll explain why shortly
        
        print('\n                                                          PASS SED                                                                     \n')
        
        sed = gta.sed('{}'.format(source), outfile='sed.fits')
        gta.write_roi("{}_{}_fit_{}_{}".format(source, period,tmin,tmax)) #save our ROI
        
        os.chdir("{}_{}".format(tmin,tmax)) # open the directory
            
        # test bin size 
        results = Table.read("{}_{}_fit_{}_{}".format(source,period,tmin,tmax) + ".fits")

        if (results['ts'][0]) >= Delta0:
            print("\n                                                          ACCEPTANCE CRITERIA SATISFIED                                                                     \n'")
        
            # save the SED values
            np.savetxt('sed.txt', np.c_[sed['dnde'],sed['e2dnde'],sed['e2dnde_err']], delimiter=';', header='dnde;e2dnde;e2dnde_err')
            
            # save the spectral index values
            np.savetxt('spectral_index.txt', np.c_[sed['param_values'][1],sed['param_errors'][1],sed['param_values'][2], sed['param_errors'][2]], delimiter=';', header='alpha;alpha_err;beta;beta_err')
                
            # save the lightcurve values
            #results = Table.read("{}_{}_fit_{}".format(source, period,j) + ".fits")
            np.savetxt('lightcurve_{}.txt'.format(j), np.c_[( np.mean([met_to_mjd(tmin), met_to_mjd(tmax)]) ), unc*count/2, results['eflux'][0], results['eflux_err'][0], count], delimiter=';', header='time;time_err;flux;flux_error;binsize')
                
            os.chdir(home) # close the directory
                
            # save the lightcurve values in the final txt file
            file2 = open('total_lightcurve.txt', 'a')
            file2.write('\n{};{};{};{};{};{};{};{};{}'.format( np.mean([met_to_mjd(tmin), met_to_mjd(tmax)]), unc*count/2, results['eflux'][0] , results['eflux_err'][0], count*unc, sed['param_values'][1], sed['param_errors'][1], sed['param_values'][2], sed['param_errors'][2] ))
            file2.close()

            j += count
            count = 1

        else:
            print('\n------------------------------------------------        {}        ------------------------------------------------'.format(count))
            print('------------------------------------------------        {}        ------------------------------------------------'.format(count))
            print('------------------------------------------------        {}        ------------------------------------------------'.format(count))
            print('------------------------------------------------        {}        ------------------------------------------------'.format(count))
            print('------------------------------------------------        {}        ------------------------------------------------\n'.format(count))
            print("Bin size: {} \n".format(count))
            count += 1

        
        print('                                                          END ANALYSIS                                                                     ')
    
    except:
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')
        print('------------------------------------------------        ERROR ANALYSIS        ------------------------------------------------')

        # save the error bin values
        file2 = open('error_bins.txt', 'a')
        file2.write('\n{};{};{}'.format( met_to_mjd( tmin ) , met_to_mjd(tmax), count ))
        file2.close()
        count += 1
    
    print('\n---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    print('---------------------------------------------        LOADING {} % >> time {}        ---------------------------------------------\n\n'.format( ( (tmin - t_min) / (t_max-t_min) )*100, tmin ))
    
    os.chdir(home)
    
    if tmax>=t_max:
        break
