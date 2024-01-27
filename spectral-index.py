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

home = os.getcwd() # get our current working directory and call it home

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

central_time      = np.loadtxt('total_lightcurve_3C454.3.txt',delimiter=';')[:,0]
central_time_err  = np.loadtxt('total_lightcurve_3C454.3.txt',delimiter=';')[:,1]

'''
for pasta in os.listdir(home):
    caminho_pasta = os.path.join(home, pasta)
    if os.path.isdir(caminho_pasta):
        arquivo_sed = os.path.join(caminho_pasta, 'sed.txt')
        if not os.path.exists(arquivo_sed):
            shutil.rmtree(caminho_pasta)
'''
 
for p in range(0, len(central_time)):
    tmin = mjd_to_met(central_time[p]-central_time_err[p])
    tmax = mjd_to_met(central_time[p]+central_time_err[p])

    #if (tmin>=t_min) and (tmax<=t_max):
    subprocess.run('sed -i \'13s/.*/  tmin: {}/\' config.yaml'.format(tmin), shell=True)
    subprocess.run('sed -i \'14s/.*/  tmax: {}/\' config.yaml'.format(tmax), shell=True)
    print('                                                          START ANALYSIS                                                                     \n')
            
    gta = GTAnalysis('config.yaml', logging={'verbosity': 3}, fileio={'outdir': "{}_{}".format(int(tmin),int(tmax))}) #define our gta object, but with the tiperiod from our list
    gta.setup() #photon selection, good time intervals, livetime cube, binning etc
    
    print('\n                                                          PASS GTA.SETUP()                                                                     \n')
    
    gta.optimize() #initial optimise
    gta.free_sources(distance=10.0,pars='norm') #frees the point sources
    gta.free_source('galdiff', pars='norm') #frees the galactic diffuse
    gta.free_source('isodiff', pars='norm') #frees the isotropic diffuse
    gta.fit() #full likelihood fit
    gta.sed(source)#, make_plots="True") #do an SED
    
    print('\n                                                          PASS SED                                                                     \n')
    
    sed = gta.sed(source, outfile='sed.fits')
    gta.write_roi("{}_{}_fit_{}_{}".format(source, period,int(tmin),int(tmax))) #save our ROI
    os.chdir("{}_{}".format(int(tmin),int(tmax))) # open the directory
            
    # test bin size 
    results = Table.read("{}_{}_fit_{}_{}".format(source, period,int(tmin),int(tmax)) + ".fits")

    # save the SED values
    np.savetxt('sed.txt', np.c_[sed['dnde'],sed['e2dnde'],sed['e2dnde_err']], delimiter=';', header='dnde;e2dnde;e2dnde_err')
            
    # save the spectral index values
    np.savetxt('spectral_index.txt', np.c_[sed['param_values'][1],sed['param_errors'][1],sed['param_values'][2], sed['param_errors'][2]], delimiter=';', header='alpha;alpha_err;beta;beta_err')
    
    os.chdir(home)
    
    
aaa = os.listdir(home)#"{}/final-sed".format(home))
directory = []

for i in range(0,len(aaa)):
    if os.path.isdir('{}/{}'.format(home,aaa[i]))==True:
        directory.append(aaa[i])


subprocess.run('mkdir sed', shell=True)


for j in range(0,len(directory)):
        os.chdir(directory[j])
        subprocess.run('mv sed.txt sed_{}.txt'.format(directory[j]), shell=True)
        subprocess.run('cp sed_{}.txt {}/sed'.format(directory[j],home), shell=True)
        os.chdir(home)


subprocess.run('mkdir spectral_index', shell=True)


for j in range(0,len(directory)):
        os.chdir(directory[j])
        subprocess.run('mv spectral_index.txt spectral_index_{}.txt'.format(directory[j]), shell=True)
        subprocess.run('cp spectral_index_{}.txt {}/spectral_index'.format(directory[j],home), shell=True)
        os.chdir(home)
