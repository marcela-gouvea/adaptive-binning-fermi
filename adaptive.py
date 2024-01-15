from fermipy.gtanalysis import GTAnalysis
import numpy as np
import subprocess
from matplotlib import pyplot as plt #matplotlib.pyplot is the package we will use to plot our light-curve.
from astropy.table import Table #astropy.table allows us to read fits tables, which is the structure of our ROI file. 
import os #os lets us easily perform file manipulation within Python.

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


def mean_err(array):
    if len(array) != 0:
        soma = 0
        for k in range(len(array)):
            soma += array[k]**2
        return np.sqrt(soma/len(array))
    

unc = 1 #uncertainty - days
Delta0 = 25
directory = []


fname = 'config.yaml'

#source = 'PKS 1510-089'
#source2 = 'PKS-1510-089'

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

mes = '2023' # organize folders


home = os.getcwd() # get our current working directory and call it home

n = (t_max-t_min)/(60*60*24*unc)       #number of bins

j = 0
count = 1

while (t_min+(j*unc*86400)) <= t_max:

    try:    
        tmin = t_min + unc*j*86400
        tmax = t_min + unc*(j+1*count)*86400
        #print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', j)
        
        subprocess.run('sed -i \'13s/.*/  tmin: {}/\' config.yaml'.format(tmin), shell=True)
        subprocess.run('sed -i \'14s/.*/  tmax: {}/\' config.yaml'.format(tmax), shell=True)
        #directory = "{}_{}_{}".format(source2,mes,j)

        print('                                                          START ANALYSIS                                                                     \n')
        
        gta = GTAnalysis('config.yaml', logging={'verbosity': 3}, fileio={'outdir': "{}_{}_{}".format(source2,mes,tmax)}) #define our gta object, but with the times from our list
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
        gta.write_roi("{}_{}_fit_{}".format(source, mes,tmax)) #save our ROI
        
        os.chdir("{}_{}_{}".format(source2,mes,tmax)) # open the directory
            
        # test bin size 
        results = Table.read("{}_{}_fit_{}".format(source, mes,tmax) + ".fits")

        if (results['ts'][0]) >= Delta0:
            print("\n                                                          ACCEPTANCE CRITERIA SATISFIED                                                                     \n'")

            directory.append("{}_{}_{}".format(source2,mes,tmax))
        
            # save the SED values
            np.savetxt('sed.txt', np.c_[sed['dnde'],sed['e2dnde'],sed['e2dnde_err']], delimiter=';', header='dnde;e2dnde;e2dnde_err')
            
                # save the spectral index values
            np.savetxt('spectral_index.txt', np.c_[sed['param_values'][1],sed['param_errors'][1],sed['param_values'][2], sed['param_errors'][2]], delimiter=';', header='alpha;alpha_err;beta;beta_err')
                
                # save the lightcurve values
            #results = Table.read("{}_{}_fit_{}".format(source, mes,j) + ".fits")
            np.savetxt('lightcurve_{}.txt'.format(j), np.c_[met_to_mjd( np.mean([tmin, tmax]) ), unc*count/2, results['eflux'][0], results['eflux_err'][0], count], delimiter=';', header='time;time_err;flux;flux_error;binsize')
                
            os.chdir(home) # close the directory
                
                # save the lightcurve values in the final txt file
            file2 = open('total_lightcurve.txt', 'a')
            file2.write('\n{};{};{};{};{}'.format( met_to_mjd( np.mean([tmin, tmax]) ), unc*count/2, results['eflux'][0] , results['eflux_err'][0], count*unc ))
            file2.close()

            j += 1
            count = 1

        else:
            print('\n------------------------------------------------        {}        ------------------------------------------------'.format(count))
            print('------------------------------------------------        {}        ------------------------------------------------'.format(count))
            print('------------------------------------------------        {}        ------------------------------------------------'.format(count))
            print('------------------------------------------------        {}        ------------------------------------------------'.format(count))
            print('------------------------------------------------        {}        ------------------------------------------------\n'.format(count))
            print("Bin size: {} \n".format(count))
            count += 1
            
            #subprocess.run('rm -r {}_{}_{}'.format(source2,mes,j) , shell=True)
            
            #if tmax >= t_max:
            #    break

        
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


# plot the lightcurve
time      = np.loadtxt('total_lightcurve.txt',delimiter=';')[:,0]
time_err  = np.loadtxt('total_lightcurve.txt',delimiter=';')[:,1]
eflux     = np.loadtxt('total_lightcurve.txt',delimiter=';')[:,2]
eflux_err = np.loadtxt('total_lightcurve.txt',delimiter=';')[:,3]

## FIXING BINS

new_flux = []
new_flux_err = []
new_time = []
new_time_err1 = []
new_time_err2 = []
central_time = []
central_time_err = []

aux_flux = []
aux_flux_err = []
aux_time = []
aux_time_err = []

aux_bins = []

i = 0 
count = 0

while i<len(eflux):
    if (i+count+1)>= len(eflux):
            break
    count = 0
    while ((eflux[i+count]/eflux[i+count+1] <= 1.2) & (eflux[i+count]/eflux[i+count+1] >= 0.8)):
        count += 1
    if count == 0:
        aux_flux.append(eflux[i])
        aux_time.append(time[i])
        aux_flux_err.append(eflux_err[i])
        aux_time_err.append(time_err[i])
        aux_bins.append(time[i]-time_err[i])
        aux_bins.append(time[i]+time_err[i])
    else:
        for j in range(count):
            aux_flux.append(eflux[i+j])
            aux_time.append(time[i+j])
            aux_flux_err.append(eflux_err[i+j])
            aux_time_err.append(time_err[i+j])
            aux_bins.append(time[i+j]-time_err[i+j])
            aux_bins.append(time[i+j]+time_err[i+j])

    new_flux.append(np.mean(aux_flux))
    new_flux_err.append(mean_err(aux_flux_err))
    new_time.append(np.mean(aux_time))
    new_time_err1.append(new_time[-1]-min(aux_bins))
    new_time_err2.append(max(aux_bins)-new_time[-1])
    central_time.append((2*new_time[-1]+new_time_err2[-1]-new_time_err1[-1])/2)
    central_time_err.append((new_time_err1[-1]+new_time_err2[-1])/2)
    
    aux_flux.clear()
    aux_flux_err.clear()
    aux_time.clear()
    aux_time_err.clear()
    aux_bins.clear()
    i += (count + 1)


np.savetxt("fixed_bins.txt", np.c_[new_time, new_time_err1, new_time_err2, central_time, central_time_err, new_flux, new_flux_err], delimiter=';', header="assym_time;assym_time_err_1;assym_time_err_2;central_time;central_time_error;flux;flux_err")

fig, ax = plt.subplots()
plt.errorbar(x = central_time, xerr = central_time_err, y = new_flux , yerr = new_flux_err , fmt='o' , ecolor='k', capsize = 5, markersize=4)
plt.legend(loc='best',numpoints=1)
#plt.xlim([met_to_mjd(t_start), met_to_mjd(t_end)])
plt.grid(linestyle = '--',linewidth = 0.5)
plt.xlabel('Time (MJD)')
plt.ylabel(r'Energy Flux ($MeV~cm^{-2}~s^{-1}$)')
plt.title("{} Fermi-LAT lightcurve - {}".format(source, mes))
ax.set_yscale('log')
#ax.set_xscale('log')
plt.savefig('lightcurve_{}_{}.png'.format(source2, mes))
plt.show()
plt.close()

all_folders = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
delete = list(set(all_folders) - set(directory))

for k in delete:
    subprocess.run('rm -r {}'.format(k), shell=True)

for p in range(0, len(central_time)):
    tmin = mjd_to_met(central_time[p]-central_time_err[p])
    tmax = mjd_to_met(central_time[p]+central_time_err[p])

    subprocess.run('sed -i \'13s/.*/  tmin: {}/\' config.yaml'.format(tmin), shell=True)
    subprocess.run('sed -i \'14s/.*/  tmax: {}/\' config.yaml'.format(tmax), shell=True)
    print('                                                          START ANALYSIS                                                                     \n')
        
    gta = GTAnalysis('config.yaml', logging={'verbosity': 3}, fileio={'outdir': "{}/final-sed/{}_{}_{}".format(home,source2,mes,tmax)}) #define our gta object, but with the times from our list
    gta.setup() #photon selection, good time intervals, livetime cube, binning etc
    
    print('\n                                                          PASS GTA.SETUP()                                                                     \n')
        
    gta.optimize() #initial optimise
    gta.free_sources(distance=10.0,pars='norm') #frees the point sources
    gta.free_source('galdiff', pars='norm') #frees the galactic diffuse
    gta.free_source('isodiff', pars='norm') #frees the isotropic diffuse
    gta.fit() #full likelihood fit
    gta.sed('{}'.format(source))#, make_plots="True") #do an SED
    
    print('\n                                                          PASS SED                                                                     \n')
    
    sed = gta.sed('{}'.format(source), outfile='sed.fits')
    gta.write_roi("{}_{}_fit_{}".format(source, mes,tmax)) #save our ROI
    
    os.chdir("{}_{}_{}".format(source2,mes,tmax)) # open the directory
        
    # test bin size 
    results = Table.read("{}_{}_fit_{}".format(source, mes,tmax) + ".fits")

    # save the SED values
    np.savetxt('sed.txt', np.c_[sed['dnde'],sed['e2dnde'],sed['e2dnde_err']], delimiter=';', header='dnde;e2dnde;e2dnde_err')
            
    # save the spectral index values
    np.savetxt('spectral_index.txt', np.c_[sed['param_values'][1],sed['param_errors'][1],sed['param_values'][2], sed['param_errors'][2]], delimiter=';', header='alpha;alpha_err;beta;beta_err')

aaa = os.listdir("{}/final-sed".format(home))
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
