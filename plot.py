import numpy as np
import matplotlib.pyplot as plt
import os

def compatibility(x1, x1err, x2, x2err, time, time_err, n):
    y1 = []
    y1err = []
    y2 = []
    y2err = []
    new_time = []
    new_time_err = []
    if n == 1: # All the non zero compatible points -- Clean the spectral evolution
        for i in range(0, len(x1)):
            if abs(x1[i]/x1err[i]) > 3:
                y1.append(x1[i])
                y1err.append(x1err[i])
                new_time.append(time[i])
                new_time_err.append(time_err[i])
                y2.append(x2[i])
                y2err.append(x2err[i])
        return y1, y1err, y2, y2err, new_time, new_time_err
    if n == 0: # All the zero compatible points -- Simple Powerlaw
        for i in range(0, len(x1)):
            if abs(x1[i]/x1err[i]) < 3:
                y1.append(x1[i])
                y1err.append(x1err[i])
                new_time.append(time[i])
                new_time_err.append(time_err[i])
                y2.append(x2[i])
                y2err.append(x2err[i])
        return y1, y1err, y2, y2err,  new_time, new_time_err


fname = 'config.yaml'

file = open('config.yaml','r')
lista = file.readlines()

if len(lista[15].split()) == 3:
    source  = str(lista[15].split()[1])+' '+str(lista[15].split()[2])
    source2 = str(lista[15].split()[1])+'_'+str(lista[15].split()[2])
else:
    source  = str(lista[15].split()[1])
    source2 = str(lista[15].split()[1])

file.close()

period = '2023' # organize folders

home = os.getcwd()

alpha     = np.loadtxt(home+'/spectral_index/result.csv', delimiter=',')[:,0]
alpha_err = np.loadtxt(home+'/spectral_index/result.csv', delimiter=',')[:,1]

beta      = np.loadtxt(home+'/spectral_index/result.csv', delimiter=',')[:,2]
beta_err  = np.loadtxt(home+'/spectral_index/result.csv', delimiter=',')[:,3]

flux      = np.loadtxt('lightcurve_2023.txt', delimiter=';')[:,1]
flux_err  = np.loadtxt('lightcurve_2023.txt', delimiter=';')[:,2]

central_time     = np.loadtxt(home+'/spectral_index/result.csv', delimiter=',')[:,4]
central_time_err = [1.5]*len(central_time)


adaptive_alpha      = np.loadtxt(home+'/adaptive/temp/spectral_index.txt',delimiter=';')[:,0]
adaptive_alpha_err  = np.loadtxt(home+'/adaptive/temp/spectral_index.txt',delimiter=';')[:,1]

adaptive_beta      = np.loadtxt(home+'/adaptive/temp/spectral_index.txt',delimiter=';')[:,2]
adaptive_beta_err  = np.loadtxt(home+'/adaptive/temp/spectral_index.txt',delimiter=';')[:,3]

adaptive_flux      = np.loadtxt(home+'/adaptive/temp/fixed_bins.txt', delimiter=';')[:,2]
adaptive_flux_err  = np.loadtxt(home+'/adaptive/temp/fixed_bins.txt', delimiter=';')[:,3]

adaptive_central_time      = np.loadtxt(home+'/adaptive/temp/fixed_bins.txt',delimiter=';')[:,0]
adaptive_central_time_err  = np.loadtxt(home+'/adaptive/temp/fixed_bins.txt',delimiter=';')[:,1]


alpha_err = [abs(x) for x in alpha_err]
beta_err  = [abs(x) for x in beta_err]

adaptive_alpha_err = [abs(x) for x in adaptive_alpha_err]
adaptive_beta_err  = [abs(x) for x in adaptive_beta_err]


# SIMPLE POWER LAW REGULAR BINNING

beta_pl, beta_err_pl, alpha_pl, alpha_err_pl, central_time_pl, central_time_err_pl = compatibility(beta, beta_err, alpha, alpha_err, central_time, central_time_err, 0)

alpha_pl, alpha_err_pl, beta_pl, beta_err_pl, central_time_pl, central_time_err_pl = compatibility(alpha_pl, alpha_err_pl, beta_pl, beta_err_pl, central_time_pl, central_time_err_pl, 1)

# SIMPLE POWER LAW ADAPTIVE BINNING

adaptive_beta_pl, adaptive_beta_err_pl, adaptive_alpha_pl, adaptive_alpha_err_pl, adaptive_central_time_pl, adaptive_central_time_err_pl = compatibility(adaptive_beta, adaptive_beta_err, adaptive_alpha, adaptive_alpha_err, adaptive_central_time, adaptive_central_time_err, 0)

adaptive_alpha_pl, adaptive_alpha_err_pl, adaptive_beta_pl, adaptive_beta_err_pl, adaptive_central_time_pl, adaptive_central_time_err_pl = compatibility(adaptive_alpha_pl, adaptive_alpha_err_pl, adaptive_beta_pl, adaptive_beta_err_pl, adaptive_central_time_pl, adaptive_central_time_err_pl, 1)



fig, ax = plt.subplots(2, sharex=True, figsize=(7,6))
plt.subplots_adjust(wspace=0, hspace=0)
fig.suptitle("{} Fermi-LAT Light curve and Spectral Index - {}".format(source, period), y=0.93)
fig.supxlabel("Time (MJD)")
fig.align_labels() 

ax[0].errorbar(x = central_time, xerr = central_time_err, y = flux , yerr = flux_err , fmt='o' , ecolor='k', capsize = 5, markersize=3)
ax[0].legend(loc='best',numpoints=1)
#ax[0].xlim([met_to_mjd(t_start), met_to_mjd(t_end)])
ax[0].grid(linestyle = '--',linewidth = 0.5)
ax[0].set_yscale('log')
ax[0].set_ylabel(r'Energy Flux ($MeV~cm^{-2}~s^{-1}$)')

ax[1].errorbar(x = central_time_pl, xerr = central_time_err_pl, y = alpha_pl , yerr = alpha_err_pl , fmt='o' , color='k', capsize = 5, markersize=4)
ax[1].legend(loc='best',numpoints=1)
#ax[0].xlim([met_to_mjd(t_start), met_to_mjd(t_end)])
ax[1].grid(linestyle = '--',linewidth = 0.5)
ax[1].set_ylabel(r'Spectral Index ($\alpha$)')


plt.savefig('lc_spectralindex_{}_{}.png'.format(source2, period))
plt.show()



fig, ax = plt.subplots(2, sharex=True, figsize=(7,6))
plt.subplots_adjust(wspace=0, hspace=0)
fig.suptitle("{} Fermi-LAT Light curve and Spectral Index - {}".format(source, period), y=0.93)
fig.supxlabel("Time (MJD)")
fig.align_labels() 

ax[0].errorbar(x = adaptive_central_time, xerr = adaptive_central_time_err, y = adaptive_flux , yerr = adaptive_flux_err , fmt='o' , ecolor='k', capsize = 5, markersize=3)
ax[0].legend(loc='best',numpoints=1)
#ax[0].xlim([met_to_mjd(t_start), met_to_mjd(t_end)])
ax[0].grid(linestyle = '--',linewidth = 0.5)
ax[0].set_yscale('log')
ax[0].set_ylabel(r'Energy Flux ($MeV~cm^{-2}~s^{-1}$)')

ax[1].errorbar(x = adaptive_central_time_pl, xerr = adaptive_central_time_err_pl, y = adaptive_alpha_pl , yerr = adaptive_alpha_err_pl , fmt='o' , color='k', capsize = 5, markersize=4)
ax[1].legend(loc='best',numpoints=1)
#ax[0].xlim([met_to_mjd(t_start), met_to_mjd(t_end)])
ax[1].grid(linestyle = '--',linewidth = 0.5)
ax[1].set_ylabel(r'Spectral Index ($\alpha$)')


plt.savefig('adaptive_lc_spectralindex_{}_{}.png'.format(source2, period))
plt.show()