import matplotlib.pyplot as plt
import numpy as np 
from scipy.optimize import curve_fit
from scipy import stats
import os
import pathlib 

path = str(pathlib.Path().absolute())
path_energies = path + '/energies'

files = os.listdir(path_energies)

OPLS_normal_min = []
OPLS_normal_av = []
OPLS_normal_bz = []
OPLS_normal_st = []
OPLS_normal_stbz = []
OPLS_normal_bz_corr = []

OPLS_fast_min = []
OPLS_fast_av = []
OPLS_fast_bz = []
OPLS_fast_st = []
OPLS_fast_stbz = []
OPLS_fast_bz_corr = []

OFF_normal_min = []
OFF_normal_av = []
OFF_normal_bz = []
OFF_normal_st = []
OFF_normal_stbz = []
OFF_normal_bz_corr = []

OFF_fast_min = []
OFF_fast_av = []
OFF_fast_bz = []
OFF_fast_st = []
OFF_fast_stbz = []
OFF_fast_bz_corr = []

for document in files:

    if 'OFfast' in document:

        with open(os.path.join(path_energies,document)) as filein:

            for line in filein:

                line = line.split(',')
                OFF_fast_min.append(float(line[0]))
                OFF_fast_av.append(float(line[1]))
                OFF_fast_bz.append(float(line[2]))
                OFF_fast_st.append(float(line[3]))
                OFF_fast_stbz.append(float(line[4]))
                OFF_fast_bz_corr.append(float(line[5]))

    elif 'OFnormal' in document:

        with open(os.path.join(path_energies,document)) as filein:

            for line in filein:

                line = line.split(',')
                OFF_normal_min.append(float(line[0]))
                OFF_normal_av.append(float(line[1]))
                OFF_normal_bz.append(float(line[2]))
                OFF_normal_st.append(float(line[3]))
                OFF_normal_stbz.append(float(line[4]))
                OFF_normal_bz_corr.append(float(line[5]))

    elif 'OPLSnormal' in document:

        with open(os.path.join(path_energies,document)) as filein:

            for line in filein:

                line = line.split(',')
                OPLS_normal_min.append(float(line[0]))
                OPLS_normal_av.append(float(line[1]))
                OPLS_normal_bz.append(float(line[2]))
                OPLS_normal_st.append(float(line[3]))
                OPLS_normal_stbz.append(float(line[4]))
                OPLS_normal_bz_corr.append(float(line[5]))

    elif 'OPLSfast' in document:

        with open(os.path.join(path_energies,document)) as filein:

            for line in filein:

                line = line.split(',')
                OPLS_fast_min.append(float(line[0]))
                OPLS_fast_av.append(float(line[1]))
                OPLS_fast_bz.append(float(line[2]))
                OPLS_fast_st.append(float(line[3]))
                OPLS_fast_stbz.append(float(line[4]))
                OPLS_fast_bz_corr.append(float(line[5]))

    else:
        continue

DG = np.array([-7.6,-12.3,-13.1,-13.3,-14.9,-15.0,-8.3,-7.8])
DH = np.array([0.0,0.2,1.3,-6.7,-8.0,-12.7,-14.4,-13.3])

protocols = [OPLS_normal_min,OPLS_normal_av,OPLS_normal_bz,OPLS_normal_st,OPLS_normal_stbz,OPLS_normal_bz_corr,OPLS_fast_min,\
    OPLS_fast_av,OPLS_fast_bz,OPLS_fast_st,OPLS_fast_stbz,OPLS_fast_bz_corr,OFF_normal_min,OFF_normal_av,\
    OFF_normal_bz,OFF_normal_st,OFF_normal_stbz,OFF_normal_bz_corr,OFF_fast_min,OFF_fast_av,OFF_fast_bz,OFF_fast_st,\
    OFF_fast_stbz,OFF_fast_bz_corr]

energies = [DG,DH]

protocols_names = ['OPLS_normal_min','OPLS_normal_av','OPLS_normal_bz','OPLS_normal_st','OPLS_normal_stbz','OPLS_normal_bz_corr','OPLS_fast_min',\
    'OPLS_fast_av','OPLS_fast_bz','OPLS_fast_st','OPLS_fast_stbz','OPLS_fast_bz_corr','OFF_normal_min','OFF_normal_av',\
    'OFF_normal_bz','OFF_normal_st','OFF_normal_stbz','OFF_normal_bz_corr','OFF_fast_min','OFF_fast_av','OFF_fast_bz','OFF_fast_st',\
    'OFF_fast_stbz','OFF_fast_bz_corr']

energies_names = ['DG','DH']
energies_names_plot = ['$\Delta G_{exp}$','$\Delta H_{exp}$']
scoring_names = ['$\Delta G^{ff}_{min}$','$\Delta G^{ff}_{av}$','$\Delta G^{ff}_{Bz}$','$\Delta G^{ff}_{step}$','$\Delta G^{ff}_{sBz}$','$\Delta G^{corr}_{Bz}$',\
    '$\Delta G^{ff}_{min}$','$\Delta G^{ff}_{av}$','$\Delta G^{ff}_{Bz}$','$\Delta G^{ff}_{step}$','$\Delta G^{ff}_{sBz}$','$\Delta G^{corr}_{Bz}$',\
    '$\Delta G^{ff}_{min}$','$\Delta G^{ff}_{av}$','$\Delta G^{ff}_{Bz}$','$\Delta G^{ff}_{step}$','$\Delta G^{ff}_{sBz}$','$\Delta G^{corr}_{Bz}$',\
    '$\Delta G^{ff}_{min}$','$\Delta G^{ff}_{av}$','$\Delta G^{ff}_{Bz}$','$\Delta G^{ff}_{step}$','$\Delta G^{ff}_{sBz}$','$\Delta G^{corr}_{Bz}$']

def lin(x,m,n):
    y = m*x + n
    return y

cont_protocol = 0
cont_energies = 0

path = str(pathlib.Path().absolute())

if  os.path.exists(path+'/images') == False:
    os.mkdir(path+'/images')

for protocol in protocols:
    for energy in energies:

        param, cov = curve_fit(lin, energy, protocol)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(energy,protocol)

        plt.figure()
        plt.title(protocols_names[cont_protocol]+'_'+energies_names[cont_energies])
        plt.scatter(energy,protocol,label='r(W):'+str("{:10.2f}".format(r_value )))
        plt.plot(energy,lin(energy,param[0],param[1]))
        plt.xlabel(energies_names_plot[cont_energies])
        plt.ylabel(scoring_names[cont_protocol])
        plt.legend(loc='best')
        plt.savefig(path+'/images/'+protocols_names[cont_protocol]+'_'+energies_names[cont_energies]+'.png', format='png')
        plt.close()

        cont_energies += 1

    cont_protocol += 1
    cont_energies = 0