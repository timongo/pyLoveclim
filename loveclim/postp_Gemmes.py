"""                                                                             
file name: postp_Gemmes.py
language: Python3.6
date: 31/01/2021
author: Hugo Martin
email: martin.hugo@live.com
This module contains functions for simply post processing Gemmes.
"""

### imports -------------------------------------------------------------------
from loveclim.gui import np, os, plt
import sys
from matplotlib.ticker import AutoMinorLocator
import matplotlib as mpl
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
pad = 0.
w_pad = 0.
h_pad = 0.
size = 10
SMALL_SIZE = size
MEDIUM_SIZE = size
BIGGER_SIZE = size
plt.rc('font', size = SMALL_SIZE) # controls default text sizes
plt.rc('axes', titlesize = BIGGER_SIZE) # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE) # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE) # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE) # fontsize of the tick labels
plt.rc('legend', fontsize = MEDIUM_SIZE) # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE) # fontsize of the figure title
mpl.rcParams['lines.linewidth'] = 0.8
mpl.rcParams['lines.markersize'] = 1.5

### functions -----------------------------------------------------------------
# load_Gemmes_data
def load_Gemmes_data(path='.', y_ini=2016, y_stop=2030):
    """
    This function can be used to extract Gemmes data
    PARAMETERS:
        path : string : path where there are init_sim_*.txt files.
        y_ini : int : first year to load
        y_stop : int : last year to load

    OUTPUTS:
       Data_dict : dict : dictionary with all Gemmes outputs.
    """

    # initialization
    years = np.arange(start=y_ini, stop=y_stop+1, step=1)    
    Data = np.zeros((years.size, 18))
    Data[:,17] = years.copy()

    # time loop
    for y, year in enumerate(years):
        # load init_sim
        name = os.path.join(path, 'init_sim_'+str(year)+'.txt')
        file1 = open(name, 'r')
        Lines = file1.readlines()
        file1.close()
        data = np.zeros(15)
        for k, l in enumerate(Lines):
            data[k] = float(l[1:-2])
        Data[y,:15] = data.copy()

        # load t2mlov
        name = os.path.join(path, 't2mlov_'+str(year)+'.txt')
        file1 = open(name, 'r')
        Lines = file1.readlines()
        file1.close()
        Data[y,15] = float(Lines[0])

        # load emissions
        name = os.path.join(path, 'emissions_'+str(year)+'.txt')
        E = np.loadtxt(name)
        Data[y,16] = E

    # make dict
    labels = ['Y', 'N', 'NG', 'lambda', 'omega', 'd', 'Eind', 'Eland', 'n', 
              'sigma', 'Gsigma', 'Pbs', 'p', 'a', 'K', 'T', 'E', 'time']

    Data_dict = {key:Data[:,j] for j, key in enumerate(labels)}
    
    print('Gemmes data extracted.')

    return Data_dict

# Gemmes_quick_view
def Gemmes_quick_view(namef, path='.', y_ini=2016, y_stop=2100, labels=None):
    """
    Function that can be used for quick view of a typicall Gemmes/Iloveclim
    simulation. A figure is created containing the fields : Y, T, omega,
    lambda, E and d. If len(path)=1 then simple plot, else, comparison. The
    parameters must of the same length.
    """
    # test on path
    try:
        testlen = len(path)>1
        if not(len(path)==len(y_ini) and len(path)==len(y_stop)
           and len(path)==len(labels)):
           raise ValueError('Arguments list must be of same length.')

    except TypeError:
        path = [path]
        y_ini = [y_ini]
        y_stop = [y_stop]
        labels = [labels]

    # loop on paths
    dataDict = []
    for Np, p in enumerate(path):
        dataDict.append(load_Gemmes_data(path=p, y_ini=y_ini[Np],
        y_stop=y_stop[Np]))

    # make figure
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(6.5, 4)) 
    for Np, lab in enumerate(labels):
        axes[0,0].plot(dataDict[Np]['time'], dataDict[Np]['Y'], label=lab)
        axes[0,1].plot(dataDict[Np]['time'], dataDict[Np]['T'], label=lab)
        axes[0,2].plot(dataDict[Np]['time'], dataDict[Np]['omega'], label=lab)
        axes[1,0].plot(dataDict[Np]['time'], dataDict[Np]['lambda'], label=lab)
        axes[1,1].plot(dataDict[Np]['time'], dataDict[Np]['E'], label=lab)
        axes[1,2].plot(dataDict[Np]['time'], dataDict[Np]['d'], label=lab)

    # labels
    axes = np.asarray(axes)
    # x labels
    for i in range(axes.shape[0]):
        for j in range(axes.shape[1]):
            axes[i,j].set_xlabel('$t$ (year)')
            axes[i,j].xaxis.set_minor_locator(AutoMinorLocator())
            axes[i,j].yaxis.set_minor_locator(AutoMinorLocator())
            axes[i,j].grid(which='both', linewidth=.0001, color='silver',
            alpha=0.4)

    # y labels
    axes[0,0].set_ylabel(r'$Y$')
    axes[0,1].set_ylabel(r'$T~(^\circ C)$')
    axes[0,2].set_ylabel(r'$\omega$')
    axes[1,0].set_ylabel(r'$\lambda$')
    axes[1,1].set_ylabel(r'$E$ (Gt CO$_2$)')
    axes[1,2].set_ylabel(r'$d$')

    # legend
    axes[0,0].legend()
    plt.tight_layout(pad, h_pad, w_pad)
    plt.savefig(namef)
    plt.close(fig)
