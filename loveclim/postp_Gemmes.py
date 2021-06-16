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
mpl.rcParams['lines.linewidth'] = 1.
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
def Gemmes_quick_view(namef, path='.', y_ini=2016, y_stop=2100, labels='', returnfig=False):
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

    # If at least one label is '', then label is the mean temperature in 2100
    labbool = False
    for lab in labels:
        if lab=='':
            labbool = True
    if labbool:
        lenlab = len(labels)
        labels = []
        for Np in range(lenlab):
            time = dataDict[Np]['time']
            temp = dataDict[Np]['T']
            ind = (time >= 2095)*(time <= 2105)
            inddd = (time >= 2065)*(time <= 2075)
            labels.append('+{:.1f}°C'.format(np.mean(temp[ind])))

    # make figure
    if False:
        # an option where only few Gemmes variables are given

        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(6.5, 4)) 
        for Np, lab in enumerate(labels):
            time = dataDict[Np]['time']
            axes[0,0].plot(time, dataDict[Np]['Y'], label=lab)
            axes[0,1].plot(time, dataDict[Np]['T'], label=lab)
            axes[1,1].plot(time, dataDict[Np]['omega'], label=lab)
            axes[1,0].plot(time, dataDict[Np]['lambda'], label=lab)
            Elab = np.trapz(dataDict[Np]['E'], dx=1)
            axes[0,2].plot(time, dataDict[Np]['E'], label='{:.0f}'.format(Elab))
            axes[1,2].plot(time, dataDict[Np]['d'], label=lab)

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
        axes[1,1].set_ylabel(r'$\omega$')
        axes[1,0].set_ylabel(r'$\lambda$')
        axes[0,2].set_ylabel(r'$E$ (Gt CO$_2$)')
        axes[1,2].set_ylabel(r'$d$')

    else:
        # an option to have all Gemmes variables
    
        fig, axes = plt.subplots(nrows=5, ncols=4, figsize=(20, 20)) 

        i = -1
        for Np, lab in enumerate(labels):
            dD = dataDict[Np]
            time = dD['time']
            dk = dD.keys()

            for k, key in enumerate(dk):
                if k%4==0:
                    i += 1
                    j = 0
                ax = axes[i,j]
                ax.plot(time, dD[key], label=lab)
                ax.set_xlabel('$t$ (year)')
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.grid(which='both', linewidth=.0001, color='silver',
                    alpha=0.4)
                ax.set_ylabel(key)
                j += 1

    if not returnfig:
        # legend
        axes[0,0].legend()
        plt.tight_layout(pad=0)
        print(namef)
        plt.savefig(namef)
        plt.close(fig)
        return None, None

    else:
        return fig, axes


def conv_for_fortran(name='gemmesR.out', path='.', y_ini=2016,
        y_stop=2100):
    """
    This function can be used to create a gemmes.out file that is in
    the same format as the fortran code.
    labels = ['Y', 'N', 'NG', 'lambda', 'omega', 'd', 'Eind', 'Eland',
       'n', 'sigma', 'Gsigma', 'Pbs', 'p', 'a', 'K', 'T', 'E', 'time']

    VARS = ['time', 'capital', 'npop', 'debt', 'wage', 'productivity',
       'price', 'eland', 'sigma', 'gsigma', 'co2at', 'co2up', 'co2lo',
       'temp', 'temp0', 'pbs', 'pcar', 'omega', 'lambda', 'debtratio',
       'gdp0', 'gdp', 'eind', 'inflation', 'abat', 'n_red_fac',
       'smallpi', 'smallpi_k', 'dam', 'dam_k', 'dam_y', 'fexo', 'find',
       'rcb']
    """
    # load data
    dD = load_Gemmes_data(path=path, y_ini=y_ini, y_stop=y_stop)
    time = dD['time']

    # vars 
    VARS = ['time', 'capital', 'npop', 'debt', 'wage', 'productivity',
       'price', 'eland', 'sigma', 'gsigma', 'co2at', 'co2up', 'co2lo',
       'temp', 'temp0', 'pbs', 'pcar', 'omega', 'lambda', 'debtratio',
       'gdp0', 'gdp', 'eind', 'inflation', 'abat', 'n_red_fac',
       'smallpi', 'smallpi_k', 'dam', 'dam_k', 'dam_y', 'fexo', 'find',
       'rcb']

    # numpy array
    data = np.zeros(shape=(time.size, len(VARS)))
    data[:,0] = time
    data[:,1] = dD['K']
    data[:,2] = dD['N']
    data[:,3] = dD['d']*dD['p']*dD['Y']
    L = dD['lambda']*dD['N']
    data[:,4] = dD['omega']*dD['p']*dD['Y']/L
    data[:,5] = dD['a']
    p = dD['p']
    data[:,6] = p
    data[:,7] = dD['Eland']
    data[:,8] = dD['sigma']
    data[:,9] = dD['Gsigma']
    data[:,13] = dD['T']
    data[:,15] = dD['Pbs']
    theta = 2.6
    print('>  theta={:.2f}'.format(theta))
    data[:,16] = dD['n']**(theta-1)*dD['Pbs']
    data[:,17] = dD['omega']
    data[:,18] = dD['lambda']
    data[:,19] = dD['d']
    data[:,20] = dD['a']*L
    data[:,21] = dD['Y']
    data[:,22] = dD['Eind']
    dotp = np.gradient(p)
    data[:,23] = dotp/p
    data[:,24] = dD['sigma']*dD['Pbs']*dD['n']**theta/theta/1000
    data[:,25] = dD['n']
    pi1, pi2, pi3, zeta3 = 0., 0.00236, 0.0000819, 6.754
    print('>  pi1={:.2e}, pi2={:.2e}, pi3={:.2e}, zeta3={:.2e}'.format(
        pi1,pi2,pi3,zeta3))
    T = dD['T']
    D = 1 - 1./(1 + pi1*T + pi2*T**2 + pi3*T**zeta3)
    data[:,28] = D
    fk = 1./3
    print('>  fk={:.2f}'.format(fk))
    Dk = fk*D
    data[:,29] = Dk
    Dy = 1 - (1 - D)/(1 - Dk)
    data[:,30] = Dy
    rstar = 0.03
    print('>  rstar={:.2f}'.format(rstar))
    data[:,33] = rstar*np.ones(time.size)
    delta, nu = 0.04, 2.7
    print('>  delta={:.2f}  nu={:.2f}'.format(delta, nu))
    deltad = delta + Dk
    pi = 1 - dD['omega'] - rstar*dD['d'] - ((dD['sigma']*data[:,24] + \
        deltad*nu)/(1 - Dy))
    data[:,26] = pi
    data[:,27] = pi*dD['Y']/dD['K']

    header = ''
    for var in VARS:
        header += var+' '
    np.savetxt(fname=name, X=data, header=header)
    print('File {} written.'.format(name))
