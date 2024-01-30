"""                                                                             
file name: postp_Gemmes.py
language: Python3.6
date: 31/01/2021
author: Hugo Martin
email: martin.hugo@live.com
This module contains functions for simply post processing Gemmes.
"""

### imports -------------------------------------------------------------------
from cycler import cycler
from loveclim.gui import np, os, plt
import sys
from matplotlib.ticker import AutoMinorLocator
import matplotlib as mpl
import pandas as pd
plt.rc('text', usetex=True)
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
mpl.rcParams['lines.linewidth'] = 0.7
mpl.rcParams['lines.markersize'] = 1.5

# global constants -----------------------------------------------------
VARS = ['all', 'capital', 'npop', 'debt', 'wage', 'productivity', 'price',
        'eland', 'sigma', 'gsigma', 'co2at', 'co2up', 'co2lo', 'temp', 
        'temp0', 'pbs', 'pcar', 'rb', 'omega', 'lambda', 'debtratio', 'gdp0', 'gdp',
        'eind', 'inflation', 'abat', 'n_red_fac', 'smallpi', 'smallpi_k',
        'dam', 'dam_k', 'dam_y', 'fexo', 'find', 'rcb', 'g0', 'investment', 'pir']
REORDER = [
    0,  #capital k
    1,  #npop k
    2,  #debt k
    3,  #wage k
    20, #gdp0 r
    21, #gdp k -> go with gdp0
    4,  #productivity k
    18, #lambda r
    17, #omega k
    19, #debtratio r
    5,  #price r
    23, #inflation k
    16, #rb k
    33, #rcb k
    26, #smallpi k
    27, #smallpi_k k
    14, #pbs k
    15, #pcar k
    25, #n_red_fac k
    24, #abat r
    7,  #sigma r
    8,  #gsigma k
    35, #investment
    36, #pir
    22, #eind k
    6,  #eind+eland k
    9,  #co2at k
    10, #co2up k
    11, #co2lo k
    12, #temp k
    13, #temp0 r
    28, #dam k
    29, #dam_k k
    30, #dam_y k
    34, #g0 -> replace gdp
]
LABELS = ['nolab', 'K', 'N', 'D', 'W', 'a', 'p', 'E_{\mathrm{ind}} + E_{\mathrm{land}}', '\sigma',
        'g_{\sigma}', 'CO_2^{\mathrm{at}}', 'CO_2^{\mathrm{up}}', 'CO_2{\mathrm{lo}}',
        'T', 'T_0', 'p_{\mathrm{bs}}', 'p_\mathrm{c}', 'r_b', '\omega', '\lambda', '(1-A)(1-\mathrm{D}^Y)d/\\nu', 'Y^0',
        'Y', 'E_{\mathrm{ind}}', 'i', 'A', 'n', '\pi', '\pi_K', '\mathrm{D}', '\mathrm{D}^K', '\mathrm{D}^Y',
        'f_{\mathrm{exo}}', 'f_{\mathrm{ind}}', 'RCB', 'g = \dot{Y}/Y', 'I', '\Pi_r']
LABELS = [r'$'+lab+'$' for lab in LABELS]
YCOLORS = ['k']*len(LABELS)
for i in [5, 7, 13, 18, 19, 20, 24]:
    YCOLORS[i]='r'
ARGD = {'name':'default', 'doption':'all', 'col':0, 'plot':False,
    'compf':False, 'tf':10000000, 'fname':'gemmes.out',
    'multi':False}
SSP85 = os.path.join('SSP_CMIP6_emissions', 'ssp58.5_emissions_Gtco2.dat')
SSP60 = os.path.join('SSP_CMIP6_emissions', 'ssp46.0_emissions_Gtco2.dat')
SSP70 = os.path.join('SSP_CMIP6_emissions', 'ssp37.0_emissions_Gtco2.dat')
FIGSIZE = (20, 20)
SIZE = 8
sta, sto = 2095, 2105
# for plot
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
linestyle = [
     'solid',
     #(0, (1, 10)),
     #(0, (1, 1)),
     'solid',
     (0, (5, 10)),
     (0, (5, 5)),
     (0, (5, 1)),
     (0, (3, 10, 1, 10)),
     (0, (3, 5, 1, 5)),
     (0, (3, 1, 1, 1)),
     (0, (3, 5, 1, 5, 1, 5)),
     (0, (3, 10, 1, 10, 1, 10))]
     #(0, (3, 1, 1, 1, 1, 1))]
default_cycler = (cycler(color=colors) + cycler(linestyle=linestyle))
plt.rc('lines', linewidth=2)
plt.rc('axes', prop_cycle=default_cycler)

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


# ----------------------------------------------------------------------
# for post prosess -----------------------------------------------------
# ----------------------------------------------------------------------
DICVAR = {var:i+1 for i, var in enumerate(VARS[1:])}
DICVAR['time'] = 0
def print_VARS_old():
    print('\n0', 'time')
    for i, var in enumerate(VARS[1:]):
        print(i+1, var, YCOLORS[i])
    print()

def print_VARS():
    for key, value in DICVAR.items():
        print(value, key)
    print()

def load_fields(path, fields, step=1):
    """
    Used for loading a set of fields:

    INPUTS:
        path : string
            path of the gemmes.out file
        fields : list (string)
            name of fields
        step : integer
            step between dates

    OUPUTS:
        dic_of_data : dictionary
            Dictionary of data
    """

    fulldata = np.loadtxt(path, skiprows=1)[::step,:]
    dic_of_data = {}
    for field in fields:
        try:
            dic_of_data[field] = fulldata[:,DICVAR[field]]
        except KeyError:
            raise KeyError(
                "Field '{}' does not exist in gemmes.out".format(field)
            )

    return dic_of_data

def load_fields_DataFrame(path, fields, step=1):
    """
    Used for loading a set of fields as a DataFrame.

    INPUTS:
        path : string
            path of the gemmes.out file
        fields : list (string)
            name of fields
        step : integer
            step between dates

    OUPUTS:
        data : DataFrame
            DataFrame of data
    """
    
    tmp = load_fields(path, fields, step)
    df = pd.DataFrame(tmp)
    
    return df


def load_data(args):
    """
    Fonction to load data from gemmes.
    """
    if not args['multi']:
        fname = args['fname']
        data1 = loadf(fname)
        if args['compf']!=False:
            data2 = loadf(args['compf'])
            datas = [data1, data2]
        else:
            datas = [data1]
    else:
        datas = []
        for fname in args['multi']:
            tmpdata = loadf(fname)
            datas.append(tmpdata)

    print('\n> Gemmes data loaded')
    return datas

def loadf(fname):
    try:
        print(fname, 1)
        data = np.loadtxt(fname=fname, skiprows=1)
    except ValueError:
        print(fname, 2)
        data = np.genfromtxt(fname=fname, skip_header=1, filling_values=1E99)

    return data

def moy(time, dat, sta, sto):                                            
    select = (time>=sta) * (time<=sto)                                   
    mean = np.mean(dat[select])                                          
    return mean

def fusage():
    """
    Function Usage
    """
    print('\n> Usage: python post-process.py [-h] [-name figure name] [-plot]')
    print('> [-fname file name] [-compf] [-tf final time] [-doption option]')
    print('>\n> Options:')
    print('>     -h :: display this message.')
    print('>     -name : string :: name of figure.')
    print('>     -fname :: string :: name of .out file')
    print('>     -compf : string :: name of the file *.out to compare with.')
    print('>     -plot :: to show the figures instead to save them.')
    print('>     -tf :: int :: final time for post-process.')
    print('>     -doption : string :: name of the considerd option in:', VARS)
    print('>     -multi : array string :: all the names that folow multi are post-processed.')
    print()
    ferror()

def change_mpl_opt(opt):
    """
    This function changes matplotlib options.
    """
    size = SIZE
    figsize = FIGSIZE
    if opt!='all':
        size = 12
        figsize = (5, 5)

    SMALL_SIZE = SIZE 
    MEDIUM_SIZE = SIZE 
    BIGGER_SIZE = SIZE 
    plt.rc('font', size = SMALL_SIZE) # controls default text sizes
    plt.rc('axes', titlesize = BIGGER_SIZE) # fontsize of the axes title
    plt.rc('axes', labelsize = MEDIUM_SIZE) # fontsize of the x and y labels
    plt.rc('xtick', labelsize = SMALL_SIZE) # fontsize of the tick labels
    plt.rc('ytick', labelsize = SMALL_SIZE) # fontsize of the tick labels
    plt.rc('legend', fontsize = MEDIUM_SIZE) # legend fontsize
    plt.rc('figure', titlesize = BIGGER_SIZE) # fontsize of the figure title
    mpl.rcParams['lines.linewidth'] = 0.5
    mpl.rcParams['lines.markersize'] = 1.5

    return figsize

def special_plot(axes, ax, var, time, datp, ts, reo, fname):
    """
    Function that plots for real.
    """
    if var=='eind+eland':
        Etot = (datp[:ts,6]+datp[:ts,22])
        ax.plot(time, Etot, label=fname)
        # ssp5-8.5
        datssp85 = np.loadtxt(SSP85)
        datssp60 = np.loadtxt(SSP60)
        datssp70 = np.loadtxt(SSP70)

        btinf = max(time[0], datssp85[0,0])
        btup = min(time[-1], datssp85[-1,0])
        tselect = (time>=btinf) * (time<=btup)
        fssp85 = np.interp(time[tselect], datssp85[:,0], datssp85[:,1])
        fssp60 = np.interp(time[tselect], datssp60[:,0], datssp60[:,1])
        fssp70 = np.interp(time[tselect], datssp70[:,0], datssp70[:,1])

        ax.plot(time[tselect], fssp85, color='k', linestyle='--')
        ax.plot(time[tselect], fssp60, color='y', linestyle='--')
        ax.plot(time[tselect], fssp70, color='r', linestyle='--')

        ax.text(x=2020, y=150, s='SSP5-8.5', color='k', fontsize='medium')
        ax.text(x=2020, y=130, s='SSP5-6.0', color='y', fontsize='medium')
        ax.text(x=2020, y=140, s='SSP5-7.0', color='r', fontsize='medium')
        Eeco = fssp70 - Etot[tselect]
        
        return tselect, Eeco, Etot

    elif var=='co2at':
        ax.plot(time, datp[:ts,reo]/2.124, label=fname)

    elif var=='pbs':
        ax.plot(time, datp[:ts,reo], label=fname)
        print('pbs/pbs0 temps final: {:.2f} \%'.format(100*datp[-1,reo]/datp[0,reo]))
    elif var=='sigma':
        ax.plot(time, datp[:ts,reo], label=fname)
        print('sigma/sigma0 temps final: {:.2f} \%'.format(100*datp[-1,reo]/datp[0,reo]))
    elif var=='gdp':
       axes[0,4].plot(time, datp[:ts,reo], linestyle='--', 
                      color='C0', label=fname)
    elif var=='debtratio':
        tc = (1. - datp[:ts,24])*(1. - datp[:ts,30])
        ax.plot(time, tc*datp[:ts,reo]/3.0, label=fname)

    elif var=='investment':
        # ici consommation
        conso = 1. - datp[:ts, reo] / datp[:ts,5] / datp[:ts,21]
        ax.plot(time, conso, label=fname)

    elif var=='pir':
        # ici epargne des ménages
        div = datp[:ts,26]*datp[:ts,5]*datp[:ts,21] - datp[:ts, reo]
        wL = datp[:ts,3] * datp[:ts,18] * datp[:ts,21]
        conso = datp[:ts,5] * datp[:ts,21] - datp[:ts, reo]
        savings = div + wL - conso
        savings = savings / datp[:ts,5] / datp[:ts,21]
        ax.plot(time, savings, label=fname)

    else:
        ax.plot(time, datp[:ts,reo], label=fname)

    return 0,0,0

def draw_figure(args, datas):
    """
    Draw the figure and save it.
    """
    if not args['multi']==False:
        multi_fig(args, datas)
        return 0

    # initialisation
    opt = args['doption']
    col = args['col']
    name = args['name']
    compf = args['compf']
    boolcompf = compf!=False
    tstop = args['tf']
    time2, datp2, t2s =  None, None, None

    if boolcompf: # there are two *.out files to compare

        if name=='default':
            name = 'gemmes_comp.svg'

        opt = 'all'
        nbs = 35
        nbr, nbc = nbs//5, nbs//7

        datp = datas[0]
        time = datp[:,0]
        time = time[time<=tstop]
        ts = time.size

        datp = datp[:,1:]
        datp2 = datas[1]
        time2 = datp2[:,0]
        time2 = time2[time2<=tstop]
        t2s = time2.size

        datp2 = datp2[:t2s,1:]
        Vars = VARS[1:]

    else: # there is only one *.out file

        if name=='default':
            name = 'gemmes_'+opt+'.pdf'

        data = datas[0]
        if opt=='all':
            nbs = 35
            nbr, nbc = nbs//5, nbs//7
            datp = data[:,1:]
            Vars = VARS[1:]
        else:
            nbr, nbc = 1, 1
            datp = data[:,col:col+1]
            Vars = VARS[col:col+1]
        time = data[:,0]
        time = time[time<=tstop]
        ts = time.size
        datp = datp[:ts,:]

    # make figure
    figsize = change_mpl_opt(opt)
    fig, axes = plt.subplots(nrows=nbr, ncols=nbc, sharex=True, figsize=figsize)
    j = -1
    if len(Vars)==1:
        axes = np.asarray([[axes]])
    for n, reo in enumerate(REORDER):
        var = Vars[reo]
        N = n+1 if (n>=15) else n
        if N%5==0 or (N==16):
            j += 1
            i = 0
        if var=='g0':
            ax = axes[1,0]
            #print('Mean prod=', np.mean(0.009+0.65*datp[:ts,reo]))
        else:
            ax = axes[j, i]
        try:
            # plot and second plot
            if var=='eind+eland':
                tselect, Eeco, Etot = special_plot(axes, ax, var, time, datp, ts, reo, args['fname']) 
            else:
                special_plot(axes, ax, var, time, datp, ts, reo, args['fname'])

            # add stuff
            additional_infos(var, reo, time2, datp2, t2s, ax, axes, boolcompf,
                args['compf'], time)

            # labels and so on
            if j==(nbr-1):
                ax.set_xlabel('t')
            if var=='gdp0':
                ax.set_ylabel(r'$Y^0$ (plain)          $Y$ (dotted)')
            elif var=='gdp':
                ax.set_ylabel(r'GDP growth')
            elif var=='pir':
                ax.set_ylabel('Cap. savings to GDP')
            elif var=='investment':
                ax.set_ylabel('Consumption to GDP')

            else:
                ax.set_ylabel(LABELS[reo+1]) # var
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            #ax.grid(which='both', linewidth=.0001, color='silver', alpha=.4)
            if YCOLORS[reo+1]!='k':
                ycolor=YCOLORS[reo+1]
                ax.tick_params(axis='y', colors=ycolor)
                ax.yaxis.label.set_color(ycolor)
                ax.spines['bottom'].set_color(ycolor)
                ax.spines['left'].set_color(ycolor)
                ax.spines['right'].set_color(ycolor)
                ax.spines['top'].set_color(ycolor)

        except IndexError:
            continue
        i += 1

    # transfer = Sfc - Tfc
    sa = 0. #np.float64(input("input value of sa: "))
    conv10to15=1.160723971/1000. # conversion factor
    tfc = datp[:,15]*conv10to15*datp[:,22]
    aY0 = datp[:,24]*datp[:,20]

    cout_publique = sa*aY0 - tfc
    cout_firmes = aY0 - cout_publique

    axes[3,4].set_ylabel(r'mitigation costs (normalized by $Y$)')
    axes[3,4].plot(time, -cout_firmes[:ts]/datp[:ts,21], color='C0', linestyle=':')
    axes[3,4].plot(time, -cout_publique[:ts]/datp[:ts,21], color='C0', linestyle='--')
    axes[3,4].plot(time, -tfc[:ts]/datp[:ts,21], color='C0', linestyle='-')
    axes[3,4].text(x=2070, y=-0.07, s='-- public', color='C0', fontsize='medium')
    axes[3,4].text(x=2070, y=-0.08, s='.. private', color='C0', fontsize='medium')
    axes[3,4].text(x=2070, y=-0.09, s='__ carbon tax', color='C0', fontsize='medium')
    axes[3,4].set_ylim(ymin=-0.1, ymax=0.05)
 
    dt = time[-1]-time[-2]
    select = (time <= 2050)
    ntselect = (time[select] <= 2050)
    tot_abatprice = -dt*np.sum(aY0[select])
    tot_cf = -dt*np.sum(cout_firmes[select])
    tot_cp = -dt*np.sum(cout_publique[select])

    #couevite = conv10to15*dt*np.sum(datp[select,15]*Eeco[ntselect])
    #tot_Eevite = dt*np.sum(Eeco[ntselect])
    #tot_E = dt*np.sum(Etot[select])

    print('\nPoids de la transition de 2016 a 2051: ')
    print(' -- abatement total = {:.2f} 10^12 $'.format(tot_abatprice)) 
    print(' -- cout sur le secteur prive = {:.2f} 10^12 $'.format(tot_cf))
    print(' -- cout sur le secteur publique = {:.2f} 10^12 $'.format(tot_cp))

    prixcar(time, datp[:,15], datp[:,14])

    #print('\nCout evite: int pc * (Essp85 - Etot) = {:.2f} 10^12 $'.format(couevite))
    #print('Total Emission evitees int Eeco = {:.2f} GtCO2-e'.format(tot_Eevite))
    #print('Total emissions Etot= {:.2f} GtCO2-e'.format(tot_E))

    select = (time <= 2100)
    ntselect = (time[select] <= 2100)
    tot_abatprice = -dt*np.sum(aY0[select])
    tot_cf = -dt*np.sum(cout_firmes[select])
    tot_cp = -dt*np.sum(cout_publique[select])

    #couevite = conv10to15*dt*np.sum(datp[select,15]*Eeco[ntselect])
    #tot_Eevite = dt*np.sum(Eeco[ntselect])
    #tot_E = dt*np.sum(Etot[select])

    print('')
    print('\nPoids de la transition de 2016 a 2100: ')
    print(' -- abatement total = {:.2f} 10^12 $'.format(tot_abatprice)) 
    print(' -- cout sur le secteur prive = {:.2f} 10^12 $'.format(tot_cf))
    print(' -- cout sur le secteur publique = {:.2f} 10^12 $'.format(tot_cp))

  
    if boolcompf:
        axes[0,0].legend()

    if args['plot']:
        plt.show()
    else:
        print('\n> Writing file {}.'.format(name))
        plt.tight_layout()
        plt.savefig(fname=name)
    plt.close(fig)

    return 0

def prixcar(time, pcar, pbs):
    print('\nMoyenne du prix du carbone:')
    t = np.arange(2015, 2101, 5)
    select = (time>=2015)*(time<=2050)
    print(' -- mean pbs de 2016 a 2050 = {:.3f}'.format(np.mean(pbs[select])))
    for n, tt in enumerate(t[:-1]):
        select = (time>tt)*(time<=t[n+1])
        print('-- de {:d} a {:d}: {:.2f} $/tCO2'.format(tt,t[n+1],
            np.mean(pcar[select])))

def additional_infos(var, reo, time, datp, ts, ax, axes, boolcompf, namecompf,
                     time_main=None):
    """
    The function additional_infos can be used to add informations to
    standard plots.
    Inputs:
        n : int : number of the ax
        time : array float : time 2nd simu
        datp : array float : data 2nd simu
        ts : int : time indice for 2nd simu
        ax : matplotlib ax : ax of the plot
        boolcompf : bool : True whether we compare two simus
        time_main : array float : time of 1fst simu
    """
    # data comparison

    if boolcompf:
        special_plot(axes, ax, var, time, datp, ts, reo, namecompf) #args['fname'])

    # atmospheric co2
    if var=='co2at':
        preind_co2 = 588./2.124 # ppm
        Dxco2 = 2.*preind_co2
        ax.plot([time_main[0], time_main[-1]], [Dxco2, Dxco2], color='red') 
        ax.text(x=2070, y=500, s='2 * CO2 pre-ind.', color='red')

    # total emissions
    elif var=='eind+eland':
        ax.set_xlim(xmin=2015, xmax=time_main.max())
        #ax.set_ylim(ymin=0, ymax=170)

    elif var=='debtratio':
       ax.set_ylim(-1.3, 1.3)
       ax.plot([time_main[0], time_main[-1]], [1., 1.], '--', color='r')

    #elif var=='g0':
        #ax.set_ylim(-0.2, 0.1)

    elif var=='temp':
        ax.plot([2095, 2105], [3.6, 3.6], color='red', linestyle='--')
        ax.text(x=2107, y=3.6, s='SSP3-7.0', color='red')

    #elif var=='capital' or var=='debt' or var=='wage' or var=='gdp0':
    #    ax.set_ylim(ymin=0.)

    #elif var=='pbs':
    #    ax.set_ylim(0., 1.)

    #elif var=='n_red_fac':
    #    ax.set_ylim(0, 1.1)

    #elif var=='lambda' or var=='omega' or var=='n_red_fac':
    #    ax.set_ylim(0, 1.1)

    #elif var=='sigma':
    #    ax.set_ylim(0, 1.2)

    elif var=='pcar':
        ax.plot([2030, 2030], [0., 100], color='r', linestyle='--')
        ax.plot([2050, 2050], [0., 400], color='r', linestyle='--')

    #elif var=='smallpi' or var=='smallpi_k':
    #    ax.set_ylim(-0.5, 1)

    #elif var=='abat':
    #    ax.set_ylim(0, 0.2)

    #elif var=='inflation' or var=='rb' or var=='rcb':
    #    ax.set_ylim(ymin=-0.3, ymax=0.2)

def makedatadic(data):
    """
    Turn np array into dict.
    """
    datadic = {}
    for n, var in enumerate(VARS):
        datadic[var] = data[:,n]

    datadic['time'] = data[:,0]

    return datadic

def ferror(key=0, msg=''):
    """
    Error function.
    """
    if key==0:
        print('\n> Post process run with sucess.\n')
    else:
        print('\n> Error {:d}'.format(key))
        if key==1:
            print('>     Need arguments in option {}, try -h'.format(msg))

        elif key==2:
            print('>     doption must be in the list, try -h.')

        print('\n> Post process run with an error.\n')

    # exit script
    print()
    sys.exit()

### Functions for multiple simulations
# multi_fig
def multi_fig(args, datas):
    """
    Do multiple post process and plot the main economical variables in 
    a same figure to compare them.
    """
    snames = args['multi']
    tstop = args['tf']

    # initialize figure
    fig, axes = plt.subplots(nrows=7, ncols=2, sharex=True, figsize=(15,10))

    # loop on simulations
    for n, sname in enumerate(snames):
        tmpdata = datas[n]
        time = tmpdata[:,0]
        time = time[time<=tstop]
        ts = time.size
        tmpdata = tmpdata[:ts,:]
        zero = sub_multi_fig(axes, time, tmpdata, sname)

    for ax in axes:
        ax[0].grid(linewidth=0.001)
        ax[1].grid(linewidth=0.001)

    axes[6,1].legend(loc='upper center',
       bbox_to_anchor=(0.5,-0.5), ncol=2)
    
    # borned to ^100
    axes[5,0].set_ylim(ymin=1., ymax=5.)

    if args['plot']:
        plt.show()
    else:
        print('\n> Writing file {}.'.format(args['name']))
        plt.tight_layout()
        plt.savefig(fname=args['name'])

    plt.close(fig)

def sub_multi_fig(axes, time, data, sname):
    """
    Where we plot the variable we want for multi.
    """
    N, ax = 8, axes[0,0] # sigma
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 9, axes[0,1] # gsigma
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 16, axes[1,0] # pc
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 28, axes[1,1] # D
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
       
    N, ax = 22, axes[2,0] # Eind
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 10, axes[2,1] # co2at
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 13, axes[3,0] # temp
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
       
    N, ax = 21, axes[3,1] # Y
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 18, axes[4,0] # lambda
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 17, axes[4,1] # omega
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 19, axes[5,0] # d
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 26, axes[5,1] # pi
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    N, ax = 23, axes[6,0] # i
    ax.plot(time, data[:,N])
    ax.set_ylabel(LABELS[N])
    ax.set_xlabel('time (years)')    
    N, ax = 33, axes[6,1] # rcb
    ax.plot(time, data[:,N], label=sname[8:-33])
    ax.set_ylabel(LABELS[N])
    ax.set_xlabel('time (years)')    

    return 0

# functions ------------------------------------------------------------
def load_args():
    """
    Function that loads data
    """
    argl = sys.argv

    # try to find help
    for arg in argl:
        if arg=='-help':
            fusage()

    # find arguments
    argdict = ARGD
    for n, arg in enumerate(argl):
        if arg[0]=='-':
            Arg = arg[1:]
            if Arg=='name':
                try:
                    argdict[Arg] = argl[n+1]
                except IndexError:
                    ferror(key=1, msg=Arg)
            elif Arg=='plot':
                argdict['plot'] = True
            elif Arg=='compf':
                try:
                    argdict['compf'] = argl[n+1]
                except IndexError:
                    ferror(key=1, msg=Arg)
            elif Arg=='tf':
                try:
                    argdict['tf'] = int(argl[n+1])
                except IndexError:
                    ferror(key=1, msg=Arg)
            elif Arg=='fname':
                try:
                    argdict['fname'] = argl[n+1]
                except IndexError:
                    ferror(key=1, msg=Arg)
            elif Arg=='path':
                try:
                    argdict['path'] = argl[n+1]
                except IndexError:
                    ferror(key=1, msg=Arg)
            elif Arg=='multi':
                multi = []
                for m in range(n+1,1000):
                    try:
                        temparg = argl[m]
                        if temparg[0]=='-':
                            break
                        else:
                            if os.path.exists(temparg):
                                multi.append(temparg)
                            else:
                                print(temparg, 'not exist > ignored')
                    except IndexError:
                        if len(multi)>0:
                            break
                        else:
                            ferror(key=1, msg=Arg)
                if len(multi)>0:
                    argdict['multi'] = multi 
            elif Arg=='doption':
                try:
                    tmp = argl[n+1]
                    varfound = False 
                    for i, var in enumerate(VARS):
                        if tmp==var:
                            varfound = True
                            argdict['col'] = i
                            break
                    argdict['doption'] = tmp if varfound else ferror(key=2)
                except IndexError:
                    ferror(key=1, msg=Arg)

    print('\n> Arguments loaded')
    return argdict

# print_costs
def print_costs(data, sa=0.5, nu=3.):
    """
    This function prints costs of transition.
    """
    conv10to15=1.160723971/1000. # conversion factor

    time = data[:,0]
    dt = time[-1]-time[-2]
    pcar = data[:,16]
    Eind = data[:,23]
    A = data[:,25]
    gdp0 = data[:,21]

    tfc = pcar*conv10to15*Eind
    aY0 = A*gdp0
    cout_publique = sa*aY0 - tfc
    cout_firmes = aY0 - cout_publique

    print(' - carbon price evolution:')
    time_int = np.arange(2017.5, 2101, 5)
    for n, t in enumerate(time_int[:-1]):
        tt = time_int[n+1]
        select = (time >= t)*(time <= tt)
        pcarmean = np.mean(pcar[select])
        print('-- [{:.1f}, {:.1f}]: {:.2f} $/tCO2'.format(t, tt,
            pcarmean))

    print(' - mitigation costs:')
    for n, t in enumerate([2050, 2100]):
        select = (time <= t)
        tot_abatprice = -dt*np.sum(aY0[select])
        tot_cf = -dt*np.sum(cout_firmes[select])
        tot_cp = -dt*np.sum(cout_publique[select])
    
        print('- from 2016 to {:d}:'.format(t))
        print(' -- total abatment = {:.2f} 10^12 $'.format(
            tot_abatprice)) 
        print(' -- private cost = {:.2f} 10^12 $'.format(
            tot_cf))
        print(' -- public cost = {:.2f} 10^12 $\n'.format(tot_cp))

def final_state(data, nu=3.):
    flambda = data[-1,19]
    fomega = data[-1,18]
    fdebtratio = data[-1,20]/nu
    print(' - final state:')
    print(' -- lambda={:.2f}'.format(flambda)+\
        '\n -- omega={:.2f}'.format(fomega)+\
        '\n -- debtratio={:.2e}\n'.format(fdebtratio))
    collapse = '{:.2f}'.format(flambda)=='0.00'

    return collapse
