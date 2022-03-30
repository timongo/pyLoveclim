"""
Functions for post-processing the globals data like mean T.
"""

### imports
from loveclim.loveclim import np, ReadGlobals
from loveclim.postp_Gemmes import plt
from scipy.ndimage import gaussian_filter1d as gf1d
import os, sys
from shutil import copy2

### Constants
Nbd = 360 # number of days in one year

### Functions
# average_yearly_T
def average_yearly_T(T, ystart=1):
    """
    Deprecated : Use get_yearly_temp_ghg.

    Return the years and the average yearly temperature from a global book file
    from iLOVECLIM.
    
    Parameters:
        T : array like
            T is the vector returned by the function ReadGlobal.
        ystart : int
            The number of the first year computed.
    
    Returns:
        years : array like
            Number of years computed in the iLOVECLIM's simulation.
        Tmoy : array like
            Average temperature vector. Same size than years.
    """
    # find number of years
    Nby = T.size//Nbd
    years = np.arange(start=ystart, stop=ystart+Nby, step=1)
    
    # loop on years
    Tmoy = [] 
    for y in range(Nby): 
        tmp = T[y*Nbd:(y+1)*Nbd] 
        Tmoy.append(np.mean(tmp)) 
    Tmoy = np.asarray(Tmoy)
    
    return years, Tmoy

# get_yearly_temp_ghg 
def get_yearly_temp_ghg(filename, ystart=1, maxrows=0):
    """
    Get yearly temperature and ghg from ipcc*

    Parameters:
        filename : string
            filename is the path of the ipcc file.
        ystart : int
            The number of the first year computed.
        maxrows : int
            Maximun number of rows in file for loading
    
    Returns:
        t : array like
            Number of years computed in the iLOVECLIM's simulation.
        tm : array like
            Average temperature vector. Same size than years.
        co2 : array like
            Average co2 concentration in atmosphere. Same size than years.
    """

    if maxrows>0:
        dat = np.loadtxt(filename, max_rows=maxrows)
    else:
        dat = np.loadtxt(filename)
    t = dat[:,0]
    tm = dat[:,2]
    co2 = dat[:,6]

    return t, tm, co2

# quick_view_T
def quick_view_T(bookname, path='./', ystart=1):
    """
    Function for quick viewing temperature

    Pamareters:
        bookname : string
            name of a book file
        path : string
            path to folder that contains bookname
        ystart : int
            The number of the first year computed.
    """
    # read data
    t,Y,D,T = ReadGlobals(path+bookname)

    # mean
    Ym, tmoy = average_yearly_T(T, ystart=ystart)

    # plot
    fig, ax = plt.subplots()
    ax.plot(t//360+ystart, T, color='C0', alpha=0.5)
    ax.plot(Ym, tmoy, color='C0', label='iloveclim')

    # legend
    namef = path+'quick_view_'+bookname+'.pdf'
    ax.set_xlabel(r't (years)')
    ax.set_ylabel(r'T (°C)')
    plt.tight_layout()
    plt.savefig(namef)
    print(namef)
    plt.close(fig)

def copy_restart(path, oldyear, newyear):
    """
    The function copy_restart can be used to change the starting year 
    of a iloveclem restart set of files.
    Inputs:
        path : string : path where there is the ic****** dir to copy
        oldyear : string of len=6 : old year of restart 
        newyear : string of len=6 : new year of restar
    """
    fnames = ['res_masks*.om', 'res*.om', 'inatphy*.dat',
              'inatdyn*.dat', 'inland*.dat', 'incoup*.dat']
    oldpath = os.getcwd()
    os.chdir(path)
    # tests
    try:
        test = int(oldyear)
        test = int(newyear)
        if (len(oldyear)!=6)or(len(newyear)!=6):
            raise ValueError("Nombre d'ints")
        os.chdir('ic{}'.format(oldyear))
        os.chdir('..')
    except ValueError:
        print('Both oldyear and newyear must be string of 6 integers')
        print('Examples: 003000, 012345, 001750')
        sys.exit()
    except FileNotFoundError:
        os.chdir(oldpath)
        raise FileNotFoundError(
            'Restart with oldyear {} does not exist'.format(oldyear))

    # copy
    os.system('cp -r ic{} ic{}'.format(oldyear, newyear))
    os.chdir('ic{}'.format(newyear))
    for n, f in enumerate(fnames):
        os.rename(f.replace('*', oldyear), f.replace('*', newyear))
        print(f)
    os.chdir(oldpath)
