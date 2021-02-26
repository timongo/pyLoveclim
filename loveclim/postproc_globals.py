"""
Functions for post-processing the globals data like mean T.
"""

### imports
from loveclim.loveclim import np, ReadGlobals
from loveclim.postp_Gemmes import plt
from scipy.ndimage import gaussian_filter1d as gf1d

### Constants
Nbd = 360 # number of days in one year

### Functions
# average_yearly_T
def average_yearly_T(T, ystart=1):
    """
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




