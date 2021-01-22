"""
Functions for post-processing the globals data like mean T.
"""

### imports
from loveclim.loveclim import np, plt

### Constants
Nbd = 360 # number of days in one year

### Functions
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




