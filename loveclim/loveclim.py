#!/usr/bin/python

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import plotxtor as px
import netCDF4 as nc

def Figure(x,y,c,
           xh=None,
           yh=None,
           fignumber=None,
           hold=False,
           title='',
           plottype='contourf',
           contourlabel=False,
           colorbar=True,
           cscale=None,
           contours=None,
           ncontours=20,
           zlim=None,
           elev=20,
           azim=60,
           theta=0,
           ih = 'h',
           sign=1.,
           xlab='',
           ylab='',
           clab='',
           equal=True,
           cmap=plt.get_cmap('coolwarm'),
           alpha=0.3,
           xtorlabel=''):
    """
    My 1d/2d/3d plotting routine
    x,y,c are supposed to have the same dimensions (2d array). c represents the data and x,y the spatial grid, of course
    """

    if (c.max()==c.min() and plottype!='1d_cuts'):
        print("Impossible to plot completely flat data")
        if fignumber is not None:
            if plt.fignum_exists(fignumber):
                plt.figure(fignumber)
                plt.clf()
        return

    if len(c.shape)==2:
        # Check if figure exists    
        # Figure
        if fignumber is None:
            fig = plt.figure()
        else:
            fig = plt.figure(fignumber)
        if hold==False:
            plt.clf()

        if c.min()==c.max():
            ncontours=0

        dx = x[1,0]-x[0,0]
        dy = y[0,1]-y[0,0]

        # Adjust plot inner size
        plt.subplots_adjust(bottom=0.11, left=.125, right=0.97, top=.94)

        # CONTOURF
        if plottype=='contourf':
            if cscale==None:
                cscale = [c.min(),c.max()]
            if contours==None:
                im = plt.contourf(x,y,c,ncontours,vmin=cscale[0],vmax=cscale[1],cmap=cmap)
                plt.contour(x,y,c,ncontours,colors='k',alpha=alpha)
            else:
                im = plt.contourf(x,y,c,contours,vmin=cscale[0],vmax=cscale[1],cmap=cmap)
                plt.contour(x,y,c,contours,colors='k',alpha=alpha)
            if equal:
                plt.axis('equal')
            if colorbar:
                cb = fig.colorbar(im)        
                cb.formatter.set_powerlimits((0, 0))
                cb.update_ticks()
        # PCOLOR
        if plottype=='pcolor':
            ax = plt.gca()
            if cscale==None:
                im = plt.pcolor(x,y,c,cmap=cmap)
            else:
                im = plt.pcolor(x,y,c,vmin=cscale[0],vmax=cscale[1],cmap=cmap)
            ax.images.append(im)
            ax.set_title(title)
            if equal:
                plt.axis('equal')
            if colorbar:
                cb = fig.colorbar(im)
                cb.formatter.set_powerlimits((0, 0))
                cb.update_ticks()
        # CONTOUR
        elif plottype=='contour':
            if contours is None:
                im=plt.contour(x,y,c,np.linspace(np.min(c),np.max(c),ncontours),cmap=cmap)
            else:
                im=plt.contour(x,y,c,contours,linewidths=1.5,cmap=cmap)
            ax=plt.gca()
            ax.set_title(title)
            if equal:
                plt.axis('equal')
            if contourlabel:
                plt.clabel(im,inline=1,fontsize=10)
        # 3D WATERFALL
        elif plottype=='3d_waterfall':
            ax=fig.gca(projection='3d')
            if cmap==None:
                cmap=plt.get_cmap('jet')
            im=ax.plot_surface(x,y,c,rstride=1,cstride=1,cmap=cmap,linewidth=0,antialiased=False)
            ax.set_title(title)
            if zlim!=None:
                ax.set_zlim(zlim)
            ax.view_init(elev=elev,azim=azim)
        # 3D CONTOUR
        elif plottype=='3d_contour':
            ax=fig.gca(projection='3d')
            if contours==None:
                im=ax.contour(x,y,c,np.linspace(np.min(c),np.max(c),ncontours))
            else:
                im=ax.contour(x,y,c,contours)
            ax.set_title(title)
            if zlim!=None:
                ax.set_zlim(zlim)
            ax.view_init(elev=elev,azim=azim)
        # 1D CUTS
        elif plottype=='1d_cuts':
            lmax = c.shape[0]-1
            mmax = c.shape[1]-1
            theta = np.mod(theta,2.*np.pi)
            it1 = int(rx.cont_2_ind(theta,np.linspace(0,2.*np.pi,mmax+1)))
            it2 = int(np.mod(it1 + mmax/2,mmax))
            # Build data
            if ih=='h':
                R0 = 0.5*(x[0,0] + x[1,0])
                x = x-R0
                r = np.zeros(2*lmax)
                r[lmax:] = np.sqrt(x[1:lmax+1,it1]**2 + y[1:lmax+1,it1]**2)
                r[:lmax] = -np.sqrt(x[lmax+1:0:-1,it2]**2 + y[lmax+1:0:-1,it2]**2)
                data = np.zeros(2*lmax)
                data[lmax:] = c[1:lmax+1,it1]
                data[:lmax] = sign*c[lmax+1:0:-1,it2]
            elif ih=='i':
                R0 = x[1,0]
                x = x-R0
                r = np.zeros(2*lmax)
                r[lmax:] = np.sqrt(x[1:lmax+1,it1]**2 + y[1:lmax+1,it1]**2)
                r[:lmax] = -np.sqrt(x[lmax+1:0:-1,it2]**2 + y[lmax+1:0:-1,it2]**2)
                data = np.zeros(2*lmax)
                data[lmax:] = c[1:lmax+1,it1]
                data[:lmax] = sign*c[lmax+1:0:-1,it2]
            # Plot
            Plot(r,data,'rs-',lw=2,ms=5,markeredgewidth=0)
        if plottype!='1d_cuts':
            plt.axis([np.min(x)-dx, np.max(x)+dx, np.min(y)-dy,np.max(y)+dy])
            plt.xlabel(xlab,fontsize=15,ha='left')
            plt.ylabel(ylab,fontsize=15,va='center',ha='center')
        else:
            absc = np.abs(c)
            minval = -max(absc.min(),absc.max())
            maxval = max(absc.min(),absc.max())
            if minval==maxval:
                minval = minval-0.5
                maxval = maxval+0.5
            plt.axis([np.min(r), np.max(r),minval,maxval])
            plt.xlabel(xlab,fontsize=14)
            plt.ylabel(clab+ylab,fontsize=14)
    elif len(c.shape)==3 and c.shape[2]==3:
        if (xh is None or yh is None):
            print("You need bigr,bigz AND bigrh,bigzh")
            return
        # sign must be flipped for negative radius in r and theta components
        sign = [-1.,-1.,1]
        compname=["$_r$","$_{\Theta}$","$_\phi$"]
        # if a vector, ih represents the first gris, so
        # ih = 'i' means ihh
        # ih = 'h' means hii
        if ih=='i':
            ihh = ['i','h','h']
            X = [x,xh,xh]
            Y = [y,yh,yh]
        elif ih=='h':
            ihh = ['h','i','i']
            X = [xh,x,x]
            Y = [yh,y,y]
        left,right,wspace=0.07,0.97,0.3;
        plt.figure(fignumber,figsize=(18,6));
        plt.clf()
        gs = matplotlib.gridspec.GridSpec(1,3);
        gs.update(wspace=wspace,left=left,right=right);
        for ncomp in range(3):
            ax=plt.subplot(gs[ncomp])
            Figure(X[ncomp],Y[ncomp],c[:,:,ncomp],
                   fignumber=fignumber,
                   hold=True,
                   title=title,
                   plottype=plottype,
                   colorbar=False,
                   cscale=cscale,
                   contours=contours,
                   ncontours=ncontours,
                   zlim=zlim,
                   elev=elev,
                   azim=azim,
                   theta=theta,
                   ih = ihh[ncomp],
                   sign=sign[ncomp],
                   xlab=xlab,
                   ylab=compname[ncomp],
                   clab=clab,
                   equal=equal,
                   cmap=cmap,
                   alpha=alpha,
                   xtorlabel=xtorlabel)

    xlim = plt.xlim()
    ylim = plt.ylim()
    xtext = xlim[0] + 0.*np.diff(xlim)
    ytext = ylim[0] + 0.95*np.diff(ylim)
    plt.text(xtext,ytext,xtorlabel)

def Plot(*args,**kwargs):
    """
    same as plot but sets different subplot parameters
    Can be used to set any kind of default options
    """
    
    fig = plt.plot(*args,**kwargs)
    plt.subplots_adjust(bottom=0.12, left=.14, right=0.94, top=.93)
    plt.grid(True)
    SetPowerLimits()

    return fig
    
def Semilogx(*args,**kwargs):
    """
    same as plot but sets different subplot parameters
    Can be used to set any kind of default options
    """
    
    fig = plt.semilogx(*args,**kwargs)
    plt.subplots_adjust(bottom=0.12, left=.13, right=0.93, top=.93)
    plt.grid(True)
    SetPowerLimits_y()

    return fig
    
def Semilogy(*args,**kwargs):
    """
    same as plot but sets different subplot parameters
    Can be used to set any kind of default options
    """
    
    fig = plt.semilogy(*args,**kwargs)
    plt.subplots_adjust(bottom=0.12, left=.175, right=0.93, top=.93)
    plt.grid(True)
    SetPowerLimits_x()

    return fig
    
def Loglog(*args,**kwargs):
    """
    same as plot but sets different subplot parameters
    Can be used to set any kind of default options
    """
    
    fig = plt.loglog(*args,**kwargs)
    plt.subplots_adjust(bottom=0.12, left=.13, right=0.93, top=.93)
    plt.grid(True)

    return fig
    
def SetPowerLimits(xmini=-2,xmaxi=3,ymini=-2,ymaxi=4):
    """
    Sets the formatting of x and y labels in the current figure
    """
    
    ax=plt.gca()
    ax.xaxis.get_major_formatter().set_powerlimits((xmini,xmaxi));
    ax.yaxis.get_major_formatter().set_powerlimits((ymini,ymaxi));
    plt.draw();

def SetPowerLimits_x(xmini=-2,xmaxi=4):
    """
    Sets the formatting of x and y labels in the current figure
    """
    
    ax=plt.gca()
    ax.xaxis.get_major_formatter().set_powerlimits((xmini,xmaxi));
    plt.draw();

def SetPowerLimits_y(ymini=-2,ymaxi=4):
    """
    Sets the formatting of x and y labels in the current figure
    """
    
    ax=plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((ymini,ymaxi));
    plt.draw();
    
def ReadGlobals(filename):
    
    with open(filename,'r') as f:
        lines = f.readlines()

    Y = []
    D = []
    T = []
    for line in lines:
        y,d,t = [float(line.split()[i]) for i in range(3)]
        Y.append(y)
        D.append(d)
        T.append(t)

    Y = np.array(Y) # Years
    D = np.array(D) # Days
    T = np.array(T) # Temperatures
    t = (Y-1)*360+D
    

    return t,Y,D,T

def ReadNames_AV(filename, AV=True):
    """
    For Atmos and Vecode
    """
    ds = nc.Dataset(filename)
    
    names = []
    long_names = []
    standard_names = []

    if AV:
        for v in ds.variables:
            if 'long_name' in dir(ds[v]):
                names.append(str(ds.variables[v].name))
                long_names.append(str(ds.variables[v].long_name))
                standard_names.append(str(ds.variables[v].standard_name))
    else:
        DwnsN = ['Surface Temperature', 'Total Precipitation', 'Relative humidity']
        for n, v in enumerate(ds.variables):
            names.append(str(ds.variables[v].name))
            if n>2:
                long_names.append(DwnsN[n-3])
            else:
                long_names.append(str(ds.variables[v].name))

    #for i in range(len(names)):
    #    print(names[i], long_names[i])

    return names,long_names,standard_names

def ReadNames_O(filename):
    """
    For Ocean (no standard name)
    """
    ds = nc.Dataset(filename)
    
    names = []
    long_names = []
        
    for v in ds.variables:
        if 'long_name' in dir(ds[v]):
            names.append(str(ds.variables[v].name))
            long_names.append(str(ds.variables[v].long_name))

    return names,long_names
