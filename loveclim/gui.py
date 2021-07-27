import tkinter as tk
from tkinter import font
import numpy as np
import loveclim as lc
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import gaussian_filter
import os, sys
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import platform
if platform.system()=='Darwin':
    matplotlib.use('MacOSX')

exitcolor = "#df0101"
green = "#9afe2e"
blue = "#33adff"

# set the colormap and centre the colorbar
class MidpointNormalize(matplotlib.colors.Normalize):
    """
    Normalise the colorbar.
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


class loveclim_GUI:
    """
    Main window
    """

    def __init__(self, master):

        self.master = master
        self.buttonfont = font.Font(family='Helvetica',
                                    size=22,
                                    weight='bold')
        self.fieldfont = font.Font(family='Helvetica',
                                   size=16,
                                   weight='bold')
        self.textfont = font.Font(family='Helvetica',
                                  size=16)
        font.families()
        self.frame = tk.Frame(self.master)
        self.buttonfont = font.Font(family='Helvetica', size=22, weight='bold')
        font.families()

        path = os.getcwd()
        self.path = path
        # files for atmos
        self.atmospath = os.path.join(path, 'atmos')
        AllFiles = os.listdir(self.atmospath)
        lcfiles = [File for File in AllFiles if File.endswith('.nc')]
        lcfiles.sort()
        nlc = len(lcfiles)
        os.chdir(path)

        # files for downscaling
        self.oceanpath = os.path.join(path, 'ocean')
        AllFiles3 = os.listdir(self.oceanpath)
        lcfiles3 = [File for File in AllFiles3 if File.endswith('.nc')]
        lcfiles3.sort()
        nlc3 = len(lcfiles3)
        os.chdir(path)

        # files for downscaling
        self.downspath = os.path.join(path, 'downscaling')
        AllFiles2 = os.listdir(self.downspath)
        lcfiles2 = [File for File in AllFiles2 if File.endswith('.nc')]
        lcfiles2.sort()
        nlc2 = len(lcfiles2)
        os.chdir(path)
 
        if nlc!=0:
            # Loveclim Atmos (A) of Vecode (V) files
            self.AV_Open_b = tk.Button(self.frame,text="Open selected\n Atmos or Vecode File",
                                       bg=green,
                                       highlightbackground=green,
                                       font=self.buttonfont,
                                       command=self.AV_Open)
           # Scrollbar for the .nc files
            self.scrollbar = tk.Scrollbar(self.frame)

            # File listbox
            maxlen=0
            for filename in lcfiles:
                l = len(filename)
                if l>maxlen:
                    maxlen = l
            self.ncFiles_lb = tk.Listbox(self.frame,
                                         selectmode=tk.SINGLE,
                                         yscrollcommand=self.scrollbar.set,
                                         font=self.fieldfont,
                                         exportselection=False,
                                         width=maxlen)
            self.scrollbar.config(command=self.ncFiles_lb.yview)
            for filename in lcfiles:
                self.ncFiles_lb.insert(tk.END,filename)
            self.ncFiles_lb.selection_set(0)

        if nlc3!=0:
            # Loveclim Ocean (O) files (CLIO)
            self.O_Open_b = tk.Button(self.frame,text="Open selected\n Ocean File",
                                      bg=green,
                                      highlightbackground=green,
                                      font=self.buttonfont,
                                      command=self.O_Open)

            # File listbox ocean
            maxlen=0
            for filename in lcfiles3:
                l = len(filename)
                if l>maxlen:
                    maxlen = l
            self.ncFiles_lb3 = tk.Listbox(self.frame,
                                         selectmode=tk.SINGLE,
                                         yscrollcommand=self.scrollbar.set,
                                         font=self.fieldfont,
                                         exportselection=False,
                                         width=maxlen)
            self.scrollbar.config(command=self.ncFiles_lb3.yview)
            for filename in lcfiles3:
                self.ncFiles_lb3.insert(tk.END,filename)
            self.ncFiles_lb3.selection_set(0)

        if nlc2!=0:
            # Loveclim Downscaling (D) files (Downscaling)
            self.D_Open_b = tk.Button(self.frame,text="Open selected\n Downscaling File",
                                      bg=green,
                                      highlightbackground=green,
                                      font=self.buttonfont,
                                      command=self.D_Open)

            # File listbox downscaling
            maxlen=0
            for filename in lcfiles2:
                l = len(filename)
                if l>maxlen:
                    maxlen = l
            self.ncFiles_lb2 = tk.Listbox(self.frame,
                                         selectmode=tk.SINGLE,
                                         yscrollcommand=self.scrollbar.set,
                                         font=self.fieldfont,
                                         exportselection=False,
                                         width=maxlen)
            self.scrollbar.config(command=self.ncFiles_lb2.yview)
            for filename in lcfiles2:
                self.ncFiles_lb2.insert(tk.END,filename)
            self.ncFiles_lb2.selection_set(0)

        self.Quit_b = tk.Button(self.frame,
                                text="Exit",
                                bg=exitcolor,
                                highlightbackground=exitcolor,
                                command=self.frame.quit,
                                font=self.buttonfont)
        # Layout
        self.frame.grid()
        if nlc!=0:
            self.ncFiles_lb.grid(row=0,rowspan=1,column=0,sticky="ew")
        if nlc3!=0:
            self.ncFiles_lb3.grid(row=1,rowspan=1,column=0,sticky="ew")
        if nlc2!=0:
            self.ncFiles_lb2.grid(row=2,rowspan=1,column=0,sticky="ew")
        self.scrollbar.grid(row=0,rowspan=1,column=1,sticky="wsn")
        self.AV_Open_b.grid(row=0,column=2)
        self.O_Open_b.grid(row=1,column=2)
        if nlc2!=0:
            self.D_Open_b.grid(row=2,column=2)
        self.Quit_b.grid(row=3,column=2)
            

    def AV_Open(self):
        """
        Open Atmos or Vecode file
        """
        os.chdir(self.atmospath)
        id = self.ncFiles_lb.curselection()
        self.filename = self.ncFiles_lb.get(id)
        self.newWindow = tk.Toplevel(self.master)
        self.app = AV_netCDF_GUI(self.newWindow, self.filename)
        os.chdir(self.path)

    def D_Open(self):
        """
        Open Downscaling and Atmos file
        """
        os.chdir(self.atmospath)
        id2 = self.ncFiles_lb.curselection()
        self.filename_atmos = self.ncFiles_lb.get(id2)
        os.chdir(self.downspath)
        id = self.ncFiles_lb2.curselection()
        self.filename = self.ncFiles_lb2.get(id)
        self.newWindow = tk.Toplevel(self.master)
        self.app = AV_netCDF_GUI(self.newWindow, self.filename,
                                 AV=False,
                                 filename_atmos=self.filename_atmos,
                                 path=self.downspath,
                                 atmospath=self.atmospath,
                                 )
        os.chdir(self.path)

    def O_Open(self):
        """
        Open Ocean file
        """
        os.chdir(self.oceanpath)
        id = self.ncFiles_lb3.curselection()
        self.filename = self.ncFiles_lb3.get(id)
        self.newWindow = tk.Toplevel(self.master)
        self.app = O_netCDF_GUI(self.newWindow,self.filename)
        os.chdir(self.path)

class AV_netCDF_GUI(MidpointNormalize):
    """
    GUI for Atmos or Vecode files
    """
    
    def __init__(self,master,filename,AV=True,filename_atmos='',path='',atmospath=''):

        self.master = master
        self.buttonfont = font.Font(family='Helvetica',
                                    size=22,
                                    weight='bold')
        self.fieldfont = font.Font(family='Helvetica',
                                   size=16,
                                   weight='bold')
        self.textfont = font.Font(family='Helvetica',
                                  size=16)
        font.families()

        self.frame = None

        self.plotexists = False
        self.fieldname = ""
        self.itime=0
        self.filename = filename
        # bool True whether it is Atmos or vecode, False for downsaling
        self.AV = AV
        self.filename_atmos = filename_atmos
        self.path = path
        self.atmospath = atmospath
        self.oceanpath = os.path.join(self.path, 'ocean.dat')
        self.anomaly = False
        self.label = 0
        # old fields values
        self.oldfieldname = ''
        self.oldtypeval = ''
        self.olditime = ''

        # Read data
        self.ds = nc.Dataset(self.filename)
        if not self.AV:
            os.chdir(self.atmospath)
            self.ds_atmos = nc.Dataset(self.filename_atmos)
            os.chdir(self.path)
        # Prepare projections
        self.projections = [ccrs.PlateCarree,
                            ccrs.Orthographic,
                            ccrs.Mercator,
                            ccrs.Robinson]
        # For Mora-like maps
        Threshold_T = np.array([26.5, 27. , 27.5, 28. , 28.5, 29. , 29.5, 30. , 30.5, 31. , 31.5,
                                32. , 32.5, 33. , 33.5, 34. , 34.5, 35. , 35.5, 36. , 36.5, 37. ,
                                37.5, 38. , 38.5, 39. , 39.5, 40. , 40.5, 41. , 41.5, 42. , 42.5,
                                43. , 43.5, 44. , 44.5, 45. , 45.5, 46. , 46.5, 47. , 47.5, 48. ,
                                48.5])
        Threshold_H = np.array([100.109,  89.015,  82.441,  76.078,  71.895,  67.712,  63.529,
                                59.493,  57.152,  54.812,  52.471,  50.131,  47.79 ,  45.449,
                                43.109,  40.807,  39.081,  37.355,  35.628,  33.902,  32.175,
                                30.449,  28.723,  26.996,  25.291,  24.016,  22.741,  21.466,
                                20.191,  18.916,  17.641,  16.366,  15.091,  13.817,  12.55 ,
                                11.614,  10.678,   9.742,   8.805,   7.869,   6.933,   5.997,
                                5.06 ,   4.124,   3.188])
        # the following function gives the temperature deadly threshold 
        # as a function of the local relative humidity
        self.MoraDeadlyThreshold = interp1d(Threshold_H,Threshold_T)
        # Setup GUI
        self._SetupGui()

    def _SetupGui(self):
        """
        Setup the Gui
        """
        if self.frame is not None:
            self.frame.destroy()
        self.frame = tk.Frame(self.master)

        self.Quit_b = tk.Button(self.frame,
                                text="Exit",
                                bg=exitcolor,
                                highlightbackground=exitcolor,
                                command=self.frame.quit,
                                font=self.buttonfont)
         
        # List of names to exclude (not fields)
        if self.AV:
            excluded_names = ["lon","lat","time","P_T3","P_T4","P_U3"]
        else:
            excluded_names = ["x","y","time"]

        # get variable names (only if they have a 'long_name')
        names,long_names,standard_names = lc.ReadNames_AV(self.filename, AV=self.AV)
        self.names = names
        self.long_names = long_names
        self.standard_names = standard_names
        units = []
        for i in range(len(names)):
            try:
                units.append(self.ds[names[i]].units)
            except AttributeError:
                if names[i]=='Tann':
                    un = '°C'
                elif names[i]=='Acc':
                    un = 'cm.yr-1'
                elif names[i]=='relhum':
                    un = '1'
                units.append(un)

        self.lb_names = ['{:} ({:})'.format(long_names[i],units[i])
                for i in range(len(names)) if names[i] not in excluded_names]

        maxlen=0
        for name in self.lb_names:
            l = len(name)
            if l>maxlen:
                maxlen = l

        # Scrollbar for this widget
        self.scrollbar = tk.Scrollbar(self.frame)
        self.Field_lb = tk.Listbox(self.frame,
                                   selectmode=tk.SINGLE,
                                   yscrollcommand=self.scrollbar.set,
                                   font=self.fieldfont,
                                   exportselection=False,
                                   width=maxlen)
        self.scrollbar.config(command=self.Field_lb.yview)
        for name in self.lb_names:
            self.Field_lb.insert(tk.END,name)
        if ('Relative Humidity' in long_names) and ('Temperature at 2 Meter' in long_names):
            self.Field_lb.insert(tk.END,'Mora Deadly Threshold (K)')
            self.lb_names.append('Temperature above threshold (K)')
        self.Field_lb.selection_set(0)
        self.Field_lb.event_generate("<<ListboxSelect>>")
        self.Field_lb.bind("<<ListboxSelect>>",self._Enter)
            
        # itime spinbox
        self.Itime_l = tk.Label(self.frame,
                                text="itime:",
                                font=self.textfont)
        self.Itime_sb = tk.Spinbox(self.frame,
                                   from_=1,
                                   to=len(self.ds['time'][:]),
                                   font=self.textfont)
        self.Itime_sb.bind('<Return>',self._EnterIZ)
        self.Itime_sb.bind('<KP_Enter>',self._EnterIZ)
        self.Time_var = tk.StringVar(self.frame)
        self.Time_var.set('')
        self.Time_l = tk.Label(self.frame,
                                textvariable=self.Time_var,
                                font=self.textfont)

        # zvar spinbox (for fields such as Air Temperature, Wind direction etc.)
        self.Zvar_l = tk.Label(self.frame,
                               text="zvar:",
                               font=self.textfont)
        self.Zvar_sb = tk.Spinbox(self.frame,
                                   from_=1,
                                   to=4,
                                   font=self.textfont)
        self.Zvar_sb.bind('<Return>',self._EnterIZ)
        self.Zvar_sb.bind('<KP_Enter>',self._EnterIZ)
        self.Zvar_var = tk.StringVar(self.frame)
        self.Zvar_var.set('')
        self.ZvarActual_l = tk.Label(self.frame,
                                     textvariable=self.Zvar_var,
                                     font=self.textfont)
        # Projection option menu
        self.Proj_l = tk.Label(self.frame,text = "Projection:",
                               font=self.textfont)
        self.Proj_options = ["PlateCarree", "Orthographic","Mercator","Robinson"]
        self.Proj_var = tk.StringVar(self.frame)
        self.Proj_var.set(self.Proj_options[0])
        self.Proj_om = tk.OptionMenu(self.frame,
                                      self.Proj_var,*self.Proj_options,
                                      command=self._Replot)
        self.Proj_om.config(font=self.fieldfont)
        menu = self.frame.nametowidget(self.Proj_om.menuname)
        menu.config(font=self.fieldfont)

        # Central Longitude scale
        res = 10
        self.Clon_s = tk.Scale(self.frame,
                               orient='horizontal',
                               command=self._Replot,
                               from_=-180,
                               to=180,
                               resolution=res,
                               label='Central longitude',
                               tickinterval=45,
                               font=self.textfont)
        self.Clon_s.set(0)

        # Central Latitude scale (for Orthographic)
        self.Clat_s = tk.Scale(self.frame,
                               orient='horizontal',
                               command=self._Replot,
                               from_=-90,
                               to=90,
                               resolution=res,
                               label='Central latitude',
                               tickinterval=45,
                               font=self.textfont)
        self.Clat_s.set(0)

        # Colorbar
        self.CCbar_l0 = tk.Label(self.frame,
                                text="Colorbar",
                                font=self.textfont)
        self.CCbar_l1 = tk.Label(self.frame,
                                text="min:",
                                font=self.textfont)
        self.CCbar_l2 = tk.Label(self.frame,
                                text="max:",
                                font=self.textfont)
        self.CCbar_sb = tk.Spinbox(self.frame,
                                   command=self._Replot,
                                   values=('auto\ for\ first\ plot'),
                                   font=self.textfont)
        self.CCbar_sb.bind('<Return>',self._Enter)
        self.CCbar_sb.bind('<KP_Enter>',self._Enter)
        self.CCbar_ssb = tk.Spinbox(self.frame,
                                   command=self._Replot,
                                   values=('auto\ for\ first\ plot'),
                                   font=self.textfont)
        self.CCbar_ssb.bind('<Return>',self._Enter)
        self.CCbar_ssb.bind('<KP_Enter>',self._Enter)

        # absolute or relatives values
        self.RItime_l = tk.Label(self.frame,
                                text="itime preind:",
                                font=self.textfont)
        self.RItime_sb = tk.Spinbox(self.frame,
                                   from_=1,
                                   to=len(self.ds['time'][:]),
                                   font=self.textfont)
        self.RItime_sb.bind('<Return>',self._EnterIZ)
        self.RItime_sb.bind('<KP_Enter>',self._EnterIZ)
        self.Time_var = tk.StringVar(self.frame)
        self.Time_var.set('')
        self.Time_l = tk.Label(self.frame,
                                textvariable=self.Time_var,
                                font=self.textfont)
        self.Typeval_l = tk.Label(self.frame,
                                  text="Type values:",
                                  font=self.textfont)
        self.Typeval = tk.Spinbox(self.frame,
                                   values=('absolute', 'anomaly'),
                                   font=self.textfont)
        self.Typeval.bind('<Return>',self._Enter)
        self.Typeval.bind('<KP_Enter>',self._Enter)

        # nb contour
        self.cont_l = tk.Label(self.frame,
                                text="Nb contour:",
                                font=self.textfont)
        self.cont_sb = tk.Spinbox(self.frame,
                                   command=self._Replot,
                                   from_=10,
                                   to=100,
                                   font=self.textfont)
        self.cont_sb.bind('<Return>',self._Enter)
        self.cont_sb.bind('<KP_Enter>',self._Enter)
 
        # Plotbutton
        self.Plot_b = tk.Button(self.frame,text="Plot",
                                bg=green,
                                highlightbackground=green,
                                font=self.buttonfont,
                                command=self._Plot)
        # Layout
        self.frame.grid()
        row=0
        self.Field_lb.grid(row=row,rowspan=13,column=0,sticky="esn")
        self.scrollbar.grid(row=row,rowspan=13,column=1,sticky="wsn")
        row=0
        self.Itime_l.grid(row=row,column=2,sticky="e")
        self.Itime_sb.grid(row=row,column=3,sticky="w")
        self.Time_l.grid(row=row,column=4)
        row = row+1
        self.Zvar_l.grid(row=row,column=2,sticky="e")
        self.Zvar_sb.grid(row=row,column=3,sticky="w")
        self.ZvarActual_l.grid(row=row,column=4)
        row = row+1
        self.Clon_s.grid(row=row,column=2,columnspan=3,sticky="ew")
        row = row+1
        self.Clat_s.grid(row=row,column=2,columnspan=3,sticky="ew")
        row=row+1
        self.CCbar_l0.grid(row=row,column=2,sticky="e")
        row = row+1
        self.CCbar_l1.grid(row=row,column=2,sticky="e")
        self.CCbar_sb.grid(row=row,column=3,sticky="w")
        self.CCrow = row
        row = row+1
        self.CCbar_l2.grid(row=row,column=2,sticky="e")
        self.CCbar_ssb.grid(row=row,column=3,sticky="w")
        row = row+1
        self.Typeval_l.grid(row=row,column=2,sticky="e")
        self.Typeval.grid(row=row,column=3,sticky="w")
        row += 1
        self.RItime_l.grid(row=row,column=2,sticky="e")
        self.RItime_sb.grid(row=row,column=3,sticky="w")
        row=row+1
        self.cont_l.grid(row=row,column=2,sticky="e")
        self.cont_sb.grid(row=row,column=3,sticky="w")
        row = row+1
        self.Proj_l.grid(row=row,column=2,sticky="e")
        self.Proj_om.grid(row=row,column=3,sticky="w")
        row=row+1
        self.Plot_b.grid(row=row,column=2,columnspan=3,sticky="ew")
        row=row+1
        self.Quit_b.grid(row=row,column=2,columnspan=3,sticky="ew")

    def load_CC(self, datmin, datmax):
        """
        loqd colorbar values
        """
        # update colorbar
        self.CCbar_sb = tk.Spinbox(self.frame,
                                   command=self._Replot,
                                   from_=datmin,
                                   to=datmax,
                                   increment=0.01,
                                   font=self.textfont)
        self.CCbar_sb.bind('<Return>',self._Enter)
        self.CCbar_sb.bind('<KP_Enter>',self._Enter)
        var = tk.StringVar(self.frame)
        var.set("{:d}".format(int(datmax+1)))
        self.CCbar_ssb = tk.Spinbox(self.frame,
                                   command=self._Replot,
                                   from_=datmin,
                                   to=datmax,
                                   increment=0.01,
                                   font=self.textfont,
                                   textvariable=var)
        self.CCbar_ssb.bind('<Return>',self._Enter)
        self.CCbar_ssb.bind('<KP_Enter>',self._Enter)
        row = self.CCrow
        self.CCbar_sb.grid(row=row,column=3,sticky="w")
        row = row+1
        self.CCbar_ssb.grid(row=row,column=3,sticky="w")

    def _Preplot(self):
        londws, latdws, datadws = [], [], []
        if not self.plotexists:
            if self.AV:
                lon,lat,data = self._GetPlotFields()
                lon,lat,data = self._PeriodicBoundaryConditions(lon,lat,data)
            else:
                # loading data
                londws,latdws,datadws,lon,lat,data = self._GetPlotFields(self.typeval)
                # interp data
                lon,lat,data = self._Perio_Interp(londws,latdws,datadws,lon,lat,data)

            # save data
            self.DATA = [lon, lat, data]
        else:
            lon = self.DATA[0]
            lat = self.DATA[1]
            data = self.DATA[2]

        return lon,lat,data

    def _Plot(self):
        """
        Plot using options
        """
        # load data
        self._ReadParams()
        lon,lat,data = self._Preplot()

        # projection
        projectionfun = self.projections[self.Proj_options.index(self.proj)]
        if self.proj=='Orthographic':
            projection = projectionfun(self.clon,self.clat)
        elif self.proj=='Mercator':
            projection = projectionfun(central_longitude=self.clon, min_latitude=-80,
                                       max_latitude=84)
        else:
            projection = projectionfun(self.clon)
        self.ax = plt.axes(projection=projection, label=self.label)
        self.ax.set_global()
        self.ax.coastlines()
        bodr = cartopy.feature.NaturalEarthFeature(category='cultural',
            name='admin_0_boundary_lines_land', scale='50m', facecolor='none', alpha=0.7)
        self.ax.add_feature(bodr, linestyle='--', edgecolor='k', alpha=1)

        # colorbar
        if self.plotexists:
            try:
                vmin = float(self.CCbar_sb.get())
                vmax = float(self.CCbar_ssb.get())
            except ValueError:
                    vmin, vmax = (np.amin(data)), (np.amax(data))
                    self.load_CC(vmin, vmax)
        else:
            vmin, vmax = (np.amin(data)), (np.amax(data))
            self.load_CC(vmin, vmax)
        if (vmin>=vmax):
            vmin, vmax = (np.amin(data)), (np.amax(data))
        cont = max(10,int(self.cont_sb.get()))
        if not self.fieldname=='mora':
            cmap = plt.get_cmap('coolwarm')
        else:
            cmap = plt.get_cmap('Greys')
        levels = np.linspace(start=vmin, stop=vmax, num=cont)

        # plot
        # load data
        tmp_atmos = data.copy()
        tmp_atmos[tmp_atmos<vmin] = vmin
        tmp_atmos[tmp_atmos>vmax] = vmax
 
        # contourf
        im = self.ax.contourf(lon-self.clon, lat-self.clat, tmp_atmos, levels=levels,
                              transform=ccrs.PlateCarree(self.clon),
                              cmap=cmap)
 
        # post process image
        self.cb = plt.colorbar(im)
        label = self.lb_names[self.Field_lb.curselection()[0]]
        if self.anomaly:
            label += ' (anomaly)'
        self.cb.ax.set_ylabel(label)
        self.Time_var.set('Time label = {:d}'.format(int(self.ds['time'][self.itime])))
        try:
            self.Zvar_var.set('{:} @ {:} = {:d} {:}'.format(
                self.Field_lb.get(self.Field_lb.curselection()),
                self.zvar_longname,
                int(self.zvar_val),
                self.zvar_units))
        except AttributeError:
            self.Zvar_var.set('')
        
        # end of function
        self.label += 1
        self.plotexists = True

    def detect_ocean(self, londws, latdws):
        ## Get whether the points are on land using Basemap.is_land
        print('File ocean.dat not found\n making ocean.dat ---> please wait')
        print(' this may take a few minutes...')
        lon_grid, lat_grid = np.meshgrid(londws-180.,latdws)
        bm = Basemap(projection='cyl', llcrnrlat=-89.75, urcrnrlat=89.75,
                     llcrnrlon=-180.0, urcrnrlon=180.0, resolution='c', area_thresh=1)
        f = np.vectorize(bm.is_land)
        xpt, ypt = bm(lon_grid, lat_grid)
        xpts = xpt.shape[1]
        xpt = np.hstack((xpt[:,xpts//2:], xpt[:,:xpts//2]))
        basemap_land_mask = f(xpt,ypt)
        np.savetxt(fname=self.oceanpath, X=np.logical_not(basemap_land_mask))
        print('File ocean.dat made')

    def _Replot(self,dummy=0):
        """
        If there is no active plot, do nothing
        Otherwise, recall Plot
        dummy is a the scale argument (nslide or iphi)
        """
        self._Plot()

    def _ReadParams(self):

        # variables that need a new data load
        self.fieldname = self._Name(self.Field_lb.get(self.Field_lb.curselection()))
        self.itime = int(self.Itime_sb.get())-1
        self.anomaly = self.Typeval.get()=='anomaly'
        self.typeval = self.anomaly
        self.plotexists = (self.oldfieldname==self.fieldname) and (self.olditime==self.itime) \
                and (self.oldtypeval==self.typeval)
        if not self.plotexists:
            self.oldfieldname = self.fieldname
            self.olditime = self.itime
            self.oldtypeval = self.typeval

        # variables with same data
        self.zvar = int(self.Zvar_sb.get())-1
        self.proj = self.Proj_var.get()
        self.clon = float(self.Clon_s.get())
        self.clat = float(self.Clat_s.get())
        self.Ritime = int(self.RItime_sb.get())-1

    def _Enter(self,even):
        self._Replot()

    def _EnterIZ(self,even):
        self.plotexists = False
        self._Replot()

    def _GetPlotFields(self, anomaly=False):
        """
        Read fields from nc file
        """

        if self.AV:
            lon = self.ds['lon'][:]
            lat = self.ds['lat'][:]
        
        else:
            lon = np.linspace(start=-179.75, stop=179.75, num=720)
            lat = np.linspace(start=-89.75, stop=89.75, num=360)

        if self.fieldname=='mora':
            # Temperature at 2 meters in degree C
            t2m = self.ds['t2m'][:][self.itime,:,:]-273.15
            # relative humidity (in percentage)
            rh = self.ds['r'][:][self.itime,:,:]*100
            # take the maximum of the data and 0 (actually -.1 to always have something to plot)
            # in order to emphasize the above threshold part of the data
            data = np.maximum(t2m.data - self.MoraDeadlyThreshold(rh.data),-1.)
            data = t2m.data - self.MoraDeadlyThreshold(rh.data)

            rawdata = np.ma.masked_array(data,mask=t2m.mask)
            return lon,lat,rawdata


        rawdata = self.ds[self.fieldname][:]
        fieldvar = self.ds[self.fieldname]

        if len(rawdata.shape)==3:
            try:
                del self.zvar_longname
                del self.zvar_val
                del self.zvar_units
            except AttributeError:
                pass

            # reconstruire les donnees avec atmphys0.f
            if self.AV:
                return lon,lat,rawdata[self.itime,:,:]

            else:
                # loat atmos data
                lon_atmos = self.ds_atmos['lon'][:]
                lat_atmos = self.ds_atmos['lat'][:]
                ds_atmos = self.ds_atmos
                if self.fieldname=='Tann':
                    fieldname_atmos = 'ts'
                elif self.fieldname=='Acc':
                    fieldname_atmos = 'pp'
                elif self.fieldname=='relhum':
                    fieldname_atmos = 'r'
                rawdata_atmos = self.ds_atmos[fieldname_atmos][:]
                if self.fieldname=='Tann':
                    rawdata_atmos -= 273.15

                # select times on 31 years
                maxtime = rawdata.shape[0]
                amaxtime = rawdata_atmos.shape[0]
                itime = self.itime
                Ritime = self.Ritime
                lb, up, mm = self.ibounds(itime, maxtime)
                lba, upa, mma = self.ibounds(Ritime, amaxtime)

                # if anomaly
                if anomaly:
                    rawdata = rawdata[itime-lb:itime+up,:,:] - \
                              rawdata[Ritime-lba:Ritime+upa,:,:]
                    rawdata_atmos = \
                              rawdata_atmos[itime-lb:itime+up,:,:] - \
                              rawdata_atmos[Ritime-lba:Ritime+upa,:,:]
                else:
                    rawdata = rawdata[itime-lb:itime+up,:,:]
                    rawdata_atmos = rawdata_atmos[itime-lb:itime+up,:,:]

                # mean data over 31 years
                rawdata = np.mean(a=rawdata, axis=0)
                rawdata_atmos = np.mean(a=rawdata_atmos, axis=0)

                return lon,lat,rawdata,lon_atmos,lat_atmos,rawdata_atmos

        elif len(rawdata.shape)==4:
            maxzvar = rawdata.shape[1]
            self.zvar_actual = min(self.zvar+1,maxzvar)-1
            self.zvar_longname = self.ds[fieldvar.dimensions[1]].long_name
            self.zvar_val = self.ds[fieldvar.dimensions[1]][self.zvar_actual]
            self.zvar_units = self.ds[fieldvar.dimensions[1]].units

            return lon,lat,rawdata[self.itime,self.zvar_actual,:,:]

    def ibounds(self, itime, maxtime):
        if (itime-15)>=0:
            if (itime+15)<=maxtime: 
                lb, up = 15, 15
            else:
                lb, up = 15, (maxtime-itime)
        else:
            if (itime+15)<=maxtime: 
                lb, up = itime, 15
            else:
                lb, up = itime, (maxtime-itime)
        mm = (up+lb+1)

        return lb, up, mm

    def _PeriodicBoundaryConditions(self,lon,lat,data,atmos=True):
        """
        Close the data set by setting additional values at the end of the arrays
        """
        clon = self.clon
        clat = self.clat
        if atmos:
            last = 360.
        else:
            last = 180.

        # Two possibilities: the mask is an array or it is just equal to False
        if type(data.mask)==np.bool_:
            bc_lon = np.zeros(len(lon)+1)
            bc_lat = np.zeros(lat.size)
            bc_data = np.zeros((data.shape[0],len(lon)+1))
            bc_lon[:-1] = lon[:]-clon
            bc_lon[-1] = last-clon
            bc_lat = lat-clat
            bc_data[:,:-1] = data.data[:,:]
            bc_data[:,-1] = data.data[:,0]

            return bc_lon,bc_lat,np.ma.array(bc_data,mask=False)        
        else:
            bc_lon = np.zeros(len(lon)+1)
            bc_data = np.zeros((data.shape[0],len(lon)+1))
            bc_lat = np.zeros(lat.size)
            bc_mask = np.zeros((data.shape[0],len(lon)+1))
            bc_lon[:-1] = lon[:]-clon
            bc_lon[-1] = last-clon
            bc_lat = lat-clat
            bc_data[:,:-1] = data.data[:,:]
            bc_mask[:,:-1] = data.mask[:,:]
            bc_data[:,-1] = data.data[:,0]
            bc_mask[:,-1] = data.mask[:,0]

            return bc_lon,bc_lat,np.ma.array(bc_data,mask=bc_mask.astype(bool))

    def _Perio_Interp(self,londws,latdws,datadws,lon,lat,data):

        # periodic dwns
        lons = londws.size
        lon_forint = np.hstack((londws[lons//2:]-360., londws, londws[:lons//2]+360.))
        data_forint = np.hstack((datadws[:,lons//2:], datadws, datadws[:,:lons//2]))
        int_data = interp2d(x=lon_forint, y=latdws, z=data_forint)
        londws = np.linspace(0, 360, londws.size)
        datadws = np.ma.asarray(int_data(londws, latdws))

        # oceans
        if not os.path.exists(self.oceanpath):
            self.detect_ocean(londws, latdws) 

        # periodic original
        lons = lon.size
        # make a line above and under the map for poles
        lat_forint = np.hstack(([-90], lat, [90]))
        lon_forint = np.hstack((lon[lons//2:]-360., lon, lon[:lons//2]+360.))
        data_forint = np.vstack((data[0,:], data, data[-1,:]))
        data_forint = np.hstack((data_forint[:,lons//2:], data_forint, data_forint[:,:lons//2]))
        int_data = interp2d(x=lon_forint, y=lat_forint, z=data_forint)
        lon = np.linspace(0, 360, londws.size)
        lat = np.linspace(-90, 90, latdws.size)
        int_data = np.ma.asarray(int_data(lon, lat))

        # treat land and ocean
        ocean = np.asarray(np.loadtxt(self.oceanpath), dtype=bool)
        outland = datadws<np.amin(int_data)
        ocean = np.logical_or(ocean, outland)
        datadws[ocean] = int_data[ocean]

        return londws,latdws,datadws


    def _Name(self,long_name):
        """
        Maps long names to names
        """
        
        if long_name == 'Mora Deadly Threshold (K)':
            return 'mora'
        else:
            return self.names[self.long_names.index(long_name.split(' (')[0])]
























# ocean
class O_netCDF_GUI():
    """
    GUI for Ocean files
    """
    
    def __init__(self,master,filename):

        self.master = master
        self.buttonfont = font.Font(family='Helvetica',
                                    size=22,
                                    weight='bold')
        self.fieldfont = font.Font(family='Helvetica',
                                   size=16,
                                   weight='bold')
        self.textfont = font.Font(family='Helvetica',
                                  size=16)
        font.families()

        self.frame = None

        self.plotexists = False
        self.fieldname = ""
        self.itime=0
        self.filename = filename
        # Read data
        self.ds = nc.Dataset(self.filename)
        # Prepare projections
        self.projections = [ccrs.Orthographic,
                            ccrs.Mercator,
                            ccrs.Robinson,
                            ccrs.PlateCarree]
        # Setup GUI
        self._SetupGui()

    def _SetupGui(self):
        """
        Setup the Gui
        """
        if self.frame is not None:
            self.frame.destroy()
        self.frame = tk.Frame(self.master)

        self.Quit_b = tk.Button(self.frame,
                                text="Exit",
                                bg=exitcolor,
                                highlightbackground=exitcolor,
                                command=self.frame.quit,
                                font=self.buttonfont)

        # List of names to exclude (not fields)
        excluded_names = ["ptlon","ptlat",
                          "pulon","pulat",
                          "tlon","tlat",
                          "tlonp","tlatp",
                          "tlon_bounds","tlat_bounds",
                          "ulon","ulat",
                          "ulonp","ulatp",
                          "ulon_bounds","ulat_bounds",
                          "tdepth","wdepth","wedges",
                          "sflat","sfdepth","sfedges",
                          "basidx",
                          "time",
                          "angle","dxs1","dxs2","dxc1","dxc2",
                          "area","tmask","umask"]
        # get variable names (only if they have a 'long_name')
        names,long_names = lc.ReadNames_O(self.filename)
        self.names = names
        self.long_names = long_names
        self.lb_names = ['{:} ({:})'.format(long_names[i],self.ds[names[i]].units)
                for i in range(len(names)) if names[i] not in excluded_names]

        maxlen=0
        for name in self.lb_names:
            l = len(name)
            if l>maxlen:
                maxlen = l

        # Scrollbar for this widget
        self.scrollbar = tk.Scrollbar(self.frame)
        self.Field_lb = tk.Listbox(self.frame,
                                   selectmode=tk.SINGLE,
                                   yscrollcommand=self.scrollbar.set,
                                   font=self.fieldfont,
                                   exportselection=False,
                                   width=maxlen)
        self.scrollbar.config(command=self.Field_lb.yview)
        for name in self.lb_names:
            self.Field_lb.insert(tk.END,name)
        self.Field_lb.selection_set(0)
        self.Field_lb.event_generate("<<ListboxSelect>>")
        self.Field_lb.bind("<<ListboxSelect>>",self._Enter)
            
        # itime spinbox
        self.Itime_l = tk.Label(self.frame,
                                text="itime:",
                                font=self.textfont)
        self.Itime_sb = tk.Spinbox(self.frame,
                                   command=self._Replot,
                                   from_=1,
                                   to=len(self.ds['time'][:]),
                                   font=self.textfont)
        self.Itime_sb.bind('<Return>',self._Enter)
        self.Itime_sb.bind('<KP_Enter>',self._Enter)
        self.Time_var = tk.StringVar(self.frame)
        self.Time_var.set('')
        self.Time_l = tk.Label(self.frame,
                                textvariable=self.Time_var,
                                font=self.textfont)

        # Depth spinbox (for fields such as Air Temperature, Wind direction etc.)
        self.Depth_l = tk.Label(self.frame,
                               text="depth:",
                               font=self.textfont)
        self.Depth_sb = tk.Spinbox(self.frame,
                                   command=self._Replot,
                                   from_=1,
                                   to=max(len(self.ds['tdepth'][:]),
                                          len(self.ds['wdepth'][:])),
                                   font=self.textfont)
        self.Depth_sb.bind('<Return>',self._Enter)
        self.Depth_sb.bind('<KP_Enter>',self._Enter)
        self.Depth_var = tk.StringVar(self.frame)
        self.Depth_var.set('')
        self.DepthActual_l = tk.Label(self.frame,
                                     textvariable=self.Depth_var,
                                     font=self.textfont)
        # Projection option menu
        self.Proj_l = tk.Label(self.frame,text = "Projection:",
                               font=self.textfont)
        self.Proj_options = ["Orthographic","Mercator","Robinson","PlateCarree"]
        self.Proj_var = tk.StringVar(self.frame)
        self.Proj_var.set(self.Proj_options[0])
        self.Proj_om = tk.OptionMenu(self.frame,
                                      self.Proj_var,*self.Proj_options,
                                      command=self._Replot)
        self.Proj_om.config(font=self.fieldfont)
        menu = self.frame.nametowidget(self.Proj_om.menuname)
        menu.config(font=self.fieldfont)

        # Central Longitude scale
        self.Clon_s = tk.Scale(self.frame,
                               orient='horizontal',
                               command=self._Replot,
                               from_=-180,
                               to=180,
                               resolution=10,
                               label='Central longitude',
                               tickinterval=45,
                               font=self.textfont)
        self.Clon_s.set(0)

        # Central Latitude scale (for Orthographic)
        self.Clat_s = tk.Scale(self.frame,
                               orient='horizontal',
                               command=self._Replot,
                               from_=-90,
                               to=90,
                               resolution=10,
                               label='Central latitude',
                               tickinterval=45,
                               font=self.textfont)
        self.Clat_s.set(0)

        # Plotbutton
        self.Plot_b = tk.Button(self.frame,text="Plot",
                                bg=green,
                                highlightbackground=green,
                                font=self.buttonfont,
                                command=self._Plot)
        # Layout
        self.frame.grid()
        row=0
        self.Field_lb.grid(row=row,rowspan=13,column=0,sticky="esn")
        self.scrollbar.grid(row=row,rowspan=13,column=1,sticky="wsn")
        row=0
        self.Itime_l.grid(row=row,column=2,sticky="e")
        self.Itime_sb.grid(row=row,column=3,sticky="w")
        self.Time_l.grid(row=row,column=4)
        row = row+1
        self.Depth_l.grid(row=row,column=2,sticky="e")
        self.Depth_sb.grid(row=row,column=3,sticky="w")
        self.DepthActual_l.grid(row=row,column=4)
        row = row+1
        self.Clon_s.grid(row=row,column=2,columnspan=3,sticky="ew")
        row = row+1
        self.Clat_s.grid(row=row,column=2,columnspan=3,sticky="ew")
        row=row+1
        self.Proj_l.grid(row=row,column=2,sticky="e")
        self.Proj_om.grid(row=row,column=3,sticky="w")
        row=row+1
        self.Plot_b.grid(row=row,column=2,columnspan=3,sticky="ew")
        row=row+1
        self.Quit_b.grid(row=row,column=2,columnspan=3,sticky="ew")

    def _Plot(self):
        """
        Plot using options
        """

        self._ReadParams()

        lon,lat,data = self._GetPlotFields()
        projectionfun = self.projections[self.Proj_options.index(self.proj)]
        if self.proj=='Orthographic':
            projection = projectionfun(self.clon,self.clat)
        elif self.proj=='Mercator':
            projection = projectionfun(central_longitude=self.clon,
                                       min_latitude=-80,
                                       max_latitude=84)
        else:
            projection = projectionfun(self.clon)
        
        self.ax = plt.axes(projection=projection)
        self.ax.set_global()
        self.ax.coastlines()
        cmap = plt.get_cmap('coolwarm')
        im = self.ax.scatter(lon, lat, c=data,marker='s',transform=ccrs.PlateCarree(),cmap=cmap)
        self.cb = plt.colorbar(im)
        self.Time_var.set('Time label = {:d}'.format(int(self.ds['time'][self.itime])))
        try:
            self.Depth_var.set('{:} @ {:} = {:d} {:}'.format(self.Field_lb.get(self.Field_lb.curselection()),
                                                            self.depth_longname,
                                                            int(self.depth_val),
                                                            self.depth_units))
        except AttributeError:
            self.Depth_var.set('')
        
        self.plotexists = True

    def _Replot(self,dummy=0):
        """
        If there is no active plot, do nothing
        Otherwise, recall Plot
        dummy is a the scale argument (nslide or iphi)
        """

        # check for figure 1
        if self.plotexists:
            self.ax.clear()
            self.cb.remove()
            self._Plot()

    def _ReadParams(self):

        self.fieldname = self._Name(self.Field_lb.get(self.Field_lb.curselection()))
        self.itime = int(self.Itime_sb.get())-1
        self.depth = int(self.Depth_sb.get())-1
        self.proj = self.Proj_var.get()
        self.clon = float(self.Clon_s.get())
        self.clat = float(self.Clat_s.get())

    def _Enter(self,even):
        self._Replot()

    def _GetPlotFields(self):
        """
        Read fields from nc file
        """

        if self.fieldname=='fcor':
            lon = self.ds['ulon'][:]
            lat = self.ds['ulat'][:]
        else:
            lon = self.ds['tlon'][:]
            lat = self.ds['tlat'][:]
        rawdata = self.ds[self.fieldname][:]
        fieldvar = self.ds[self.fieldname]

        N = np.prod(lon.shape)

        if len(rawdata.shape)==2:
            try:
                del self.depth_longname
                del self.depth_val
                del self.depth_units
            except AttributeError:
                pass

            lonf = lon.data.flatten()
            latf = lat.data.flatten()
            datf = rawdata.flatten()

            lonfok = lonf[np.invert(datf.mask)]
            latfok = latf[np.invert(datf.mask)]
            datfok = datf.data[np.invert(datf.mask)]

            return lonfok,latfok,datfok
        elif len(rawdata.shape)==3:
            try:
                del self.depth_longname
                del self.depth_val
                del self.depth_units
            except AttributeError:
                pass

            lonf = lon.data.flatten()
            latf = lat.data.flatten()
            datf = rawdata[self.itime,:,:].flatten()

            lonfok = lonf[np.invert(datf.mask)]
            latfok = latf[np.invert(datf.mask)]
            datfok = datf.data[np.invert(datf.mask)]

            return lonfok,latfok,datfok
        elif len(rawdata.shape)==4:
            maxdepth = rawdata.shape[1]
            self.depth_actual = min(self.depth+1,maxdepth)-1
            self.depth_longname = self.ds[fieldvar.dimensions[1]].long_name
            self.depth_val = self.ds[fieldvar.dimensions[1]][-(self.depth_actual+1)]
            self.depth_units = self.ds[fieldvar.dimensions[1]].units

            lonf = lon.data.flatten()
            latf = lat.data.flatten()
            datf = rawdata[self.itime,-(self.depth_actual+1),:,:].flatten()

            lonfok = lonf[np.invert(datf.mask)]
            latfok = latf[np.invert(datf.mask)]
            datfok = datf.data[np.invert(datf.mask)]

            return lonfok,latfok,datfok

    def _Name(self,long_name):
        """
        Maps long names to names
        """
        
        return self.names[self.long_names.index(long_name.split(' (')[0])]


def GI_netcdf(datapath='.'):
    """
    Function that launches the Graphics Interface for netcdf outputs.
    PARAMETER:
        datapath : string : path where there are the atmym*.nc files.

    """    
    matplotlib.use('TkAgg')
    matplotlib.interactive('t')
    cwd = os.getcwd()
    os.chdir(datapath)
    window = tk.Tk()
    gui = loveclim_GUI(window)
    window.mainloop()
    os.chdir(cwd)

if __name__ == '__main__':
    GI_netcdf()
