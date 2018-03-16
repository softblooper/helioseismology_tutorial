import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.animation as animation
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

import numpy as np
from scipy import ndimage
from pylab import find

from PIL import Image, ImageTk

'''Notes to self:
-I think we're good
'''
#overplotting (maybe?)

#------------------------------------------------------------------------------#

#---Classes---#

class FitsImage(object): # Allows creating of plots with data from fits files. Also used to manage datasets created from within the program, so name is a bit misleading.
    
    # Initializes object, converts fits data into Python array. Detects shape for corresponding plot.
    
    def __init__(self,file,title,colorscale,colorlabel,shortname,dimensions):
        
        self.file = file # File from which data is acquired
        self.title = title # Name to be given to the dataset
        self.color = colorscale # Colorscale from pyplot to be used to see image
        self.colorlabel = colorlabel # Units of what color represents
        self.shortname = shortname # A shorthand name for the dataset (not really used)
        self.dimensions = dimensions # Size of the dataset for user reference (# x #)
        
        # Attempts to use data from a fits file. If it isn't a fits file, then whatever input for the file attribute is given is used as the data itself.
        try:
            self.imgdata = fits.getdata(self.file)
        except:
            self.imgdata = self.file
        
        self.shape = np.shape(self.imgdata) #Shape of data set in (z,y,x) or (y,x)
        self.dim = len(self.shape) # Number of dimensions (2 or 3)

class PowerSpectraAvg(object): # Seperate object for the averaged PS, which requires 3 data sets in order to plot (contour plot)
    
    def __init__(self,data1,data2,data3,title,ratio_x,ratio_y):
        self.data1 = data1
        self.data2 = data2
        self.data3 = data3
        self.title = title
        self.dim = 1 # Abritrary dimension, just to distinguish it from typical 2D and 3D data sets.
        self.xticksmin = np.arange(0, 64, 8)
        self.xticksmax = np.round(np.arange(0, 64, 8)*ratio_x)
        self.yticksmin = np.arange(0, 256, 25)
        self.yticksmax = np.round(np.arange(0, 256, 25)*ratio_y)
        self.dimensions = '(128x128)'

class FitsFiles(object): # Object that will be used as the main data handler. Keeps a list of all data sets, default and created, available to be used.

    def __init__(self):
        
        self.files = [] # List that will contain all of the instances of the other objects. Where actual data is.
        self.options = [] # List that only contains the titles of the datasets.
        self.dataopts = {} # Dictionary of all available objects. Assigns a number to each dataset name for indexing purposes.
    
    def add(self, fitsimage): # Function that allows addition of datasets to the handler.
        
        self.files.append(fitsimage)
        self.options.append(fitsimage.title)
        self.dataopts[fitsimage.title] = len(self.files) - 1

#---Data---#

Data = FitsFiles() # Create instance of data handler.

# Initiliaze all "default" datasets, which are the actual fits files used.

Data.add(FitsImage('fd_Ic_6h_01d.fits','Intensity','gray','Continuum Intensity','intensity','(1024x1024)'))
Data.add(FitsImage('fd_M_96m_01.fits','Magnetogram','gray','Guass (G)','magnetogram','(1024x1024)'))
Data.add(FitsImage('fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)','dopplergram','(1024x1024)'))
Data.add(FitsImage('data1.fits','Data 1','gray','Velocity (m/s)','data1','(128x128x512)'))
Data.add(FitsImage('data2.fits','Data 2','gray','Velocity (m/s)','data2','(128x128x512)'))


#---GUI Modules---#

class StatusBar(tk.Frame): # Status bar located on the top-right side of the window. Used to give some feedback and info to the user. Maybe remove it?
    
    def __init__(self,parent):
        
        tk.Frame.__init__(self,parent)
        
        text = 'Welcome to the Helioseismology Tutorial. Select data and have fun plotting!'
        
        self.statusbar = tk.Label(self,text=text,relief=tk.SUNKEN,bg='white',width=80,anchor='w')
        self.statusbar.grid(row=0)

class OptionMenu(tk.LabelFrame): # Section of sidebar containing main option and buttons to view the other menus.
    
    def __init__(self,parent):
        
        tk.LabelFrame.__init__(self,parent, text = 'Menu')
        
        self.choice = tk.StringVar()
        self.choice.set('(Select a data set)')
        
        self.menu = tk.OptionMenu(self,self.choice,*Data.options)
        self.menu.grid(row = 0, pady = 5)
        self.menu.configure(width = 15)
        
        tabs = tk.Frame(self)
        tabs.grid(row = 1, sticky = tk.E+tk.W)
        
        tabs.columnconfigure(0, weight = 1)
        tabs.columnconfigure(1, weight = 1)
        
        self.showPT = tk.Label(tabs, text = 'Range', relief = tk.RAISED, width = 7)
        self.showPT.grid(row = 0, column = 0)
        
        self.showSlM = tk.Label(tabs, text = 'Slice', relief = tk.RAISED, width = 7)
        self.showSlM.grid(row = 0, column = 1)
        
        self.showScM = tk.Label(tabs, text = 'Scatter', relief = tk.RAISED, width = 7)
        self.showScM.grid(row = 1, column = 1)
        
        self.showAM = tk.Label(tabs, text = 'Animate', relief = tk.RAISED, width = 7)
        self.showAM.grid(row = 2, column = 1)
        
        self.showCM = tk.Label(tabs, text = 'Compute', relief = tk.RAISED, width = 7)
        self.showCM.grid(row = 1, column = 0)
        
    def updateMenu(self): # Function that is used everytime a new data set is created, in order to update dropdown options.
        
        m = self.menu.children['menu']
        m.delete(0, "end")
        for i in range(len(Data.files)):
            m.add_command(label=Data.options[i], command=lambda value=Data.options[i]: self.choice.set(value))

class PlotTools(tk.LabelFrame):
    
    def __init__(self,parent):
        
        tk.LabelFrame.__init__(self,parent, text = 'Plotting Tools')
        
        #Labels for 'Min' and 'Max'
        minlabel = tk.Label(self,text = 'Min')
        minlabel.grid(row=0,column=1)
        maxlabel = tk.Label(self, text = 'Max')
        maxlabel.grid(row=0,column = 2)      
        
        #This label simply marks X as to identify its corresponding entry boxes
        xtools = tk.Label(self,text='X')
        xtools.grid(row = 1, column = 0,padx = 5)
        
        #Entry box for X Minimum
        self.minrangex = tk.Entry(self,width = 6)
        self.minrangex.grid(row=1,column=1,padx=5)
        
        #Entry box for X Maximum
        self.maxrangex = tk.Entry(self,width = 6)
        self.maxrangex.grid(row=1,column=2,padx=5)
        
        #Previous descriptions repeat accordingtly to Y and Z as well
        ytools = tk.Label(self,text='Y')
        ytools.grid(row = 2, column = 0, padx = 5)
        
        self.minrangey = tk.Entry(self,width = 6)
        self.minrangey.grid(row = 2, column = 1, padx = 5)
        
        self.maxrangey = tk.Entry(self,width = 6)
        self.maxrangey.grid(row = 2, column = 2, padx = 5)
        
        ztools = tk.Label(self,text='Z')
        ztools.grid(row = 3, column = 0, padx = 5)
        
        self.minrangez = tk.Entry(self,width = 6)
        self.minrangez.grid(row = 3, column = 1, padx = 5)
        
        self.maxrangez = tk.Entry(self,width = 6)
        self.maxrangez.grid(row = 3, column = 2, padx = 5)
    
    def clearEntries(self):
        
        self.minrangex.delete(0,tk.END)
        self.maxrangex.delete(0,tk.END)
        self.minrangey.delete(0,tk.END)
        self.maxrangey.delete(0,tk.END)
        self.minrangez.delete(0,tk.END)
        self.maxrangez.delete(0,tk.END)

class SliceMenu(tk.LabelFrame):
    
    def __init__(self,parent,command):
        
        tk.LabelFrame.__init__(self,parent,text='Slice Plots')
        
        sliceentries = tk.Frame(self)
        sliceentries.grid(row=0)
        
        xslicelabel = tk.Label(sliceentries,text='X')
        xslicelabel.grid(row=0,column=0,padx=5)
        
        yslicelabel = tk.Label(sliceentries,text='Y')
        yslicelabel.grid(row=0,column=1,padx=5)
        
        tslicelabel = tk.Label(sliceentries,text='t')
        tslicelabel.grid(row=0,column=2,padx=5)
        
        self.xslice = tk.Entry(sliceentries,width = 3)
        self.xslice.grid(row=1,column=0,padx=5)
        
        self.yslice = tk.Entry(sliceentries,width = 3)
        self.yslice.grid(row=1,column=1,padx=5)
        
        self.tslice = tk.Entry(sliceentries,width = 3)
        self.tslice.grid(row=1,column=2,padx=5)
        
        #Checkbox that allows user to open the desired image in a new window.
        self.keep = tk.BooleanVar()
        keepimgopt = tk.Checkbutton(self, text ='Open in new window',variable=self.keep)
        keepimgopt.grid(row=2)
        
        slicebutton = tk.Button(self, text = 'View', command = command)
        slicebutton.grid(row = 3)
    
    def clearEntries(self):
        
        self.tslice.delete(0,tk.END)
        self.xslice.delete(0,tk.END)
        self.yslice.delete(0,tk.END)

class ScatterMenu(tk.LabelFrame):
    
    def __init__(self,parent,command):
        
        tk.LabelFrame.__init__(self,parent,text='Scatter Plots')
        
        #Choose data set to be used for the x-axis
        xaxis = tk.Label(self, text='X-Axis')
        xaxis.grid(row = 0)
        self.choice = tk.Label(self, relief = tk.SUNKEN, text = '(Select a data set)', width = 20)
        self.choice.grid(row = 1)
        
        #Choose data set to be used for the y-axis
        yaxis = tk.Label(self,text='Y-Axis')
        yaxis.grid(row = 2)
        self.datay = tk.StringVar()
        self.datay.set('')
        self.ploty = tk.OptionMenu(self,self.datay,*Data.options)
        self.ploty.grid(row = 3)
        
        #Checkbox that allows user to open the desired image in a new window.
        self.keep = tk.BooleanVar()
        self.keep.set(0)
        keepplotopt = tk.Checkbutton(self, text = "Open in new window",variable=self.keep) #Implement feature
        keepplotopt.grid(row = 4)
        
        plotbutton = tk.Button(self, text = 'Plot', command = command)
        plotbutton.grid(row = 5)
        
        #Button that allows user to plot both data sets in a scatter plat on the canvas.

    def updateMenu(self):
        
        n = self.ploty.children['menu']
        n.delete(0, "end")
        for i in range(len(Data.files)):
            n.add_command(label=Data.options[i], command=lambda value=Data.options[i]: self.datay.set(value))

class AnimationMenu(tk.LabelFrame):
    
    def __init__(self,parent,command):
        
        tk.LabelFrame.__init__(self,parent,text = 'Animate')
        
        frame1 = tk.Frame(self)
        frame1.grid(row = 1, column = 0)
        
        fpstext = tk.Label(frame1, text = 'FPS')
        fpstext.grid(row = 0, column = 0, sticky = tk.S)
        
        self.fps = tk.Entry(frame1,width = 3)
        self.fps.insert(tk.END,'5')
        self.fps.grid(row = 1, column = 0)
        
        axistest = tk.Label(self, text = 'Animation Axis')
        axistest.grid(row = 0, column = 1)
        
        frame2 = tk.Frame(self)
        frame2.grid(row = 1, column = 1)
        
        self.aniaxis = tk.IntVar()
        
        tani = tk.Radiobutton(frame2, text = 't-Axis', variable = self.aniaxis, value = 0)
        xani = tk.Radiobutton(frame2, text = 'X-Axis', variable = self.aniaxis, value = 2)
        yani = tk.Radiobutton(frame2, text = 'Y-Axis', variable = self.aniaxis, value = 1)
        tani.grid(row = 0)
        xani.grid(row = 1)
        yani.grid(row = 2)
        
        animateoption = tk.Button(self, text='Animate',command=command)
        animateoption.grid(row=2,columnspan = 2)
        
        self.columnconfigure(0, weight = 1)
        self.columnconfigure(1, weight = 1)
        self.rowconfigure(0, weight = 1)
        self.rowconfigure(1, weight = 1)
        self.rowconfigure(2, weight = 1)

class ComputationMenu(tk.LabelFrame):
    
    def __init__(self,parent,command):
        
        tk.LabelFrame.__init__(self,parent,text = 'Generate')
        
        #Frame to organize checkbuttons for Average and Difference computations.
        avgdif = tk.Frame(self)
        avgdif.grid(row = 1)
        
        #Checkbutton that allows the user to generate an average between two slices of a 3-D data set.
        self.imgavg = tk.BooleanVar()
        imgavgopt = tk.Checkbutton(avgdif,text='2 Min. Average',variable=self.imgavg)
        imgavgopt.grid(row=0,column=0)
        
        #Checkbutton that allows the user to generate a difference between two slices of a 3-D data set.
        self.imgdif = tk.BooleanVar()#
        imgdifopt = tk.Checkbutton(avgdif,text='2 Min. Difference',variable=self.imgdif)
        imgdifopt.grid(row=1,column=0)
        
        '''temporallabel = tk.Label(self,text = 'Temporal')
        temporallabel.grid(row = 2)'''
        
        temporal = tk.Frame(self)
        temporal.grid(row = 3)
        
        self.tempavg = tk.BooleanVar()#
        tempavgopt = tk.Checkbutton(temporal, text = 'Mean', variable = self.tempavg)
        tempavgopt.grid(row = 0, column = 0)
        
        self.tempdiff = tk.BooleanVar()
        tempdiffopt = tk.Checkbutton(temporal, text = 'Residuals', variable = self.tempdiff)
        tempdiffopt.grid(row = 0, column = 1)
        
        self.var = tk.BooleanVar()#
        varianceopt = tk.Checkbutton(self, text = 'Variance',variable=self.var)
        varianceopt.grid(row = 4)
        
        #Checkbutton that allows the user to generate a power spectra of a data set.
        self.ps = tk.BooleanVar()#
        powerspectraopt = tk.Checkbutton(self, text = 'Power Spectra',variable=self.ps)
        powerspectraopt.grid(row = 5)
        
        #Button that allows the user to compute the selected computations.
        computebutton = tk.Button(self, text='Compute', command=command)
        computebutton.grid(row = 6)

class PlotCanvas(tk.Frame):
    
    def __init__(self,parent):
        
        tk.Frame.__init__(self,parent)
        
        self.f, self.ax = plt.subplots()
        
        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
        
        #Displays Matplotlib figure toolbar.
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class Window(tk.Frame):
    
    def __init__(self,parent):
        
        tk.Frame.__init__(self,parent)
        
        self.statusBar = StatusBar(self)
        self.statusBar.pack(anchor = 'nw', side = tk.TOP, pady = (3,3))
        
        self.canvas = PlotCanvas(self)
        self.canvas.pack(fill=tk.BOTH, expand=True)

class FileInfo(tk.LabelFrame):
    
    def __init__(self,parent):
        
        tk.LabelFrame.__init__(self,parent, text = 'File Info')
        
        self.name = 'Name: '
        self.dim = 'Size: '
        self.file = 'File: '
        
        self.namelabel = tk.Label(self, text = self.name, anchor = tk.W)
        self.namelabel.grid(row = 0, sticky = tk.W)
        
        self.dimlabel = tk.Label(self, text = self.dim, anchor = tk.W)
        self.dimlabel.grid(row = 1, sticky = tk.W)
        
        self.filelabel = tk. Label(self, text = self.file, anchor = tk.W)
        self.filelabel.grid(row = 2, sticky = tk.W)
    
    def update(self,fitsimage):
        self.namelabel.config(text = self.name + ' ' + fitsimage.title, justify = tk.LEFT)
        self.dimlabel.config(text = self.dim + ' ' + fitsimage.dimensions, justify = tk.LEFT)
        try:
            self.filelabel.config(text = self.file + ' ' + fitsimage.file, justify = tk.LEFT)
        except:
            self.filelabel.config(text = self.file + ' Generated Data')

#---Main App---#

class Helioseismology(tk.Tk):
    
    def __init__(self,*args,**kwargs):
        
        tk.Tk.__init__(self,*args,**kwargs)
        
        tk.Tk.wm_title(self,'Helioseismology Tutorial')
        
        self.mainCanvas = PlotCanvas(self)
        self.mainCanvas.grid(row = 1, column = 0, columnspan = 2, rowspan = 7, sticky = tk.W+tk.E+tk.N+tk.S)
        
        sideMenu = tk.Frame()
        sideMenu.grid(row = 1, column = 3,sticky = tk.N+tk.S)
        
        self.optionMenu = OptionMenu(sideMenu)
        self.optionMenu.grid(row = 0, sticky = tk.E+tk.W)
        
        self.computationMenu = ComputationMenu(sideMenu,self.compute)
        self.computationMenu.grid(row = 1, sticky = tk.E+tk.W, ipady = 3)
        self.computationMenu.grid_remove()
        
        self.plotTools = PlotTools(sideMenu)
        self.plotTools.grid(row = 2, sticky = tk.E+tk.W, ipady = 3)
        self.plotTools.grid_remove()
        
        self.sliceMenu = SliceMenu(sideMenu,self.viewslice)
        self.sliceMenu.grid(row = 3, sticky = tk.E+tk.W, ipady = 3)
        
        self.scatterMenu = ScatterMenu(sideMenu,self.plot)
        self.scatterMenu.grid(row = 3, sticky = tk.E+tk.W, ipady = 3)
        self.scatterMenu.grid_remove()
        
        self.animationMenu = AnimationMenu(sideMenu,self.animation)
        self.animationMenu.grid(row = 3, sticky = tk.E+tk.W, ipady = 3)
        self.animationMenu.grid_remove()
        
        self.fileInfo = FileInfo(sideMenu)
        self.fileInfo.grid(row = 4, sticky = tk.E+tk.W, ipady = 3)
        
        self.optionMenu.showSlM.config(relief = tk.SUNKEN)
        self.optionMenu.showPT.bind('<Button-1>',self.showPT)
        self.optionMenu.showSlM.bind('<Button-1>',self.showSlM)
        self.optionMenu.showScM.bind('<Button-1>',self.showScM)
        self.optionMenu.showCM.bind('<Button-1>',self.showCM)
        self.optionMenu.showAM.bind('<Button-1>',self.showAM)
        
        self.columnconfigure(0, weight = 10)
        self.rowconfigure(1, weight = 10)
        self.PT = False
        self.CM = False
        self.optionMenu.choice.trace('w', self.xScatterChange)
        
        nsf = Image.open('nsf1.gif')
        nsf = nsf.resize((96, 96))
        self.nsf = ImageTk.PhotoImage(nsf)
        nsflogo = tk.Label(sideMenu,image = self.nsf)
        nsflogo.grid(row=6,sticky=tk.S)
        
        nmsupng = Image.open('NM_State_logo.png')
        nmsupng = nmsupng.resize((96, 96))
        self.nmsu = ImageTk.PhotoImage(nmsupng)
        nmsulogo = tk.Label(sideMenu,image = self.nmsu)
        nmsulogo.grid(row=7,stick=tk.S)
        
        
        
        self.aspect(1,1,1,1)
        #self.attributes('-zoomed', True)
        #self.minsize(0,800)
        #self.iconbitmap('icon.ico')

    def viewslice(self):
        keep = self.sliceMenu.keep.get()
        try:
            if keep:
                sliceWindow = tk.Toplevel()
                sliceWindow.withdraw()
                
                sliceFrame = Window(sliceWindow)
                sliceFrame.pack(fill=tk.BOTH, expand=True)
            else:
                self.mainCanvas.f.clear()
                plt.figure(1)
            
            index = Data.dataopts[self.optionMenu.choice.get()]
            minx = self.plotTools.minrangex.get()
            maxx = self.plotTools.maxrangex.get()
            miny = self.plotTools.minrangey.get()
            maxy = self.plotTools.maxrangey.get()
            minz = self.plotTools.minrangez.get()
            maxz = self.plotTools.maxrangez.get()
            slicex = self.sliceMenu.xslice.get()
            slicey = self.sliceMenu.yslice.get()
            slicet = self.sliceMenu.tslice.get()
            textx = str(slicex)
            texty = str(slicey)
            textt = str(slicet)
            
            image = Data.files[index]
            
            labels = ('Velocity (m/s)','Time (s)','X-Pix','Y-Pix')
            labelindx = None
            labelindy = None
            
            if (minz and maxz):
                minz = int(minz)
                maxz = int(maxz)
            else:
                minz = None
                maxz = None
            
            if image.dim == 3:
                sizet = np.arange(image.shape[0])
                sizey = np.arange(image.shape[1])
                sizex = np.arange(image.shape[2])
                if slicex:
                    slicex = int(slicex)
                    if slicey:
                        slicey = int(slicey)
                        labelindx = 1
                        slicetext = 'x = '+textx+' and y = '+texty
                        data = image.imgdata[:,slicey,slicex]
                        plt.plot(sizet,data,lw=0.5,color='black')
                    elif slicet:
                        slicet = int(slicet)
                        labelindx = 3
                        slicetext = 't = '+textt+' and x = '+textx
                        data = image.imgdata[slicet,:,slicex]
                        plt.plot(sizey,data,lw=0.5,color='black')
                    else:
                        labelindx = 1
                        labelindy = 3
                        slicetext = 'x = '+textx
                        data = ndimage.rotate(image.imgdata[:,:,slicex], 270)
                        plt.imshow(data, cmap = image.color, vmin=minz, vmax=maxz)
                elif slicey:
                    slicey = int(slicey)
                    if slicet:
                        slicet = int(slicet)
                        labelindx = 2
                        slicetext = 't = '+textt+' and y = '+texty
                        data = image.imgdata[slicet,slicey,:]
                        plt.plot(sizex,data,lw=0.5,color='black')
                    else:
                        labelindx = 2
                        labelindy = 1
                        slicetext = 'y = '+texty
                        data = image.imgdata[:,slicey,:]
                        plt.imshow(data, cmap = image.color, vmin=minz, vmax=maxz,origin = 'lower')
                elif slicet:
                    slicet = int(slicet)
                    labelindx = 2
                    labelindy = 3
                    slicetext = 't = '+textt
                    data = image.imgdata[slicet,...]
                    plt.imshow(data, cmap = image.color, vmin=minz, vmax=maxz)
                else:
                    data = image.imgdata[0]
                    plt.imshow(data, cmap = image.color, vmin=minz, vmax=maxz)
                    labelindx = 2
                    labelindy = 3
                    slicetext = 't = 0'
                
                if (slicex and slicey) or (slicex and slicet) or (slicey and slicet):
                    plt.title(labels[0]+' vs '+labels[labelindx]+' at '+slicetext)
                    plt.ylabel(labels[0])
                    plt.xlabel(labels[labelindx])
                else:
                    plt.title('Slice of '+image.title+' at '+slicetext)
                    plt.xlabel(labels[labelindx])
                    plt.ylabel(labels[labelindy])
                    ibar = plt.colorbar()
                    ibar.set_label(image.colorlabel)
                plt.gca().invert_yaxis()
            
            elif image.dim == 2:
                sizey = np.arange(image.shape[0])
                sizex = np.arange(image.shape[1])
                if slicex:
                    slicex=int(slicex)
                    plt.plot(sizey,image.imgdata[:,slicex],lw=0.5,color='black')
                    labelindx = 3
                    slicetext = 'x ='+textx
                elif slicey:
                    slicey=int(slicey)
                    plt.plot(sizex,image.imgdata[slicey,:],lw=0.5,color='black')
                    labelindx = 2
                    slicetext = 'y ='+texty
                elif slicet:
                    messagebox.showerror('Error','Invalid axis.')
                else:
                    plt.imshow(image.imgdata, cmap = image.color, vmin=minz, vmax=maxz)
                    labelindx = 2
                    labelindy = 3
                if (slicex or slicey):
                    plt.title(image.colorlabel+' at '+slicetext)
                    plt.ylabel(image.colorlabel)
                    plt.xlabel(labels[labelindx])
                else:
                    plt.title(image.title)
                    plt.xlabel(labels[labelindx])
                    plt.ylabel(labels[labelindy])
                    ibar = plt.colorbar()
                    ibar.set_label(image.colorlabel)
                plt.gca().invert_yaxis()
            
            else:
                plt.title(image.title)
                plt.xlabel('Wavenumber (1/Mm)')
                plt.ylabel('Frequency (mHz)')
                #plt.xticks(image.xticksmin, image.xticksmax)
                #plt.yticks(image.yticksmin, image.yticksmax)
                plt.contourf(image.data1, image.data2, image.data3, 100)
                ibar = plt.colorbar()
                ibar.set_label('Power')
        
            if (minx and maxx):
                minx = int(minx)
                maxx = int(maxx)
                plt.xlim(minx,maxx)
            if (miny and maxy):
                miny = int(miny)
                maxy = int(maxy)
                plt.ylim(maxy,miny)
            if image.dim == 1:
                text = 'Showing: ' + image.title
            else:
                text = "Showing: "+image.title+'. Size: '+image.dimensions
            
            if keep:
                sliceFrame.statusBar.statusbar.config(text = text)
                sliceWindow.deiconify()
            else:
                self.mainCanvas.f.canvas.draw()
            
            self.plotTools.clearEntries()
            self.sliceMenu.clearEntries()
            self.fileInfo.update(image)
        
        except:
            messagebox.showerror('Error','Select a dataset.')
            if keep:
                sliceWindow.destroy()
    
    def animation(self):
        #keep = self.sliceMenu.keep.get()
        try:
            aniWindow = tk.Toplevel()
            aniWindow.withdraw()
            
            aniFrame = Window(aniWindow)
            aniFrame.pack(fill=tk.BOTH, expand=True)
            
            index = Data.dataopts[self.optionMenu.choice.get()]
            minx = self.plotTools.minrangex.get()
            maxx = self.plotTools.maxrangex.get()
            miny = self.plotTools.minrangey.get()
            maxy = self.plotTools.maxrangey.get()
            minz = self.plotTools.minrangez.get()
            maxz = self.plotTools.maxrangez.get()
            
            labels = ('Time (s)','X-Pix','Y-Pix')
            labelindx = None
            labelindy = None
            
            aniaxis = int(self.animationMenu.aniaxis.get())
            image = Data.files[index]
            length = image.shape[aniaxis]
            try:
                speed = 1000/(int(self.animationMenu.fps.get()))
            except:
                self.statusBar.statusbar.config(text = 'Error! FPS must be an integer')
            
            if (minz and maxz):
                minz = int(minz)
                maxz = int(maxz)
            else:
                minz = None
                maxz = None
            
            if aniaxis == 0:
                anitext = 't'
                labelindx = 1
                labelindy = 2
                im = plt.imshow(image.imgdata[0], cmap = image.color, vmin=minz, vmax=maxz)
            elif aniaxis == 2:
                anitext = 'X'
                labelindx = 0
                labelindy = 2
                data = ndimage.rotate(image.imgdata[:,0,:], 270)
                im = plt.imshow(data, cmap = image.color, vmin=minz, vmax=maxz)
            elif aniaxis == 1:
                anitext = 'Y'
                labelindx = 1
                labelindy = 0
                im = plt.imshow(image.imgdata[...,0], cmap = image.color, vmin=minz, vmax=maxz, origin = 'lower')
            
            if (minx and maxx):
                minx = int(minx)
                maxx = int(maxx)
                plt.xlim(minx,maxx)
            if (miny and maxy):
                miny = int(miny)
                maxy = int(maxy)
                plt.ylim(maxy,miny)
            
            title = plt.title('Animation of '+image.title+' through '+anitext+'-Axis')
            plt.xlabel(labels[labelindx])
            plt.ylabel(labels[labelindy])
            ibar = plt.colorbar()
            ibar.set_label(image.colorlabel)
            
            def ani(i):
                frame = i + 1
                if aniaxis == 0:
                    im.set_array(image.imgdata[i])
                    title.set_text('Time Animation, t = '+str(frame))
                    return [im]
                    return title
                if aniaxis == 2:
                    d = ndimage.rotate(image.imgdata[:,i,:], 270)
                    im.set_array(d)
                    title.set_text('X-Axis Animation, x = '+str(frame))
                    return [im]
                    return title
                if aniaxis == 1:
                    im.set_array(image.imgdata[...,i])
                    title.set_text('Y-Axis Animation, y = '+str(frame))
                    return [im]
                    return title
            
            aniFrame.statusBar.statusbar.config(text = 'Animation of '+image.title+' through '+anitext+'-Axis')
            aniWindow.deiconify()
            self.plotTools.clearEntries()
            
            ani = animation.FuncAnimation(aniFrame.canvas.f,ani,frames=range(length),interval = speed)
            ani._start()
        except:
            messagebox.showerror('Error',"Can't animate this dataset.")
            
            aniWindow.destroy()
    
    def plot(self):
        keep = self.scatterMenu.keep.get()
        
        try:
            if keep:
                plotWindow = tk.Toplevel()
                plotWindow.withdraw()
                
                plotFrame = Window(plotWindow)
                plotFrame.pack()
            
            else:
                self.mainCanvas.f.clear()
                plt.figure(1)
            
            index1 = Data.dataopts[self.optionMenu.choice.get()]
            index2 = Data.dataopts[self.scatterMenu.datay.get()]
            minx = self.plotTools.minrangex.get()
            maxx = self.plotTools.maxrangex.get()
            miny = self.plotTools.minrangey.get()
            maxy = self.plotTools.maxrangey.get()
            data1 = Data.files[index1]
            data2 = Data.files[index2]
        
            if (data1.dim == 3) and (data2.dim == 3):
                plt.plot(data1.imgdata[0],data2.imgdata[0],',',color='black')
            else:
                plt.plot(data1.imgdata,data2.imgdata,',',color='black')
            plt.xlabel(data1.colorlabel)
            plt.ylabel(data2.colorlabel)
            if (minx and maxx):
                minx = int(minx)
                maxx = int(maxx)
                plt.xlim(minx,maxx)
            if (miny and maxy):
                miny = int(miny)
                maxy = int(maxy)
                plt.ylim(miny,maxy)
            
            text = "Showing plot of "+data1.title+" against "+data2.title
            
            if keep:
                plotFrame.statusBar.statusbar.config(text = text)
                plotWindow.deiconify()
            else:
                self.mainCanvas.f.canvas.draw()
            
            self.plotTools.clearEntries()
            
        except Exception:
            messagebox.showerror('Error',"Can't plot different sized datasets.")
            if keep:
                plotWindow.destroy()
    
    def imgaverage(self):
        index = Data.dataopts[self.optionMenu.choice.get()]
        if Data.files[index].dim == 3:
            average = (Data.files[index].imgdata[0] + Data.files[index].imgdata[2])/2
            title = Data.files[index].title + ' Average'
            color = Data.files[index].color
            colorlabel = Data.files[index].colorlabel
            shortname = 'avg'+Data.files[index].shortname
            dimensions = Data.files[index].dimensions[0:8] + ')'
            Data.add(FitsImage(average,title,color,colorlabel,shortname,dimensions))
        else:
            messagebox.showerror('Error',"Can't compute with this dataset!")
    
    def imgdifference(self):
        index = Data.dataopts[self.optionMenu.choice.get()]
        if Data.files[index].dim == 3:
            difference = Data.files[index].imgdata[2] - Data.files[index].imgdata[1]
            title = Data.files[index].title + ' Difference'
            color = Data.files[index].color
            colorlabel = Data.files[index].colorlabel
            shortname = 'dif'+Data.files[index].shortname
            dimensions = Data.files[index].dimensions[0:8] + ')'
            Data.add(FitsImage(difference,title,color,colorlabel,shortname,dimensions))
        else:
            messagebox.showerror('Error',"Can't compute with this dataset!")
    
    def temporalavg(self):
        index = Data.dataopts[self.optionMenu.choice.get()]
        if Data.files[index].dim == 3:
            tempavg = np.mean(Data.files[index].imgdata,axis=0)
            title = Data.files[index].title + ' Mean'
            color = Data.files[index].color
            colorlabel = Data.files[index].colorlabel
            shortname = 'tempavg'+Data.files[index].shortname
            dimensions = Data.files[index].dimensions[0:8] + ')'
            Data.add(FitsImage(tempavg,title,color,colorlabel,shortname,dimensions))
        else:
            messagebox.showerror('Error',"Can't compute with this dataset!")
    
    def temporaldiff(self):
        index = Data.dataopts[self.optionMenu.choice.get()]
        if Data.files[index].dim == 3:
            tempavg = np.mean(Data.files[index].imgdata,axis=0)
            length = Data.files[index].shape[0]
            tempdiff = []
            for i in range(length):
                td = Data.files[index].imgdata[i] - tempavg
                tempdiff.append(td)
            title = Data.files[index].title + ' Residuals'
            color = Data.files[index].color
            colorlabel = Data.files[index].colorlabel
            shortname = 'residual'+Data.files[index].shortname
            dimensions = Data.files[index].dimensions
            Data.add(FitsImage(tempdiff,title,color,colorlabel,shortname,dimensions))
        else:
            messagebox.showerror('Error',"Can't compute with this dataset!")
    
    def powerspectra(self):
        index = Data.dataopts[self.optionMenu.choice.get()]
        try:
            ofrq1 = np.fft.fftn(Data.files[index].imgdata)
            frq1 = np.fft.fftshift(ofrq1)
            power = np.log(np.abs(frq1**2))
            
            title = Data.files[index].title + ' Power Spectra'
            
            #PS ends here.
            
            Data.add(FitsImage(power,title,None,'Power',title+'power','(128x128x512)'))
            
            begin_space = 0
            end_time = power.shape[0]
            end_space = power.shape[1]
            cen_time = power.shape[0]/2
            cen_space = power.shape[1]/2
            
            x = np.linspace(-cen_space, cen_space-1, end_space)
            y = np.linspace(-cen_space, cen_space-1, end_space)
            X,Y = np.meshgrid(x,y)
            dist = np.hypot(X,Y)
            
            cen_time = int(cen_time)
            cen_space = int(cen_space)
            a = np.zeros([cen_time, cen_space])
            w = 0.5
            
            for i in range(cen_time, end_time):
                pp = power[i,...]
                for r in range(begin_space, cen_space):
                    inds = find(np.logical_and(dist > r-w, dist < r+w))
                    flatpower = pp.flatten()
                    avg = np.mean(flatpower[inds])
                    a[i-cen_time, r] = a[i-cen_time, r] +avg
            
            dx = 1.39 # length per pixel in Mm
            pix = 128 # number of pixels
            kx = np.pi/dx
            r_sun = 696 # radius of the sun in Mm
            l = kx * r_sun # wavenumber in unit of degree
            ratio_x = kx/64 #l/(pix/2)
            
            dt = 60 # in sec
            omega = np.pi/dt # temporal frquency
            v = omega/(2*np.pi) # cyclic frequency
            frq = v*1000 # in unit of mHz
            ratio_y = frq/256
            
            m = np.linspace(0, kx, 64)
            n = np.linspace(0, frq, 256)
            M, N = np.meshgrid(m, n)
            A = a
            
            title = title + ' Avg.'
            
            Data.add(PowerSpectraAvg(M,N,A,title,ratio_x,ratio_y)) #Power Spectrum Average
                
        except:
            messagebox.showerror('Error',"Can't calculate power spectrum! Wrong dataset.")
    
    def variance(self):
        index = Data.dataopts[self.optionMenu.choice.get()]
        if Data.files[index].dim == 3:
            tempavg = np.mean(Data.files[index].imgdata,axis=0)
            length = Data.files[index].shape[0]
            tempdiff = []
            for i in range(length):
                td = Data.files[index].imgdata[i] -tempavg
                tempdiff.append(td)
            variance = (1/len(tempdiff)) * (np.sum(np.array(tempdiff)**2, axis = 0))
            title = Data.files[index].title + ' Variance'
            color = Data.files[index].color
            colorlabel = r'Velocity$^2$ ($m^2/s^2$)'
            shortname = 'variance'+Data.files[index].shortname
            dimensions = Data.files[index].dimensions[0:8] + ')'
            Data.add(FitsImage(variance,title,color,colorlabel,shortname,dimensions))
        else:
            messagebox.showerror('Error',"Can't compute with this dataset!")
    
    def compute(self):
        avg = self.computationMenu.imgavg.get()
        dif = self.computationMenu.imgdif.get()
        tempavg = self.computationMenu.tempavg.get()
        tempdiff = self.computationMenu.tempdiff.get()
        ps = self.computationMenu.ps.get()
        var = self.computationMenu.var.get()
        
        if avg:
            self.imgaverage()
        if dif:
            self.imgdifference()
        if ps:
            self.powerspectra()
        if tempavg:
            self.temporalavg()
        if tempdiff:
            self.temporaldiff()
        if var:
            self.variance()
        
        self.optionMenu.updateMenu()
        self.scatterMenu.updateMenu()
    
    def xScatterChange(self,*args):
        self.scatterMenu.choice.config(text = self.optionMenu.choice.get())
    
    def showPT(self,event):
        if not self.PT:
            self.PT = True
            self.optionMenu.showPT.config(relief = tk.SUNKEN)
            self.plotTools.grid()
        else:
            self.PT = False
            self.optionMenu.showPT.config(relief = tk.RAISED)
            self.plotTools.grid_remove()            
    
    def showSlM(self,event):
        self.optionMenu.showSlM.config(relief = tk.SUNKEN)
        self.optionMenu.showScM.config(relief = tk.RAISED)
        self.optionMenu.showAM.config(relief = tk.RAISED)
        self.sliceMenu.grid()
        self.scatterMenu.grid_remove() 
        self.animationMenu.grid_remove()
    
    def showScM(self,event):
        self.optionMenu.showSlM.config(relief = tk.RAISED)
        self.optionMenu.showScM.config(relief = tk.SUNKEN)
        self.optionMenu.showAM.config(relief = tk.RAISED)
        self.sliceMenu.grid_remove()
        self.scatterMenu.grid() 
        self.animationMenu.grid_remove()
    
    def showCM(self,event):
        if not self.CM:
            self.CM = True
            self.optionMenu.showCM.config(relief = tk.SUNKEN)
            self.computationMenu.grid()
        else:
            self.CM = False
            self.optionMenu.showCM.config(relief = tk.RAISED)
            self.computationMenu.grid_remove() 
    
    def showAM(self,event):
        self.optionMenu.showSlM.config(relief = tk.RAISED)
        self.optionMenu.showScM.config(relief = tk.RAISED)
        self.optionMenu.showAM.config(relief = tk.SUNKEN)
        self.sliceMenu.grid_remove()
        self.scatterMenu.grid_remove() 
        self.animationMenu.grid()
    
app = Helioseismology()
app.mainloop()