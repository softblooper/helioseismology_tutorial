import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import style
import matplotlib.ticker as mticker

import tkinter as tk
from tkinter import ttk

from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

import numpy as np
from scipy import ndimage
from pylab import find

#------------------------------------------------------------------------------#
#--Tutorial-Required Classes

class FitsImage(object): #Allows creating of plots with data from fits files
    
    #Initializes object, converts fits data into Python array. Detects shape for corresponding plot.
    def __init__(self,file,title,colorscale,colorlabel,shortname,dimensions):
        self.file = file
        self.title = title
        self.color = colorscale
        self.colorlabel = colorlabel
        self.shortname = shortname
        self.dimensions = dimensions
        
        try:
            self.imgdata = fits.getdata(self.file)
        except:
            self.imgdata = self.file
        
        self.shape = np.shape(self.imgdata)
        self.dim = len(self.shape)
    
    #Animation Method: Allows animation of image throughout time. NOTE: Must state variable when used, eg: a = Example.animateplot(...).
    def animateplot(self,numofslices,rowdim,coldim,row,col,rowspan,colspan): #Fix animations. Not working?
        plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
        plt.title(self.title)
        plt.gca().invert_yaxis()
        plt.xticks([])
        plt.yticks([])
        ims = []
        for i in range(numofslices):
            im = plt.imshow(self.imgdata[i], cmap = self.color)
            ims.append([im])
        return ims
    
    def plot(self,rowdim,coldim,row,col,rowspan,colspan):
            #plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
            plt.figure
            plt.title(self.title)
            if self.dim == 3:
                plt.imshow(self.imgdata[0], cmap = self.color)
            else:
                plt.imshow(self.imgdata, cmap = self.color)
            ibar = plt.colorbar()
            ibar.set_label(self.colorlabel)
            plt.gca().invert_yaxis()
            plt.xticks([])
            plt.yticks([])
            plt.grid(False)
    
class PowerSpectra(object):
    
    def __init__(self,data1,data2,data3,title,ratio_x,ratio_y):
        self.data1 = data1
        self.data2 = data2
        self.data3 = data3
        self.title = title
        self.dim = 1
        self.xticksmin = np.arange(0, 64, 8)
        self.xticksmax = np.round(np.arange(0, 64, 8)*ratio_x)
        self.yticksmin = np.arange(0, 256, 25)
        self.yticksmax = np.round(np.arange(0, 256, 25)*ratio_y)

class FitsFiles(object):

    def __init__(self):
        
        self.files = []
        self.options = []
        self.dataopts = {}
    
    def add(self, fitsimage):
        
        self.files.append(fitsimage)
        self.options.append(fitsimage.title)
        self.dataopts[fitsimage.title] = len(self.files) - 1

#------------------------------------------------------------------------------#
#--Tutorials--
Intensity = FitsImage('fd_Ic_6h_01d.fits','Intensity','gray','Continuum Intensity','intensity','(1024x1024)')
Magnetogram = FitsImage('fd_M_96m_01.fits','Magnetogram','gray','Guass (G)','magnetogram','(1024x1024)')
Dopplergram = FitsImage('fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)','dopplergram','(1024x1024)')
Data1 = FitsImage('data1.fits','Data 1','gray','Velocity (m/s)','data1','(128x128x512)')
Data2 = FitsImage('data2.fits','Data 2','gray','Velocity (m/s)','data2','(128x128x512)')

Data = [Intensity, Magnetogram, Dopplergram, Data1, Data2] #D1vD2, Average1, Difference1,Average2,Difference2)

DATA = FitsFiles()

DATA.add(Intensity)
DATA.add(Magnetogram)
DATA.add(Dopplergram)
DATA.add(Data1)
DATA.add(Data2)

Options = ['Intensity', 'Magnetogram', 'Dopplergram', 'Data 1', 'Data 2']
DataIndex = [0, 1, 2, 3, 4]

DataOpts = dict(zip(Options,DataIndex))

#------------------------------------------------------------------------------#
#--GUI--

'''class RangeLabels(tk.Frame): #Revise this. Frame has to be argument.

    def __init__(self,frame,*args,**kwargs):

        tk.Frame.__init__(self,*args,**kwargs)

        self.minlabel = tk.Label(frame,text = 'Min')
        self.minlabel.grid(row=0,column=0)
        self.maxlabel = tk.Label(frame, text = 'Max')
        self.maxlabel.grid(row=0,column=1)'''

class Helioseismology(tk.Tk):
    
    def __init__(self,*args,**kwargs):
        
        tk.Tk.__init__(self,*args,**kwargs)
        
        tk.Tk.wm_title(self,'Helioseismology Tutorial')
        
        #Scatter plot function
        def plot():
            keep = self.keepplot.get()
            if keep:
                g = plt.figure(2)
                g.clear()
            else:
                self.f.clear()
                plt.figure(1)
            index1 = DATA.dataopts[datax.get()]
            index2 = DATA.dataopts[datay.get()]
            minx = minrangex.get()
            maxx = maxrangex.get()
            miny = minrangey.get()
            maxy = maxrangey.get()
            try:
                if (Data[index1].dim == 3) and (Data[index2].dim == 3):
                    plt.plot(Data[index1].imgdata[0],Data[index2].imgdata[0],',',color='black')
                else:
                    plt.plot(Data[index1].imgdata,Data[index2].imgdata,',',color='black')
                plt.xlabel(Data[index1].colorlabel)
                plt.ylabel(Data[index2].colorlabel)
                if (minx and maxx):
                    minx = int(minx)
                    maxx = int(maxx)
                    plt.xlim(minx,maxx)
                if (miny and maxy):
                    miny = int(miny)
                    maxy = int(maxy)
                    plt.ylim(miny,maxy)
                if keep:
                    kept = tk.Toplevel()
                    
                    text="Showing plot of "+Data[index1].title+" against "+Data[index2].title
                    self.statusbar = tk.Label(kept,text=text,relief=tk.SUNKEN,bg='white',width=60,anchor='w')
                    self.statusbar.pack(anchor='nw',side=tk.TOP,pady=(3,3))
                    
                    keptcanvas = FigureCanvasTkAgg(g, kept)
                    keptcanvas.show()
                    keptcanvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
                    
                    kepttoolbar = NavigationToolbar2TkAgg(keptcanvas, kept)
                    kepttoolbar.update()
                    keptcanvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
                else:
                    self.f.canvas.draw()
                    self.statusbar.config(text="Showing plot of "+Data[index1].title+" against "+Data[index2].title)
                minrangex.delete(0,tk.END)
                maxrangex.delete(0,tk.END)
                minrangey.delete(0,tk.END)
                maxrangey.delete(0,tk.END)
            except Exception:
                self.statusbar.config(text="Can't plot these data sets! Differente sizes.")
        
        def viewslice():
            keep = self.keepimg.get()
            if keep:
                g = plt.figure(2)
                g.clear()
            else:
                self.f.clear()
                plt.figure(1)
            index = DATA.dataopts[self.slicechoice.get()]
            minx = minrangex.get()
            maxx = maxrangex.get()
            miny = minrangey.get()
            maxy = maxrangey.get()
            minz = minrangez.get()
            maxz = maxrangez.get()
            slicex = xslice.get()
            slicey = yslice.get()
            slicet = tslice.get()
            
            image = DATA.files[index]
            
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
                        slicetext = 'x = '+xslice.get()+' and y = '+yslice.get()
                        data = image.imgdata[:,slicey,slicex]
                        plt.plot(sizet,data,lw=0.5,color='black')
                    elif slicet:
                        slicet = int(slicet)
                        labelindx = 3
                        slicetext = 't = '+tslice.get()+' and x = '+xslice.get()
                        data = image.imgdata[slicet,:,slicex]
                        plt.plot(sizey,data,lw=0.5,color='black')
                    else:
                        labelindx = 1
                        labelindy = 3
                        slicetext = 'x = '+xslice.get()
                        data = ndimage.rotate(image.imgdata[:,:,slicex], 270)
                        plt.imshow(data, cmap = image.color, vmin=minz, vmax=maxz)
                elif slicey:
                    slicey = int(slicey)
                    if slicet:
                        slicet = int(slicet)
                        labelindx = 2
                        slicetext = 't = '+tslice.get()+' and y = '+yslice.get()
                        data = image.imgdata[slicet,slicey,:]
                        plt.plot(sizex,data,lw=0.5,color='black')
                    else:
                        labelindx = 2
                        labelindy = 1
                        slicetext = 'y = '+yslice.get()
                        data = image.imgdata[:,slicey,:]
                        plt.imshow(data, cmap = image.color, vmin=minz, vmax=maxz)
                elif slicet:
                    slicet = int(slicet)
                    labelindx = 2
                    labelindy = 3
                    slicetext = 't = '+tslice.get()
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
            
            elif image.dim == 2:
                sizey = np.arange(image.shape[0])
                sizex = np.arange(image.shape[1])
                if slicex:
                    slicex=int(slicex)
                    plt.plot(sizey,image.imgdata[:,slicex],lw=0.5,color='black')
                    labelindx = 3
                    slicetext = 'x ='+xslice.get()
                elif slicey:
                    slicey=int(slicey)
                    plt.plot(sizex,image.imgdata[slicey,:],lw=0.5,color='black')
                    labelindx = 2
                    slicetext = 'y ='+yslice.get()
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
            
            else:
                plt.xlabel('wavenumber l')
                plt.ylabel('frequency (mHz)')
                plt.xticks(image.xticksmin, image.xticksmax)
                plt.yticks(image.yticksmin, image.yticksmax)
                plt.contourf(image.data1, image.data2, image.data3, 100)
        
            if (minx and maxx):
                minx = int(minx)
                maxx = int(maxx)
                plt.xlim(minx,maxx)
            if (miny and maxy):
                miny = int(miny)
                maxy = int(maxy)
                plt.ylim(maxy,miny)
            
            if keep:
                kept = tk.Toplevel()
                
                if image.dim == 1:
                    text = 'Showing: '
                else:
                    text="Showing: "+image.title+'. Size: '+image.dimensions
                self.statusbar = tk.Label(kept,text=text,relief=tk.SUNKEN,bg='white',width=60,anchor='w')
                self.statusbar.pack(anchor='nw',side=tk.TOP,pady=(3,3))
                
                keptcanvas = FigureCanvasTkAgg(g, kept)
                keptcanvas.show()
                keptcanvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,expand=True)
                
                kepttoolbar = NavigationToolbar2TkAgg(keptcanvas, kept)
                kepttoolbar.update()
                keptcanvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            else:
                self.f.canvas.draw()
                if image.dim == 1:
                    self.statusbar.config(text = 'Showing: ' + image.title)
                else:
                    self.statusbar.config(text="Showing: "+image.title+'. Size: '+image.dimensions)
            
            minrangex.delete(0,tk.END)
            maxrangex.delete(0,tk.END)
            minrangey.delete(0,tk.END)
            maxrangey.delete(0,tk.END)
            minrangez.delete(0,tk.END)
            maxrangez.delete(0,tk.END)
            tslice.delete(0,tk.END)
            xslice.delete(0,tk.END)
            yslice.delete(0,tk.END)
        
        '''def animate(): #FIX
            h = plt.figure(3)
            index = DATA.dataopts[self.slicechoice.get()]
            minz = minrangez.get()
            maxz = maxrangez.get()
            if DATA.files[index].dim == 3:
                ims = []
                for i in range(DATA.files[index].shape[0]):
                    im = plt.imshow(DATA.files[index].imgdata[i], cmap = DATA.files[index].color, vmin=minz, vmax=maxz)
                    ims.append([im])
                plt.title(DATA.files[index].title)
                
                Ani = animation.ArtistAnimation(h, ims, interval=100, blit=True,
                    repeat_delay=1000)
                
                aniwin = tk.Toplevel()
                anicanvas = FigureCanvasTkAgg(h,aniwin)
                anicanvas.show()
                anicanvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,expand=True)
                
                anitoolbar = NavigationToolbar2TkAgg(anicanvas, aniwin)
                anitoolbar.update()
                anicanvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)'''
        
        def imgaverage(): #IT LIVES!!!
            index = DATA.dataopts[cmpchoice.get()]
            average = (DATA.files[index].imgdata[0] + DATA.files[index].imgdata[2])/2
            title = DATA.files[index].title + ' Average'
            color = DATA.files[index].color
            colorlabel = DATA.files[index].colorlabel
            shortname = 'avg'+DATA.files[index].shortname
            dimensions = DATA.files[index].dimensions
            DATA.add(FitsImage(average,title,color,colorlabel,shortname,dimensions))
            #updatemenu()
        
        def imgdifference():
            index = DATA.dataopts[cmpchoice.get()]
            difference = DATA.files[index].imgdata[2] - DATA.files[index].imgdata[1]
            title = DATA.files[index].title + ' Difference'
            color = DATA.files[index].color
            colorlabel = DATA.files[index].colorlabel
            shortname = 'dif'+DATA.files[index].shortname
            dimensions = DATA.files[index].dimensions
            DATA.add(FitsImage(difference,title,color,colorlabel,shortname,dimensions))
            #updatemenu()
        
        def powerspectra():
            index = DATA.dataopts[cmpchoice.get()]
            
            ofrq1 = np.fft.fftn(DATA.files[index].imgdata)
            frq1 = np.fft.fftshift(ofrq1)
            power = np.log(np.abs(frq1**2))
            
            begin_time = 0
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
            
            m = np.linspace(1, cen_space, cen_space)
            n = np.linspace(1, cen_time, cen_time)
            M, N = np.meshgrid(m, n)
            A = a
            
            dx = 1.39 # length per pixel in Mm
            pix = 128 # number of pixels
            kx = ky = np.pi/dx
            r_sun = 696 # radius of the sun in Mm
            l = kx * r_sun # wavenumber in unit of degree
            ratio_x = l/(pix/2)
            
            dt = 60 # in sec
            omega = np.pi/dt # temporal frquency
            v = omega/(2*np.pi) # cyclic frequency
            frq = v*1000 # in unit of mHz
            ratio_y = frq/256
            
            title = DATA.files[index].title + ' Power Spectra'
            
            DATA.add(PowerSpectra(M,N,A,title,ratio_x,ratio_y))
            #updatemenu()
        
        def updatemenu():
            m = self.slicemenu.children['menu']
            m.delete(0, "end")
            for i in range(len(DATA.files)):
                m.add_command(label=DATA.options[i], command=lambda value=DATA.options[i]: self.slicechoice.set(value))
        
        def compute():
            avg = self.imgavg.get()
            dif = self.imgdif.get()
            ps = self.ps.get()
            
            if avg:
                imgaverage()
            if dif:
                imgdifference()
            if ps:
                powerspectra()
            
            updatemenu()
        
        #This is the main frame. It contains and organizes all other frames in the window.
        container=tk.Frame(self)
        container.pack(side='top',fill='both',expand=True)
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0,weight=1)
        
        #This frame is the top one, holding the status bar that follows.
        menu = tk.Frame(container)
        menu.pack(side=tk.TOP, fill=tk.X)
        
        #Status bar at top right of window. Gives feedback to user.
        self.statusbar = tk.Label(menu,text='Welcome to the Helioseismology Tutorial. Select data and have fun plotting!',relief=tk.SUNKEN,bg='white',width=60,anchor='w')
        self.statusbar.pack(anchor='nw',side=tk.LEFT,pady=(3,3))
        
        #This frame is the frame containing the canvas and all plotting tools.
        frame = tk.Frame(container)
        frame.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
        
        #This frame is the one used to organize plotting tools.
        sideframe = tk.Frame(frame)
        sideframe.pack(side=tk.RIGHT,anchor='n')
        
        #--Menu Frames--#
        #All the main frames in order, as for ease of re-arrangement.
        
        #Tools for changing display of image or plot, such as setting X Y Z limits.
        toolsframe = tk.LabelFrame(sideframe,text='Plot Tools')
        toolsframe.pack(side=tk.TOP,fill=tk.X)
        
        #Frame for slice plots
        sliceframe = tk.LabelFrame(sideframe, text='Slice Plots')
        sliceframe.pack(side=tk.TOP,fill=tk.X)
        
        #Frame for scatter plot options
        scatterframe = tk.LabelFrame(sideframe, text="Scatter Plots")
        scatterframe.pack(side=tk.TOP,fill=tk.X)
        
        #Frame for computation of various functions to FITS files.
        computeframe = tk.LabelFrame(sideframe, text='Generate')
        computeframe.pack(side=tk.TOP,fill=tk.X)
        
        #--Tools Frame Widgets--#
        
        #This label simply marks X as to identify its corresponding entry boxes
        xtools = tk.Label(toolsframe,text='X')
        xtools.pack(side=tk.TOP)
        
        #Small frame to neatly organize min and max entry boxes and labels
        rangeframex = tk.Frame(toolsframe)
        rangeframex.pack(side=tk.TOP)
        
        #Labels for 'Min' and 'Max'
        minlabelx = tk.Label(rangeframex,text = 'Min')
        minlabelx.grid(row=0,column=0)
        maxlabelx = tk.Label(rangeframex, text = 'Max')
        maxlabelx.grid(row=0,column=1)
        
        #Entry box for X Minimum
        minrangex = tk.Entry(rangeframex,width = 6)
        minrangex.grid(row=1,column=0,padx=5)
        
        #Entry box for X Maximum
        maxrangex = tk.Entry(rangeframex,width = 6)
        maxrangex.grid(row=1,column=1,padx=5)
        
        #Previous descriptions repeat accordingtly to Y and Z as well
        ytools = tk.Label(toolsframe,text='Y')
        ytools.pack(side=tk.TOP)
        
        rangeframey = tk.Frame(toolsframe)
        rangeframey.pack(side=tk.TOP)
        
        minlabely = tk.Label(rangeframey,text = 'Min')
        minlabely.grid(row=0,column=0)
        maxlabely = tk.Label(rangeframey, text = 'Max')
        maxlabely.grid(row=0,column=1)
        
        minrangey = tk.Entry(rangeframey,width = 6)
        minrangey.grid(row=1,column=0,padx=5)
        
        maxrangey = tk.Entry(rangeframey,width = 6)
        maxrangey.grid(row=1,column=1,padx=5)
        
        ztools = tk.Label(toolsframe,text='Z')
        ztools.pack(side=tk.TOP)
        
        rangeframez = tk.Frame(toolsframe)
        rangeframez.pack(side=tk.TOP)
        
        minlabelz = tk.Label(rangeframez,text = 'Min')
        minlabelz.grid(row=0,column=0)
        maxlabelz = tk.Label(rangeframez, text = 'Max')
        maxlabelz.grid(row=0,column=1)
        
        minrangez = tk.Entry(rangeframez,width = 6)
        minrangez.grid(row=2,column=0,padx=5)
        
        maxrangez = tk.Entry(rangeframez,width = 6)
        maxrangez.grid(row=2,column=1,padx=5,pady=(0,5))
        
        #--Slice Frame Widgets--#
        
        self.slicechoice = tk.StringVar()
        self.slicechoice.set('')
        self.slicemenu = tk.OptionMenu(sliceframe,self.slicechoice,*DATA.options)
        self.slicemenu.pack(side=tk.TOP)
        
        sliceentries = tk.Frame(sliceframe)
        sliceentries.pack(side=tk.TOP)
        
        xslicelabel = tk.Label(sliceentries,text='X')
        xslicelabel.grid(row=0,column=0,padx=5)
        
        yslicelabel = tk.Label(sliceentries,text='Y')
        yslicelabel.grid(row=0,column=1,padx=5)
        
        tslicelabel = tk.Label(sliceentries,text='t')
        tslicelabel.grid(row=0,column=2,padx=5)
        
        xslice = tk.Entry(sliceentries,width = 3)
        xslice.grid(row=1,column=0,padx=5)
        
        yslice = tk.Entry(sliceentries,width = 3)
        yslice.grid(row=1,column=1,padx=5)
        
        tslice = tk.Entry(sliceentries,width = 3)
        tslice.grid(row=1,column=2,padx=5)
        
        #Checkbox that allows user to open the desired image in a new window.
        self.keepimg = tk.BooleanVar()
        keepimgopt = tk.Checkbutton(sliceframe, text = "Open in new window",variable=self.keepimg)
        keepimgopt.pack(side=tk.TOP)
        
        #self.animate = tk.BooleanVar()
        animateoption = tk.Button(sliceframe, text='Animate',command=animate)
        animateoption.pack(side=tk.TOP,pady=(0,5))
        
        slicebutton = tk.Button(sliceframe, text='Plot',command=viewslice)
        slicebutton.pack(side=tk.TOP,pady=(0,5))
        
        #--Scatter Plot Frame--#
        
        #Choose data set to be used for the x-axis
        xaxis = tk.Label(scatterframe,text='X-Axis')
        xaxis.pack(side=tk.TOP)
        datax = tk.StringVar()
        datax.set('')
        plotx = tk.OptionMenu(scatterframe,datax,*DATA.options)
        plotx.pack(side=tk.TOP)
        
        #Choose data set to be used for the y-axis
        yaxis = tk.Label(scatterframe,text='Y-Axis')
        yaxis.pack(side=tk.TOP)
        datay = tk.StringVar()
        datay.set('')
        ploty = tk.OptionMenu(scatterframe,datay,*DATA.options)
        ploty.pack(side=tk.TOP)
        
        #Checkbox that allows user to open the desired image in a new window.
        self.keepplot = tk.BooleanVar()
        self.keepplot.set(0)
        keepplotopt = tk.Checkbutton(scatterframe, text = "Open in new window",variable=self.keepplot) #Implement feature
        keepplotopt.pack(side=tk.TOP)
        
        #Button that allows user to plot both data sets in a scatter plat on the canvas.
        plotbutton = tk.Button(scatterframe, text='Plot', command=plot)
        plotbutton.pack(side=tk.TOP,pady=(0,5))
        
        #--Computation Frame Widgets--#
        
        #Option menu that allows user to choose data set to apply some sort of calculation.
        cmpchoice = tk.StringVar()
        cmpchoice.set('')
        cmpmenu = tk.OptionMenu(computeframe,cmpchoice,*DATA.options)
        cmpmenu.pack(side=tk.TOP)
        
        #These checkbuttons will most likely change to buttons, so temporary.
        
        #Frame to organize checkbuttons for Average and Difference computations.
        avgdif = tk.Frame(computeframe)
        avgdif.pack(side=tk.TOP)
        
        #Checkbutton that allows the user to generate an average between two slices of a 3-D data set.
        self.imgavg = tk.BooleanVar()
        imgavgopt = tk.Checkbutton(avgdif,text='Average',variable=self.imgavg)
        #imgavgopt = tk.Button(avgdif, text = 'Average', command=imgaverage)
        imgavgopt.grid(row=0,column=0)
        
        #Checkbutton that allows the user to generate a difference between two slices of a 3-D data set.
        self.imgdif = tk.BooleanVar()#
        imgdifopt = tk.Checkbutton(avgdif,text='Difference',variable=self.imgdif)
        #imgdifopt = tk.Button(avgdif, text = 'Difference', command=imgdifference)
        imgdifopt.grid(row=0,column=1)
        
        #Checkbutton that allows the user to generate a power spectra of a data set.
        self.ps = tk.BooleanVar()#
        powerspectraopt = tk.Checkbutton(computeframe, text = 'Power Spectra',variable=self.ps)
        #powerspectraopt = tk.Button(computeframe, text='Power Spectra', command=powerspectra)
        powerspectraopt.pack(side=tk.TOP)
        
        #Button that allows the user to compute the selected computations.
        computebutton = tk.Button(computeframe, text='Compute', command=compute)
        computebutton.pack(side=tk.TOP,pady=(0,5))
        
        #Button to open tutorial documentation
        documentation = tk.Button(sideframe, text='Documentation',command=None)
        documentation.pack(side=tk.TOP,pady=(10,0))
        
        #plot1.bind('<ButtonRelease-1>', tutorialSelect1)
        
        #Main figure in order to display generated plots.
        self.f=plt.figure(1)
        
        #Canvas widget which contains the figure to show plots.
        canvas = FigureCanvasTkAgg(self.f, frame)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
        
        #Displays Matplotlib figure toolbar.
        toolbar = NavigationToolbar2TkAgg(canvas, frame)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
app = Helioseismology()
#app.wm_attributes('-zoomed', True)
app.mainloop()
