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

#------------------------------------------------------------------------------#
#--Tutorial-Required Classes

class FitsImage(object): #Allows creating of plots with data from fits files
    
    #Initializes object, converts fits data into Python array. Detects shape for corresponding plot.
    def __init__(self,file,title,colorscale,colorlabel,shortname):
        self.file = file
        self.title = title
        self.color = colorscale
        self.colorlabel = colorlabel
        self.shortname = shortname
        
        self.imgdata = fits.getdata(self.file)
        
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

class Average(object): #Finds the average between two data sets and plots it.
    
    def __init__(self,data,title,colorscale,colorlabel):
        self.data = data
        self.title = title
        self.color = colorscale
        self.colorlabel = colorlabel
        
    def plot(self,rowdim,coldim,row,col,rowspan,colspan):
        plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
        plt.title(self.title)
        plt.imshow((self.data[0] + self.data[2])/2, cmap = self.color)
        plt.title(self.title)
        ibar1 = plt.colorbar()
        plt.gca().invert_yaxis()
        plt.xticks([])
        plt.yticks([])
        plt.grid(False)

class Difference(object): #Finds the difference between two data sets and plots it.
    
    def __init__(self,data,title,colorscale,colorlabel):
        self.data = data
        self.title = title
        self.color = colorscale
        self.colorlabel = colorlabel
        
    def plot(self,rowdim,coldim,row,col,rowspan,colspan):
        plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
        plt.title(self.title)
        plt.imshow((self.data[0] - self.data[2]), cmap = self.color)
        plt.title(self.title)
        ibar1 = plt.colorbar()
        plt.gca().invert_yaxis()
        plt.xticks([])
        plt.yticks([])
        plt.grid(False)

#------------------------------------------------------------------------------#
#--Tutorials--
Intensity = FitsImage('fd_Ic_6h_01d.fits','Intensity','gray','Continuum Intensity','intensity')
Magnetogram = FitsImage('fd_M_96m_01.fits','Magnetogram','gray','Guass (G)','magnetogram')
Dopplergram = FitsImage('fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)','dopplergram')
Data1 = FitsImage('data1.fits','Data 1','gray','Velocity (m/s)','data1')
Data2 = FitsImage('data2.fits','Data 2','gray','Velocity (m/s)','data2')
'''D1vD2 = ScatterPlot(Data1.imgdata[0],Data2.imgdata[0],'Data 1 vs Data 2','Data 1 (m/s)','Data 2 (m/s)')
Average1 = Average(Data1.imgdata,'Average of Data 1','gray','Velocity (m/s)')
Difference1 = Difference(Data1.imgdata,'Difference of Data 1','gray','Velocity (m/s)')
Average2 = Average(Data2.imgdata,'Average of Data 2','gray','Velocity (m/s)')
Difference2 = Difference(Data2.imgdata,'Difference of Data 2','gray','Velocity (m/s)')'''

Data = (Intensity, Magnetogram, Dopplergram, Data1, Data2) #D1vD2, Average1, Difference1,Average2,Difference2)
Options = ('Intensity', 'Magnetogram', 'Dopplergram', 'Data 1', 'Data 2')
DataIndex = (0, 1, 2, 3, 4)

DataOpts = dict(zip(Options,DataIndex))

gridoption = 1

def setgrid(number):
    global gridoption
    gridoption = number
    print(gridoption)
    return gridoption
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
        
        def plot():
            keep = self.keepplot.get()
            if keep:
                g = plt.figure(2)
                g.clear()
            else:
                self.f.clear()
                plt.figure(1)
            index1 = DataOpts[datax.get()]
            index2 = DataOpts[datay.get()]
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
            
        def viewimg():
            keep = self.keepimg.get()
            if keep:
                g = plt.figure(2)
                g.clear()
            else:
                self.f.clear()
                plt.figure(1)
            index = DataOpts[imgchoice.get()]
            minx = minrangex.get()
            maxx = maxrangex.get()
            miny = minrangey.get()
            maxy = maxrangey.get()
            minz = minrangez.get()
            maxz = maxrangez.get()
            if Data[index].dim == 3:
                if (minz and maxz):
                    minz = int(minz)
                    maxz = int(maxz)
                    plt.imshow(Data[index].imgdata[0], cmap = Data[index].color, vmin=minz, vmax=maxz)
                else:
                    plt.imshow(Data[index].imgdata[0], cmap = Data[index].color)
            else:
                if (minz and maxz):
                    minz = int(minz)
                    maxz = int(maxz)
                    plt.imshow(Data[index].imgdata, cmap = Data[index].color, vmin=minz, vmax=maxz)
                else:
                    plt.imshow(Data[index].imgdata, cmap = Data[index].color)
            if (minx and maxx):
                minx = int(minx)
                maxx = int(maxx)
                plt.xlim(minx,maxx)
            if (miny and maxy):
                miny = int(miny)
                maxy = int(maxy)
                plt.ylim(maxy,miny)
            ibar = plt.colorbar()
            ibar.set_label(Data[index].colorlabel)
            plt.gca().invert_yaxis()
            #plt.xticks([])
            #plt.yticks([])
            plt.xlabel('pix')
            plt.ylabel('pix')
            #plt.grid(False)
            if keep:
                kept = tk.Toplevel()
                
                keptcanvas = FigureCanvasTkAgg(g, kept)
                keptcanvas.show()
                keptcanvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
                
                kepttoolbar = NavigationToolbar2TkAgg(keptcanvas, kept)
                kepttoolbar.update()
                keptcanvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            else:
                self.f.canvas.draw()
            self.statusbar.config(text="Showing image of "+Data[index].title)
            minrangex.delete(0,tk.END)
            maxrangex.delete(0,tk.END)
            minrangey.delete(0,tk.END)
            maxrangey.delete(0,tk.END)
            minrangez.delete(0,tk.END)
            maxrangez.delete(0,tk.END)
        
        def imgaverage():
            index = DataOpts[cmpchoice.get()]
            average = (self.data[0] + self.data[2])/2
            hdu = fits.PrimaryHDU(average)
            name = 'avg_' + self.shortname
            hdu.writeto(name)
            
        def compute():
            avg = self.imgavg.get()
            if avg:
                imgaverage()
        
        container=tk.Frame(self)
        container.pack(side='top',fill='both',expand=True)
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0,weight=1)
        
        menu = tk.Frame(container)
        menu.pack(side=tk.TOP, fill=tk.X)
        
        self.statusbar = tk.Label(menu,text='Welcome to the Helioseismology Tutorial. Select data and have fun plotting!',relief=tk.SUNKEN,bg='white',width=60,anchor='w')
        self.statusbar.pack(anchor='nw',side=tk.LEFT,pady=(3,3))
        
        frame = tk.Frame(container)
        frame.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
        
        sideframe = tk.Frame(frame)
        sideframe.pack(side=tk.LEFT,anchor='n')
        
        scatterframe = tk.LabelFrame(sideframe, text="Scatter Plots")
        scatterframe.pack(side=tk.TOP,fill=tk.X)
        
        xaxis = tk.Label(scatterframe,text='X-Axis')
        xaxis.pack(side=tk.TOP)
        datax = tk.StringVar()
        datax.set('')
        plotx = tk.OptionMenu(scatterframe,datax,*Options)
        plotx.pack(side=tk.TOP)
        
        yaxis = tk.Label(scatterframe,text='Y-Axis')
        yaxis.pack(side=tk.TOP)
        datay = tk.StringVar()
        datay.set('')
        ploty = tk.OptionMenu(scatterframe,datay,*Options)
        ploty.pack(side=tk.TOP)
        
        self.keepplot = tk.BooleanVar()
        self.keepplot.set(0)
        keepplotopt = tk.Checkbutton(scatterframe, text = "Open in new window",variable=self.keepplot) #Implement feature
        keepplotopt.pack(side=tk.TOP)
        
        plotbutton = tk.Button(scatterframe, text='Plot', command=plot)
        plotbutton.pack(side=tk.TOP,pady=(0,5))
        
        imgframe = tk.LabelFrame(sideframe, text="Image Plots")
        imgframe.pack(side=tk.TOP,fill=tk.X)
        
        imgtext = tk.Label(imgframe,text='Image')
        imgtext.pack(side=tk.TOP)
        
        imgchoice = tk.StringVar()
        imgchoice.set('')
        imgmenu = tk.OptionMenu(imgframe,imgchoice,*Options)
        imgmenu.pack(side=tk.TOP)
        
        self.imgani = tk.BooleanVar()#
        imganiopt = tk.Checkbutton(imgframe, text = 'Animate',variable=self.imgani,command=None)
        imganiopt.pack(side=tk.TOP)
        
        self.keepimg = tk.BooleanVar()#
        keepimgopt = tk.Checkbutton(imgframe, text = "Open in new window",variable=self.keepimg,command=None) #Implement feature
        keepimgopt.pack(side=tk.TOP)
        
        imgbutton = tk.Button(imgframe, text='View', command=viewimg)
        imgbutton.pack(side=tk.TOP,pady=(0,5))
        
        computeframe = tk.LabelFrame(sideframe, text='Generate')
        computeframe.pack(side=tk.TOP,fill=tk.X)
        
        cmpchoice = tk.StringVar()
        cmpchoice.set('')
        cmpmenu = tk.OptionMenu(computeframe,cmpchoice,*Options)
        cmpmenu.pack(side=tk.TOP)
        
        avgdif = tk.Frame(computeframe)
        avgdif.pack(side=tk.TOP)
        
        self.imgavg = tk.BooleanVar()#
        imgavgopt = tk.Checkbutton(avgdif,text='Average',variable=self.imgavg,command=None)
        imgavgopt.grid(row=0,column=0)
        
        self.imgdif = tk.BooleanVar()#
        imgdifopt = tk.Checkbutton(avgdif,text='Difference',variable=self.imgdif,command=None)
        imgdifopt.grid(row=0,column=1)
        
        self.powerspectra = tk.BooleanVar()#
        powerspectraopt = tk.Checkbutton(computeframe, text = 'Power Spectra',variable=self.powerspectra,command = None)
        powerspectraopt.pack(side=tk.TOP)
        
        computebutton = tk.Button(computeframe, text='Compute', command=None)
        computebutton.pack(side=tk.TOP,pady=(0,5))
        
        toolsframe = tk.LabelFrame(sideframe,text='Plot Tools')
        toolsframe.pack(side=tk.TOP,fill=tk.X)
        
        xtools = tk.Label(toolsframe,text='X')
        xtools.pack(side=tk.TOP)
        
        rangeframex = tk.Frame(toolsframe)
        rangeframex.pack(side=tk.TOP)
        
        minlabelx = tk.Label(rangeframex,text = 'Min')
        minlabelx.grid(row=0,column=0)
        maxlabelx = tk.Label(rangeframex, text = 'Max')
        maxlabelx.grid(row=0,column=1)
        
        minrangex = tk.Entry(rangeframex,width = 6)
        minrangex.grid(row=1,column=0,padx=5)
        
        maxrangex = tk.Entry(rangeframex,width = 6)
        maxrangex.grid(row=1,column=1,padx=5)
        
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
        
        minrangez = tk.Entry(rangeframez,width = 6) #Implement range feature
        minrangez.grid(row=2,column=0,padx=5)
        
        maxrangez = tk.Entry(rangeframez,width = 6)
        maxrangez.grid(row=2,column=1,padx=5,pady=(0,5))
        
        documentation = tk.Button(sideframe, text='Documentation',command=None)
        documentation.pack(side=tk.TOP,pady=(10,0))
        
        #plot1.bind('<ButtonRelease-1>', tutorialSelect1)
        
        self.f=plt.figure(1)
        
        canvas = FigureCanvasTkAgg(self.f, frame)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
        
        toolbar = NavigationToolbar2TkAgg(canvas, frame)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

app = Helioseismology()
#app.wm_attributes('-zoomed', True)
app.mainloop()
