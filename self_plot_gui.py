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
    def __init__(self,file,title,colorscale,colorlabel):
        self.file = file
        self.title = title
        self.color = colorscale
        self.colorlabel = colorlabel
        
        self.imgdata = fits.getdata(self.file)
        
        self.shape = np.shape(self.imgdata)
        self.dim = len(self.shape)
    
    #Plotting method: Plots data and determines location and size of subplot. Grid option boolean.
    def plot(self,rowdim,coldim,row,col,rowspan,colspan):
        if self.dim == 3:
            plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
            plt.title(self.title)
            plt.imshow(self.imgdata[0], cmap = self.color)
            ibar = plt.colorbar()
            ibar.set_label(self.colorlabel)
            plt.gca().invert_yaxis()
            plt.xticks([])
            plt.yticks([])
            plt.grid(False)
        else:
            plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
            plt.title(self.title)
            plt.imshow(self.imgdata, cmap = self.color)
            ibar = plt.colorbar()
            ibar.set_label(self.colorlabel)
            plt.gca().invert_yaxis()
            plt.xticks([])
            plt.yticks([])
            plt.grid(False)
    
    #Plot Time Slices Method: Plots the layout of a slice of the datacube in order to show change over time.
    def sliceplot(self,rowdim,coldim,row,col,rowspan,colspan):
        plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
        plt.title(self.title+" Time Slices")
        plt.imshow(self.imgdata[:,64,:], cmap=self.color)
        ibar= plt.colorbar()
        ibar.set_label(self.colorlabel)
        plt.gca().invert_yaxis()
        plt.xticks([])
        plt.ylabel('Frame #')
        plt.grid(False)
    
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

class ScatterPlot(object): #Creates a scatterplot with given data
    
    #Initializes object, set paramaters for data used and other plot info
    def __init__(self,xdata,ydata,title,xlabel,ylabel):
        self.xdata = xdata
        self.ydata = ydata
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
    
    #Plotting method: Plots data and determines where to locate subplot
    def plot(self,rowdim,coldim,row,col,rowspan,colspan,xlimmin,xlimmax,ylimmin,ylimmax):
        plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
        plt.title(self.title)
        plt.plot(self.xdata, self.ydata,',', color='black')
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.xlim(xlimmin, xlimmax)
        plt.ylim(ylimmin, ylimmax)

    #Line method: Draws a line with set initial and end points.
    def line(self,xi,xf,yi,yf):
        x = np.linspace(xi,xf)
        y = np.linspace(yi,yf)
        plt.plot(x, y, color='black', linewidth = 1.0, linestyle = '-')
    
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
Intensity = FitsImage('fd_Ic_6h_01d.fits','Intensity','gray','Continuum Intensity')
Magnetogram = FitsImage('fd_M_96m_01.fits','Magnetogram','gray','Guass (G)')
Dopplergram = FitsImage('fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)')
Data1 = FitsImage('data1.fits','Data 1','gray','Velocity (m/s)')
Data2 = FitsImage('data2.fits','Data 2','gray','Velocity (m/s)')
D1vD2 = ScatterPlot(Data1.imgdata[0],Data2.imgdata[0],'Data 1 vs Data 2','Data 1 (m/s)','Data 2 (m/s)')
Average1 = Average(Data1.imgdata,'Average of Data 1','gray','Velocity (m/s)')
Difference1 = Difference(Data1.imgdata,'Difference of Data 1','gray','Velocity (m/s)')
Average2 = Average(Data2.imgdata,'Average of Data 2','gray','Velocity (m/s)')
Difference2 = Difference(Data2.imgdata,'Difference of Data 2','gray','Velocity (m/s)')

Data = (Intensity, Magnetogram, Dopplergram, Data1, Data2, D1vD2, Average1, Difference1,Average2,Difference2)
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

class Helioseismology(tk.Tk):
    
    def __init__(self,*args,**kwargs):
        
        tk.Tk.__init__(self,*args,**kwargs)
        
        tk.Tk.wm_title(self,'Helioseismology Tutorial')
        
        def plot():
            f.clear()
            #index1 = plot1.curselection()[0]
            index1 = DataOpts[data1x.get()]
            if gridoption == 1:
                plt.subplot2grid((1, 1), (0, 0))
                print(gridoption)
            elif gridoption == 2:
                plt.subplot2grid((1, 2), (0, 0))
                print(gridoption)
            elif gridoption == 3:
                plt.subplot2grid((1, 3), (0, 0))
                print(gridoption)
            elif gridoption == 4:
                plt.subplot2grid((2, 2), (0, 0))
                print(gridoption)
            plt.title(Data[index1].title)
            if self.imgcheck1:
                if Data[index1].dim == 3:
                    plt.imshow(Data[index1].imgdata[0], cmap = Data[index1].color)
                else:
                    plt.imshow(Data[index1].imgdata, cmap = Data[index1].color)
                ibar = plt.colorbar()
                ibar.set_label(Data[index1].colorlabel)
                plt.gca().invert_yaxis()
                plt.xticks([])
                plt.yticks([])
                plt.grid(False)
                f.canvas.draw()
                self.statusbar.config(text="Showing image of "+Data[index1].title)
            else:
                #index2 = plot2.curselection()[0]
                index2 = DataOpts[data1y.get()]
                try:
                    if (Data[index1].dim == 3) and (Data[index2].dim == 3):
                        plt.plot(Data[index1].imgdata[0],Data[index2].imgdata[0],',',color='black')
                    else:
                        plt.plot(Data[index1].imgdata,Data[index2].imgdata,',',color='black')
                    plt.xlabel(Data[index1].colorlabel)
                    plt.ylabel(Data[index2].colorlabel)
                    f.canvas.draw()
                    self.statusbar.config(text="Showing plot of "+Data[index1].title+" against "+Data[index2].title)
                except Exception:
                    self.statusbar.config(text="Can't plot these data sets! Differente sizes.")
        
        def imgcheckoption():
            self.imgcheck1 = not self.imgcheck1
        
        container=tk.Frame(self)
        container.pack(side='top',fill='both',expand=True)
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0,weight=1)
        
        self.menu = tk.Frame(container)
        self.menu.pack(side=tk.TOP, fill=tk.X)
        
        self.statusbar = tk.Label(self.menu,text='Welcome to the Helioseismology Tutorial. Select data and have fun plotting!',relief=tk.SUNKEN,bg='white',width=60,anchor='w')
        self.statusbar.pack(anchor='nw',side=tk.LEFT,pady=(3,3))
        
        self.multiplot4 = tk.Button(self.menu, text='2x2 Grid', command=setgrid(4))
        self.multiplot4.pack(side=tk.RIGHT)
        
        self.multiplot3 = tk.Button(self.menu, text='1x3 Grid', command=setgrid(3))
        self.multiplot3.pack(side=tk.RIGHT)
        
        self.multiplot2 = tk.Button(self.menu, text='1x2 Grid', command=setgrid(2))
        self.multiplot2.pack(side=tk.RIGHT)
        
        self.multiplot1 = tk.Button(self.menu, text='1x1 Grid', command=setgrid(1))
        self.multiplot1.pack(side=tk.RIGHT)
        
        self.frame = tk.Frame(container)
        self.frame.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
        
        self.sideframe = tk.Frame(self.frame)
        self.sideframe.pack(side=tk.LEFT,anchor='n')#,fill=tk.BOTH,expand=True)
        
        plotframe1 = tk.LabelFrame(self.sideframe, text="Plot 1 Tools")
        #plotframe1.pack(side=tk.TOP)
        plotframe1.grid(row=0)
        
        xaxis1 = tk.Label(plotframe1,text='Image/X-Axis')
        xaxis1.pack(side=tk.TOP)
        data1x = tk.StringVar()
        data1x.set('Intensity')
        plot1x = tk.OptionMenu(plotframe1,data1x,*Options)
        plot1x.pack(side=tk.TOP)
        
        self.imgcheck1 = False
        imgoption1 = tk.Checkbutton(plotframe1,text="Show FITS Image",command=imgcheckoption)
        imgoption1.pack(side=tk.TOP)
        
        yaxis1 = tk.Label(plotframe1,text='Y-Axis')
        yaxis1.pack(side=tk.TOP)
        data1y = tk.StringVar()
        data1y.set('Intensity')
        plot1y = tk.OptionMenu(plotframe1,data1y,*Options)
        plot1y.pack(side=tk.TOP)
        
        plotbutton1 = tk.Button(plotframe1, text='Plot', command=plot)
        plotbutton1.pack(side=tk.TOP)
        
        self.plotframe2 = tk.LabelFrame(self.sideframe, text="Plot 2 Tools")
        #self.plotframe2.pack(side=tk.TOP)
        self.plotframe2.grid(row=1)
        
        xaxis2 = tk.Label(self.plotframe2,text='Image/X-Axis')
        xaxis2.pack(side=tk.TOP)
        data2x = tk.StringVar()
        data2x.set('Intensity')
        plot2x = tk.OptionMenu(self.plotframe2,data2x,*Options)
        plot2x.pack(side=tk.TOP)
        
        self.imgcheck2 = False
        imgoption2 = tk.Checkbutton(self.plotframe2,text="Show FITS Image",command=imgcheckoption) #fix
        imgoption2.pack(side=tk.TOP)
        
        yaxis2 = tk.Label(self.plotframe2,text='Y-Axis')
        yaxis2.pack(side=tk.TOP)
        data2y = tk.StringVar()
        data2y.set('Intensity')
        plot2y = tk.OptionMenu(self.plotframe2,data2y,*Options)
        plot2y.pack(side=tk.TOP)
        
        plotbutton2 = tk.Button(self.plotframe2, text='Plot', command=plot)
        plotbutton2.pack(side=tk.TOP)
        
        plotframe3 = tk.LabelFrame(self.sideframe, text="Plot 3 Tools")
        #plotframe3.pack(side=tk.TOP)
        plotframe3.grid(row=2)
        
        xaxis3 = tk.Label(plotframe3,text='Image/X-Axis')
        xaxis3.pack(side=tk.TOP)
        data3x = tk.StringVar()
        data3x.set('Intensity')
        plot3x = tk.OptionMenu(plotframe3,data3x,*Options)
        plot3x.pack(side=tk.TOP)
        
        self.imgcheck3 = False
        imgoption3 = tk.Checkbutton(plotframe3,text="Show FITS Image",command=imgcheckoption) #fix
        imgoption3.pack(side=tk.TOP)
        
        yaxis3 = tk.Label(plotframe3,text='Y-Axis')
        yaxis3.pack(side=tk.TOP)
        data3y = tk.StringVar()
        data3y.set('Intensity')
        plot3y = tk.OptionMenu(plotframe3,data1y,*Options)
        plot3y.pack(side=tk.TOP)
        
        plotbutton3 = tk.Button(plotframe3, text='Plot', command=plot)
        plotbutton3.pack(side=tk.TOP)
        
        plotframe4 = tk.LabelFrame(self.sideframe, text="Plot 4 Tools")
        #plotframe4.pack(side=tk.TOP)
        plotframe4.grid(row=3)
        
        xaxis4 = tk.Label(plotframe4,text='Image/X-Axis')
        xaxis4.pack(side=tk.TOP)
        data4x = tk.StringVar()
        data4x.set('Intensity')
        plot4x = tk.OptionMenu(plotframe4,data4x,*Options)
        plot4x.pack(side=tk.TOP)
        
        self.imgcheck4 = False
        imgoption4 = tk.Checkbutton(plotframe4,text="Show FITS Image",command=imgcheckoption) #fix
        imgoption4.pack(side=tk.TOP)
        
        yaxis4 = tk.Label(plotframe4,text='Y-Axis')
        yaxis4.pack(side=tk.TOP)
        data4y = tk.StringVar()
        data4y.set('Intensity')
        plot4y = tk.OptionMenu(plotframe4,data1y,*Options)
        plot4y.pack(side=tk.TOP)
        
        plotbutton4 = tk.Button(plotframe4, text='Plot', command=plot)
        plotbutton4.pack(side=tk.TOP)
        
        #plot1.bind('<ButtonRelease-1>', tutorialSelect1)
        
        f=plt.figure()
        
        canvas = FigureCanvasTkAgg(f, self.frame)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
        
        toolbar = NavigationToolbar2TkAgg(canvas, self.frame)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

app = Helioseismology()
app.wm_attributes('-zoomed', True)
app.mainloop()