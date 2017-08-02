#--Helioseismology Tutorial--
#--Dr. Jason Jackiewicz, Kimberly Navarro, Jorge Garcia

#Import statements
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

#import os

#States
Main_Font=('Arial',12)

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
            plt.subplot2grid((rowdim, coldim), (row, col), rowspan = rowspan, colspan = colspan)
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
Magnetogram = FitsImage('fd_M_96m_01.fits','Magnetogram','gray','Guass')
Dopplergram = FitsImage('fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)')
Data1 = FitsImage('data1.fits','Data 1','gray','Velocity (m/s)')
Data2 = FitsImage('data2.fits','Data 2','gray','Velocity (m/s)')
D1vD2 = ScatterPlot(Data1.imgdata[0],Data2.imgdata[0],'Data 1 vs Data 2','Data 1 (m/s)','Data 2 (m/s)')
Average1 = Average(Data1.imgdata,'Average of Data 1','gray','Velocity (m/s)')
Difference1 = Difference(Data1.imgdata,'Difference of Data 1','gray','Velocity (m/s)')
Average2 = Average(Data2.imgdata,'Average of Data 2','gray','Velocity (m/s)')
Difference2 = Difference(Data2.imgdata,'Difference of Data 2','gray','Velocity (m/s)')

f1=plt.figure()
Intensity.plot(2,2,0,0,1,1)
Magnetogram.plot(2,2,0,1,1,1)
Dopplergram.plot(2,2,1,0,1,1)

f2=plt.figure()
Data1.plot(2,2,0,0,1,1)
Data2.plot(2,2,1,0,1,1)
D1vD2.plot(2,2,0,1,2,1,-600,600,-600,600)
D1vD2.line(-600,600,-600,600)

f3=plt.figure()
Average1.plot(2,2,0,0,1,1)
Difference1.plot(2,2,0,1,1,1)
Average2.plot(2,2,1,0,1,1)
Difference2.plot(2,2,1,1,1,1)

f4 = plt.figure()
Data1.sliceplot(2,3,0,0,2,1)    
Data2.sliceplot(2,3,0,2,2,1)
ani1 = Data1.animateplot(512,2,3,0,1,1,1)
ani2 = Data2.animateplot(512,2,3,1,1,1,1)

#------------------------------------------------------------------------------#
#--GUI--

class Helioseismology(tk.Tk):
    
    def __init__(self,*args,**kwargs):
        
        tk.Tk.__init__(self,*args,**kwargs)
        
        tk.Tk.wm_title(self,'Helioseismology Tutorial')
        
        container=tk.Frame(self)
        container.pack(side='top',fill='both',expand=True)
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0,weight=1)
        
        menubar=tk.Menu(container)
        filemenu=tk.Menu(menubar,tearoff=0)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=quit)
        menubar.add_cascade(label="Options",menu=filemenu)
        
        tk.Tk.config(self,menu=menubar)
        
        self.frames = {}
            
        for F in (StartPage,Tutorial1,Tutorial2,Tutorial3,Tutorial4,Documentation):
            
            frame = F(container,self)
            self.frames[F]=frame
            frame.grid(row=0,column=0,sticky='nsew')
        
        self.show_frame(StartPage)
        
        #tk.Tk.iconbitmap(self,default='exampleicon.ico') Should be icon. Not working.

    def show_frame(self, cont):
        frame=self.frames[cont]
        frame.tkraise()
        

class StartPage(tk.Frame):
    
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
                
        button1 = tk.Button(self, text='Home',relief=tk.SUNKEN,
            command=lambda: controller.show_frame(StartPage))
        button1.pack(padx=5,pady=5,anchor='nw',side=tk.TOP)
        
        tutframe = tk.LabelFrame(self, text="Tutorials")
        tutframe.pack(padx=5,pady=5,anchor='nw',side=tk.LEFT)
        
        tut1 = tk.Button(tutframe, text='Tutorial 1',
            command=lambda: controller.show_frame(Tutorial1))
        tut1.pack(padx=2,pady=5)
        
        tut2 = tk.Button(tutframe, text='Tutorial 2',
            command=lambda: controller.show_frame(Tutorial2))
        tut2.pack(pady=5)
        
        tut3 = tk.Button(tutframe, text='Tutorial 3',
            command=lambda: controller.show_frame(Tutorial3))
        tut3.pack(pady=5)
        
        tut4 = tk.Button(tutframe, text='Tutorial 4',
            command=lambda: controller.show_frame(Tutorial4))
        tut4.pack(pady=5)
        
        tut5 = tk.Button(tutframe, text='Tutorial 5',
            command=lambda: controller.show_frame(StartPage))
        tut5.pack(pady=5)
        
        tut6 = tk.Button(tutframe, text='Tutorial 6',
            command=lambda: controller.show_frame(StartPage))
        tut6.pack(pady=5)
        
        tut7 = tk.Button(tutframe, text='Tutorial 7',
            command=lambda: controller.show_frame(StartPage))
        tut7.pack(pady=5)
        
        docbutton=tk.Button(tutframe,text='Documentation',
            command=lambda: controller.show_frame(Documentation))
        docbutton.pack(pady=5)
        
        greeting = tk.Label(self, text='Welcome to the Helioseismology Tutorial!',font=Main_Font)
        greeting2 = tk.Label(self, text='Select a tutorial to start.',font=Main_Font)
        greeting3 = tk.Label(self, text='Maximize window for best display.',font=Main_Font)
        greeting.pack(pady=20)
        greeting2.pack(pady=20)
        greeting3.pack()


class Tutorial1(tk.Frame):
    
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        
        button1 = tk.Button(self, text='Home',
            command=lambda: controller.show_frame(StartPage))
        button1.pack(padx=5,pady=5,anchor='nw',side=tk.TOP)
        
        tutframe = tk.LabelFrame(self, text="Tutorials")
        tutframe.pack(padx=5,pady=5,anchor='nw',side=tk.LEFT)
        
        tut1 = tk.Button(tutframe, text='Tutorial 1', relief=tk.SUNKEN,
            command=lambda: controller.show_frame(Tutorial1))
        tut1.pack(padx=2,pady=5)
        
        tut2 = tk.Button(tutframe, text='Tutorial 2',
            command=lambda: controller.show_frame(Tutorial2))
        tut2.pack(pady=5)
        
        tut3 = tk.Button(tutframe, text='Tutorial 3',
            command=lambda: controller.show_frame(Tutorial3))
        tut3.pack(pady=5)
        
        tut4 = tk.Button(tutframe, text='Tutorial 4',
            command=lambda: controller.show_frame(Tutorial4))
        tut4.pack(pady=5)
        
        tut5 = tk.Button(tutframe, text='Tutorial 5',
            command=lambda: controller.show_frame(StartPage))
        tut5.pack(pady=5)
        
        tut6 = tk.Button(tutframe, text='Tutorial 6',
            command=lambda: controller.show_frame(StartPage))
        tut6.pack(pady=5)
        
        tut7 = tk.Button(tutframe, text='Tutorial 7',
            command=lambda: controller.show_frame(StartPage))
        tut7.pack(pady=5)
        
        docbutton=tk.Button(tutframe,text='Documentation',
            command=lambda: controller.show_frame(Documentation))
        docbutton.pack(pady=5)
        
        canvas = FigureCanvasTkAgg(f1, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
        
        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class Tutorial2(tk.Frame):
    
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        
        button1 = tk.Button(self, text='Home',
            command=lambda: controller.show_frame(StartPage))
        button1.pack(padx=5,pady=5,anchor='nw',side=tk.TOP)
        
        tutframe = tk.LabelFrame(self, text="Tutorials")
        tutframe.pack(padx=5,pady=5,anchor='nw',side=tk.LEFT)
        
        tut1 = tk.Button(tutframe, text='Tutorial 1',
            command=lambda: controller.show_frame(Tutorial1))
        tut1.pack(padx=2,pady=5)
        
        tut2 = tk.Button(tutframe, text='Tutorial 2', relief=tk.SUNKEN,
            command=lambda: controller.show_frame(Tutorial2))
        tut2.pack(pady=5)
        
        tut3 = tk.Button(tutframe, text='Tutorial 3',
            command=lambda: controller.show_frame(Tutorial3))
        tut3.pack(pady=5)
        
        tut4 = tk.Button(tutframe, text='Tutorial 4',
            command=lambda: controller.show_frame(Tutorial4))
        tut4.pack(pady=5)
        
        tut5 = tk.Button(tutframe, text='Tutorial 5',
            command=lambda: controller.show_frame(StartPage))
        tut5.pack(pady=5)
        
        tut6 = tk.Button(tutframe, text='Tutorial 6',
            command=lambda: controller.show_frame(StartPage))
        tut6.pack(pady=5)
        
        tut7 = tk.Button(tutframe, text='Tutorial 7',
            command=lambda: controller.show_frame(StartPage))
        tut7.pack(pady=5)
        
        docbutton=tk.Button(tutframe,text='Documentation',
            command=lambda: controller.show_frame(Documentation))
        docbutton.pack(pady=5)
        
        canvas = FigureCanvasTkAgg(f2, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
        
        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class Tutorial3(tk.Frame):
    
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)

        button1 = tk.Button(self, text='Home',
            command=lambda: controller.show_frame(StartPage))
        button1.pack(padx=5,pady=5,anchor='nw',side=tk.TOP)
        
        tutframe = tk.LabelFrame(self, text="Tutorials")
        tutframe.pack(padx=5,pady=5,anchor='nw',side=tk.LEFT)
        
        tut1 = tk.Button(tutframe, text='Tutorial 1',
            command=lambda: controller.show_frame(Tutorial1))
        tut1.pack(padx=2,pady=5)
        
        tut2 = tk.Button(tutframe, text='Tutorial 2',
            command=lambda: controller.show_frame(Tutorial2))
        tut2.pack(pady=5)
        
        tut3 = tk.Button(tutframe, text='Tutorial 3', relief=tk.SUNKEN,
            command=lambda: controller.show_frame(Tutorial3))
        tut3.pack(pady=5)
        
        tut4 = tk.Button(tutframe, text='Tutorial 4',
            command=lambda: controller.show_frame(Tutorial4))
        tut4.pack(pady=5)
        
        tut5 = tk.Button(tutframe, text='Tutorial 5',
            command=lambda: controller.show_frame(StartPage))
        tut5.pack(pady=5)
        
        tut6 = tk.Button(tutframe, text='Tutorial 6',
            command=lambda: controller.show_frame(StartPage))
        tut6.pack(pady=5)
        
        tut7 = tk.Button(tutframe, text='Tutorial 7',
            command=lambda: controller.show_frame(StartPage))
        tut7.pack(pady=5)
        
        docbutton=tk.Button(tutframe,text='Documentation',
            command=lambda: controller.show_frame(Documentation))
        docbutton.pack(pady=5)
        
        canvas = FigureCanvasTkAgg(f3, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
        
        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class Tutorial4(tk.Frame):
    
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)

        button1 = tk.Button(self, text='Home',
            command=lambda: controller.show_frame(StartPage))
        button1.pack(padx=5,pady=5,anchor='nw',side=tk.TOP)
        
        tutframe = tk.LabelFrame(self, text="Tutorials")
        tutframe.pack(padx=5,pady=5,anchor='nw',side=tk.LEFT)
        
        tut1 = tk.Button(tutframe, text='Tutorial 1',
            command=lambda: controller.show_frame(Tutorial1))
        tut1.pack(padx=2,pady=5)
        
        tut2 = tk.Button(tutframe, text='Tutorial 2',
            command=lambda: controller.show_frame(Tutorial2))
        tut2.pack(pady=5)
        
        tut3 = tk.Button(tutframe, text='Tutorial 3',
            command=lambda: controller.show_frame(Tutorial3))
        tut3.pack(pady=5)
        
        tut4 = tk.Button(tutframe, text='Tutorial 4', relief=tk.SUNKEN,
            command=lambda: controller.show_frame(Tutorial4))
        tut4.pack(pady=5)
        
        tut5 = tk.Button(tutframe, text='Tutorial 5',
            command=lambda: controller.show_frame(StartPage))
        tut5.pack(pady=5)
        
        tut6 = tk.Button(tutframe, text='Tutorial 6',
            command=lambda: controller.show_frame(StartPage))
        tut6.pack(pady=5)
        
        tut7 = tk.Button(tutframe, text='Tutorial 7',
            command=lambda: controller.show_frame(StartPage))
        tut7.pack(pady=5)
        
        docbutton=tk.Button(tutframe,text='Documentation',
            command=lambda: controller.show_frame(Documentation))
        docbutton.pack(pady=5)
        
        canvas = FigureCanvasTkAgg(f4, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=True)
        
        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class Documentation(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        
        button1 = tk.Button(self, text='Back',
            command=lambda: controller.show_frame(StartPage))
        button1.pack()

app = Helioseismology()
Ani1 = animation.ArtistAnimation(f4, ani1, interval=100, blit=True,
    repeat_delay=1000)
Ani2 = animation.ArtistAnimation(f4, ani2, interval=100, blit=True,
    repeat_delay=1000)
app.wm_attributes('-zoomed', True)
app.mainloop()