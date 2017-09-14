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
Intensity = FitsImage('fd_Ic_6h_01d.fits','Intensity','gray','Continuum Intensity','intensity','(1024x1024)')
Magnetogram = FitsImage('fd_M_96m_01.fits','Magnetogram','gray','Guass (G)','magnetogram','(1024x1024)')
Dopplergram = FitsImage('fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)','dopplergram','(1024x1024)')
Data1 = FitsImage('data1.fits','Data 1','gray','Velocity (m/s)','data1','(128x128x512)')
Data2 = FitsImage('data2.fits','Data 2','gray','Velocity (m/s)','data2','(128x128x512)')

Data = (Intensity, Magnetogram, Dopplergram, Data1, Data2) #D1vD2, Average1, Difference1,Average2,Difference2)
Options = ('Intensity', 'Magnetogram', 'Dopplergram', 'Data 1', 'Data 2')
DataIndex = (0, 1, 2, 3, 4)

DataOpts = dict(zip(Options,DataIndex))

ani = None

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
        
        #Open FITS file as image function
        def viewimg():
            ani = None
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
            
            if (minz and maxz):
                minz = int(minz)
                maxz = int(maxz)
            else:
                minz = None
                maxz = None
            
            if Data[index].dim == 3:
                plt.imshow(Data[index].imgdata[0], cmap = Data[index].color, vmin=minz, vmax=maxz)
            else:
                    plt.imshow(Data[index].imgdata, cmap = Data[index].color, vmin=minz, vmax=maxz)
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
        
        def viewslice():
            ani = None
            keep = self.keepimg.get()
            animate = self.animate.get()
            if keep:
                g = plt.figure(2)
                g.clear()
            else:
                self.f.clear()
                plt.figure(1)
            index = DataOpts[slicechoice.get()]
            minx = minrangex.get()
            maxx = maxrangex.get()
            miny = minrangey.get()
            maxy = maxrangey.get()
            minz = minrangez.get()
            maxz = maxrangez.get()
            slicex = xslice.get()
            slicey = yslice.get()
            slicet = tslice.get()
            
            labels = ('Velocity (m/s)','Time (s)','X-Pix','Y-Pix')
            labelindx = None
            labelindy = None
            
            if (minz and maxz):
                minz = int(minz)
                maxz = int(maxz)
            else:
                minz = None
                maxz = None
            
            if Data[index].dim == 3:
                sizet = np.arange(Data[index].shape[0])
                sizey = np.arange(Data[index].shape[1])
                sizex = np.arange(Data[index].shape[2])
                if slicex:
                    slicex = int(slicex)
                    if slicey:
                        slicey = int(slicey)
                        labelindx = 1
                        slicetext = 'x = '+xslice.get()+' and y = '+yslice.get()
                        data = Data[index].imgdata[:,slicey,slicex]
                        plt.plot(sizet,data,lw=0.5,color='black')
                    elif slicet:
                        slicet = int(slicet)
                        labelindx = 3
                        slicetext = 't = '+tslice.get()+' and x = '+xslice.get()
                        data = Data[index].imgdata[slicet,:,slicex]
                        plt.plot(sizey,data,lw=0.5,color='black')
                    else:
                        labelindx = 1
                        labelindy = 3
                        slicetext = 'x = '+xslice.get()
                        data = ndimage.rotate(Data[index].imgdata[:,:,slicex], 270)
                        plt.imshow(data, cmap = Data[index].color, vmin=minz, vmax=maxz)
                elif slicey:
                    slicey = int(slicey)
                    if slicet:
                        slicet = int(slicet)
                        labelindx = 2
                        slicetext = 't = '+tslice.get()+' and y = '+yslice.get()
                        data = Data[index].imgdata[slicet,slicey,:]
                        plt.plot(sizex,data,lw=0.5,color='black')
                    else:
                        labelindx = 2
                        labelindy = 1
                        slicetext = 'y = '+yslice.get()
                        data = Data[index].imgdata[:,slicey,:]
                        plt.imshow(data, cmap = Data[index].color, vmin=minz, vmax=maxz)
                elif slicet:
                    slicet = int(slicet)
                    labelindx = 2
                    labelindy = 3
                    slicetext = 't = '+tslice.get()
                    data = Data[index].imgdata[slicet,...]
                    plt.imshow(data, cmap = Data[index].color, vmin=minz, vmax=maxz)
                else:
                    data = Data[index].imgdata[0]
                    plt.imshow(data, cmap = Data[index].color, vmin=minz, vmax=maxz)
                    labelindx = 2
                    labelindy = 3
                    slicetext = 't = 0'
                
                if (slicex and slicey) or (slicex and slicet) or (slicey and slicet):
                    plt.title(labels[0]+' vs '+labels[labelindx]+' at '+slicetext)
                    plt.ylabel(labels[0])
                    plt.xlabel(labels[labelindx])
                else:
                    plt.title('Slice of '+Data[index].title+' at '+slicetext)
                    plt.xlabel(labels[labelindx])
                    plt.ylabel(labels[labelindy])
                    ibar = plt.colorbar()
                    ibar.set_label(Data[index].colorlabel)
            
            else: #I think these slice backwards
                sizey = np.arange(Data[index].shape[0])
                sizex = np.arange(Data[index].shape[1])
                if slicex:
                    slicex=int(slicex)
                    plt.plot(sizey,Data[index].imgdata[:,slicex],lw=0.5,color='black')
                    labelindx = 3
                    slicetext = 'x ='+xslice.get()
                elif slicey:
                    slicey=int(slicey)
                    plt.plot(sizex,Data[index].imgdata[slicey,:],lw=0.5,color='black')
                    labelindx = 2
                    slicetext = 'y ='+yslice.get()
                else:
                    plt.imshow(Data[index].imgdata, cmap = Data[index].color, vmin=minz, vmax=maxz)
                    labelindx = 2
                    labelindy = 3
                if (slicex or slicey):
                    plt.title(Data[index].colorlabel+' at '+slicetext)
                    plt.ylabel(Data[index].colorlabel)
                    plt.xlabel(labels[labelindx])
                else:
                    plt.title(Data[index].title)
                    plt.xlabel(labels[labelindx])
                    plt.ylabel(labels[labelindy])
                    ibar = plt.colorbar()
                    ibar.set_label(Data[index].colorlabel)
            
            '''if animate:
                if Data[index].dim == 3:
                    ims = []
                    for i in range(Data[index].shape[0]):
                        im = plt.imshow(Data[index].imgdata[i], cmap = Data[index].color, vmin=minz, vmax=maxz)
                        ims.append([im])
                    plt.title(Data[index].title)
                    ani = animation.ArtistAnimation(self.f, ims, interval=100, blit=True, repeat_delay=1000)'''
        
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
                
                text="Showing: "+Data[index].title+'. Size: '+Data[index].dimensions
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
                self.statusbar.config(text="Showing: "+Data[index].title+'. Size: '+Data[index].dimensions)
            
            self.f.canvas.mpl_connect('key_press_event', process_key)
            
            minrangex.delete(0,tk.END)
            maxrangex.delete(0,tk.END)
            minrangey.delete(0,tk.END)
            maxrangey.delete(0,tk.END)
            minrangez.delete(0,tk.END)
            maxrangez.delete(0,tk.END)
            tslice.delete(0,tk.END)
            xslice.delete(0,tk.END)
            yslice.delete(0,tk.END)
        
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
        
        #This is the main frame. It contains and organizes all other frames in the window.
        container=tk.Frame(self)
        container.pack(side='top',fill='both',expand=True)
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0,weight=1)
        
        #This frame is the top one, holding the status bar that follows.
        menu = tk.Frame(container)
        menu.pack(side=tk.TOP, fill=tk.X)
        
        #Status bar at top right of window. Gives feedback to user.
        self.statusbar = tk.Label(menu,text='Welcome to the Helioseismology Tutorial. Select data and have fun plotting!',relief=tk.SUNKEN,bg='maroon',fg='white',font='Arial 14 bold',width=220,anchor='center')
        self.statusbar.pack(anchor='center',pady=(3,3))
        
        #This frame is the frame containing the canvas and all plotting tools.
        frame = tk.Frame(container)
        frame.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
        
        #This frame is the one used to organize plotting tools.
        sideframe = tk.Frame(frame)
        sideframe.pack(side=tk.RIGHT,anchor='n')
        
        #--Menu Frames--#
        #All the main frames in order, as for ease of re-arrangement.
        
        #Tools for changing display of image or plot, such as setting X Y Z limits.
        toolsframe = tk.LabelFrame(sideframe,text='Plot Tools',fg='white',bg='grey',font='Arial 10 bold')
        toolsframe.pack(side=tk.TOP,fill=tk.X)
        
        #Frame for slice plots
        sliceframe = tk.LabelFrame(sideframe, text='Slice Plots',fg='white',bg='grey',font='Arial 10 bold')
        sliceframe.pack(side=tk.TOP,fill=tk.X)
        
        #Frame for image viewing options.
        #imgframe = tk.LabelFrame(sideframe, text="Image Plots")
        #imgframe.pack(side=tk.TOP,fill=tk.X)
        
        #Frame for scatter plot options
        scatterframe = tk.LabelFrame(sideframe, text="Scatter Plots",fg='white',bg='grey',font='Arial 10 bold')
        scatterframe.pack(side=tk.TOP,fill=tk.X)
        
        #Frame for computation of various functions to FITS files.
        computeframe = tk.LabelFrame(sideframe, text='Generate',fg='white',bg='grey',font='Arial 10 bold')
        computeframe.pack(side=tk.TOP,fill=tk.X)
        
        #--Tools Frame Widgets--#
        
        #This label simply marks X as to identify its corresponding entry boxes
        xtools = tk.Label(toolsframe,text='X',fg='white',bg='grey',font='Arial 10 bold')
        xtools.pack(side=tk.TOP)
        
        #Small frame to neatly organize min and max entry boxes and labels
        rangeframex = tk.Frame(toolsframe)
        rangeframex.pack(side=tk.TOP)
        
        #Labels for 'Min' and 'Max'
        minlabelx = tk.Label(rangeframex,text = 'Min',fg='white',bg='grey',font='Arial 10 bold')
        minlabelx.grid(row=0,column=0)
        maxlabelx = tk.Label(rangeframex, text = 'Max',fg='white',bg='grey',font='Arial 10 bold')
        maxlabelx.grid(row=0,column=1)
        
        #Entry box for X Minimum
        minrangex = tk.Entry(rangeframex,width = 6,bg='white')
        minrangex.grid(row=1,column=0,padx=5)
        
        #Entry box for X Maximum
        maxrangex = tk.Entry(rangeframex,width = 6,bg='white')
        maxrangex.grid(row=1,column=1,padx=5)
        
        #Previous descriptions repeat accordingtly to Y and Z as well
        ytools = tk.Label(toolsframe,text='Y',fg='white',bg='grey',font='Arial 10 bold')
        ytools.pack(side=tk.TOP)
        
        rangeframey = tk.Frame(toolsframe)
        rangeframey.pack(side=tk.TOP)
        
        minlabely = tk.Label(rangeframey,text = 'Min',fg='white',bg='grey',font='Arial 10 bold')
        minlabely.grid(row=0,column=0)
        maxlabely = tk.Label(rangeframey, text = 'Max',fg='white',bg='grey',font='Arial 10 bold')
        maxlabely.grid(row=0,column=1)
        
        minrangey = tk.Entry(rangeframey,width = 6)
        minrangey.grid(row=1,column=0,padx=5)
        
        maxrangey = tk.Entry(rangeframey,width = 6)
        maxrangey.grid(row=1,column=1,padx=5)
        
        ztools = tk.Label(toolsframe,text='Z',fg='white',bg='grey',font='Arial 10 bold')
        ztools.pack(side=tk.TOP)
        
        rangeframez = tk.Frame(toolsframe)
        rangeframez.pack(side=tk.TOP)
        
        minlabelz = tk.Label(rangeframez,text = 'Min',fg='white',bg='grey',font='Arial 10 bold')
        minlabelz.grid(row=0,column=0)
        maxlabelz = tk.Label(rangeframez, text = 'Max',fg='white',bg='grey',font='Arial 10 bold')
        maxlabelz.grid(row=0,column=1)
        
        minrangez = tk.Entry(rangeframez,width = 6)
        minrangez.grid(row=2,column=0,padx=5)
        
        maxrangez = tk.Entry(rangeframez,width = 6)
        maxrangez.grid(row=2,column=1,padx=5,pady=(0,5))
        
        '''#--Image Frame Widgets--#
        
        #Label to indicate the purpose of the following option menu.
        imgtext = tk.Label(imgframe,text='Image')
        imgtext.pack(side=tk.TOP)
        
        #Option menu to choose which image the user wants to view.
        imgchoice = tk.StringVar()
        imgchoice.set('')
        imgmenu = tk.OptionMenu(imgframe,imgchoice,*Options)
        imgmenu.pack(side=tk.TOP)
        
        #Checkbox that gives user option to animate the image. (Tempoarary, will be removed).
        self.imgani = tk.BooleanVar()#
        imganiopt = tk.Checkbutton(imgframe, text = 'Animate',variable=self.imgani,command=None)
        imganiopt.pack(side=tk.TOP)
        
        #Button that allows user to view selected image in the canvas.
        imgbutton = tk.Button(imgframe, text='View', command=viewimg)
        imgbutton.pack(side=tk.TOP,pady=(0,5))'''
        
        #--Scatter Plot Frame--#
        
        #Choose data set to be used for the x-axis
        xaxis = tk.Label(scatterframe,text='X-Axis',fg='white',bg='grey',font='Arial 10 bold')
        xaxis.pack(side=tk.TOP)
        datax = tk.StringVar()
        datax.set('')
        plotx = tk.OptionMenu(scatterframe,datax,*Options)
        plotx.config(bg='white')
        plotx.pack(side=tk.TOP)
        
        #Choose data set to be used for the y-axis
        yaxis = tk.Label(scatterframe,text='Y-Axis',fg='white',bg='grey',font='Arial 10 bold')
        yaxis.pack(side=tk.TOP)
        datay = tk.StringVar()
        datay.set('')
        ploty = tk.OptionMenu(scatterframe,datay,*Options)
        ploty.config(bg='white')
        ploty.pack(side=tk.TOP)
        
        #Checkbox that allows user to open the desired image in a new window.
        self.keepplot = tk.BooleanVar()
        self.keepplot.set(0)
        keepplotopt = tk.Checkbutton(scatterframe, text = "Open in new window",fg='white',bg='grey',font='Arial 10 bold',variable=self.keepplot) #Implement feature
        keepplotopt.pack(side=tk.TOP)
        
        #Button that allows user to plot both data sets in a scatter plat on the canvas.
        plotbutton = tk.Button(scatterframe, text='Plot',fg='white',bg='grey',font='Arial 10 bold', command=plot)
        plotbutton.pack(side=tk.TOP,pady=(0,5))
        
        #--Computation Frame Widgets--#
        
        #Option menu that allows user to choose data set to apply some sort of calculation.
        cmpchoice = tk.StringVar()
        cmpchoice.set('')
        cmpmenu = tk.OptionMenu(computeframe,cmpchoice,*Options)
        cmpmenu.config(bg='white')
        cmpmenu.pack(side=tk.TOP)
        
        #These checkbuttons will most likely change to buttons, so temporary.
        
        #Frame to organize checkbuttons for Average and Difference computations.
        avgdif = tk.Frame(computeframe)
        avgdif.pack(side=tk.TOP)
        
        #Checkbutton that allows the user to generate an average between two slices of a 3-D data set.
        self.imgavg = tk.BooleanVar()#
        imgavgopt = tk.Checkbutton(avgdif,text='Average',fg='white',bg='grey',font='Arial 10 bold',variable=self.imgavg,command=None)
        imgavgopt.grid(row=0,column=0)
        
        #Checkbutton that allows the user to generate a difference between two slices of a 3-D data set.
        self.imgdif = tk.BooleanVar()#
        imgdifopt = tk.Checkbutton(avgdif,text='Difference',fg='white',bg='grey',font='Arial 10 bold',variable=self.imgdif,command=None)
        imgdifopt.grid(row=0,column=1)
        
        #Checkbutton that allows the user to generate a power spectra of a data set.
        self.powerspectra = tk.BooleanVar()#
        powerspectraopt = tk.Checkbutton(computeframe, text = 'Power Spectra',fg='white',bg='grey',font='Arial 10 bold',variable=self.powerspectra,command = None)
        powerspectraopt.pack(side=tk.TOP)
        
        #Button that allows the user to compute the selected computations.
        computebutton = tk.Button(computeframe, text='Compute',fg='white',bg='grey',font='Arial 10 bold', command=None)
        computebutton.pack(side=tk.TOP,pady=(0,5))
        
        #--Slice Frame Widgets--#
        
        slicechoice = tk.StringVar()
        slicechoice.set('')
        slicemenu = tk.OptionMenu(sliceframe,slicechoice,*Options)
        slicemenu.config(bg='white')
        slicemenu.pack(side=tk.TOP)
        
        sliceentries = tk.Frame(sliceframe)
        sliceentries.pack(side=tk.TOP)
        
        xslicelabel = tk.Label(sliceentries,text='X',fg='white',bg='grey',font='Arial 10 bold')
        xslicelabel.grid(row=0,column=0,padx=5)
        
        yslicelabel = tk.Label(sliceentries,text='Y',fg='white',bg='grey',font='Arial 10 bold')
        yslicelabel.grid(row=0,column=1,padx=5)
        
        tslicelabel = tk.Label(sliceentries,text='t',fg='white',bg='grey',font='Arial 10 bold')
        tslicelabel.grid(row=0,column=2,padx=5)
        
        xslice = tk.Entry(sliceentries,width = 3)
        xslice.grid(row=1,column=0,padx=5)
        
        yslice = tk.Entry(sliceentries,width = 3)
        yslice.grid(row=1,column=1,padx=5)
        
        tslice = tk.Entry(sliceentries,width = 3)
        tslice.grid(row=1,column=2,padx=5)
        
        #Checkbox that allows user to open the desired image in a new window.
        self.keepimg = tk.BooleanVar()
        keepimgopt = tk.Checkbutton(sliceframe, text = "Open in new window",fg='white',bg='grey',font='Arial 10 bold',variable=self.keepimg)
        keepimgopt.pack(side=tk.TOP)
        
        self.animate = tk.BooleanVar()
        animateoption = tk.Checkbutton(sliceframe, text='Animate',fg='white',bg='grey',font='Arial 10 bold',variable=self.animate)
        animateoption.pack(side=tk.TOP,pady=(0,5))
        
        slicebutton = tk.Button(sliceframe, text='Plot',fg='white',bg='grey',font='Arial 10 bold',command=viewslice)
        slicebutton.pack(side=tk.TOP,pady=(0,5))
        
        #Button to open tutorial documentation
        documentation = tk.Button(sideframe, text='Documentation',fg='white',bg='grey',font='Arial 10 bold',command=None)
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
