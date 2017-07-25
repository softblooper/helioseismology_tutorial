#02 - A comparison of single image frames
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import matplotlib.image as mpimg

#Create Fits Image Object
class FitsImage(object):
    
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

class ScatterPlot(object):
    
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

Data1 = FitsImage('data1.fits','Data 1','gray','Velocity (m/s)')
Data2 = FitsImage('data2.fits','Data 2','gray','Velocity (m/s)')
D1vD2 = ScatterPlot(Data1.imgdata[0],Data2.imgdata[0],'Data 1 vs Data 2','Data 1 (m/s)','Data 2 (m/s)')

Data1.plot(2,2,0,0,1,1)
Data2.plot(2,2,1,0,1,1)
D1vD2.plot(2,2,0,1,2,1,-600,600,-600,600)
D1vD2.line(-600,600,-600,600)

mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()