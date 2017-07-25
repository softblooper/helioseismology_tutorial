#03 - Acoustic signals
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

class Average(object):
    
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

class Difference(object):
    
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

Data1 = FitsImage('data1.fits','Data 1','gray','Velocity (m/s)')
Data2 = FitsImage('data2.fits','Data 2','gray','Velocity (m/s)')

Average1 = Average(Data1.imgdata,'Average of Data 1','gray','Velocity (m/s)')
Difference1 = Difference(Data1.imgdata,'Difference of Data 1','gray','Velocity (m/s)')

Average2 = Average(Data2.imgdata,'Average of Data 2','gray','Velocity (m/s)')
Difference2 = Difference(Data2.imgdata,'Difference of Data 2','gray','Velocity (m/s)')


Average1.plot(2,2,0,0,1,1)
Difference1.plot(2,2,0,1,1,1)
Average2.plot(2,2,1,0,1,1)
Difference2.plot(2,2,1,1,1,1)

mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()