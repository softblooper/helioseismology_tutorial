#04 - Time slices
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import matplotlib.image as mpimg
import matplotlib.animation as animation

#Create Fits Image Object
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
    def animateplot(self,numofslices,rowdim,coldim,row,col,rowspan,colspan):
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

#Figure information
f=plt.figure(1, facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('4: Time Slices')

Data1 = FitsImage('data1.fits','Data 1','gray','Velocity (m/s)')
Data2 = FitsImage('data2.fits','Data 2','gray','Velocity (m/s)')

Data1.sliceplot(2,3,0,0,2,1)
Data2.sliceplot(2,3,0,2,2,1)

ani1 = Data1.animateplot(512,2,3,0,1,1,1)
ani2 = Data2.animateplot(512,2,3,1,1,1,1)

Ani1 = animation.ArtistAnimation(f, ani1, interval=100, blit=True,
    repeat_delay=1000)
Ani2 = animation.ArtistAnimation(f, ani2, interval=100, blit=True,
    repeat_delay=1000)

#Open Figure Maximized
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()