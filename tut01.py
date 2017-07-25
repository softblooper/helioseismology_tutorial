#01-Reading in the data (general properties of the datacubes)
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

#Matplotlib Figure
plt.figure(1, facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('1: Reading in the data')

#Define images
Intensity = FitsImage('fd_Ic_6h_01d.fits','Intensity','gray','Continuum Intensity')
Magnetogram = FitsImage('fd_M_96m_01.fits','Magnetogram','gray','Guass')
Dopplergram = FitsImage('fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)')

#Plot images
Intensity.plot(2,2,0,0,1,1)
Magnetogram.plot(2,2,0,1,1,1)
Dopplergram.plot(2,2,1,0,1,1)

#Open Figure Maximized
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()