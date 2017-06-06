#01-Reading In The Data
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import matplotlib.image as mpimg
#from plotfits import plotfits???

#Functions
def plotfits(subsize,title,imgdata,color,cbarlabel):
    plt.subplot(subsize)
    plt.title(title)
    plt.imshow(imgdata, cmap=color)
    ibar= plt.colorbar()
    ibar.set_label(cbarlabel)
    plt.gca().invert_yaxis()

#Data Info
fitsnplot=np.array([['fd_Ic_6h_01d.fits','Intensity','gray','Continuum Intensity'],['fd_M_96m_01.fits','Magnetogram','gray','Gauss'],['fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)']])

#Data and Text
fitsimage1=fits.open(fitsnplot[0,0])
imgdata1=fitsimage1[0].data
fitsimage2=fits.open(fitsnplot[1,0])
imgdata2=fitsimage2[0].data
fitsimage3=fits.open(fitsnplot[2,0])
imgdata3=fitsimage3[0].data

tuttext1='''Here we have two datacubes, as well as read in and display \nsample images from the MDI instrument taken at during the \nsame interval (2002 March 31 15:00 - 23:32 UT). \nThe data cubes wll be the focus of the remainder of the exrcises. \nHowever, for now, examine the sample Dopplergram, Magnetogram \nand Continuum-Intensity image from the MDI instrument. \nIt is probable that you have seen simlar looking images as the \nintensity image, and magnetogram. The Dopplergram shows \nthe line-of-sight velocity of the photosphere at each pixel. \nA positive velocity indicates motion away from the observer (redshift). \nNotice that the other side of the Dopplergram is blue. Also notice \nthe patterns of red and blue structures which increase in contrast \ntowards the limb. This is supergranulation.'''
tuttext1q1='''If the units of the Dopplegram are in m/s, and the Sun has a radius of 696 Mm, can you determine th solar rotaion period?'''
tuttext1q2='''Why do you think supergranulation is more visible near the limb in these Dopplergrams?'''
tuttext1q3='''Can you see evidence of supergranulation in the continuum image or magnetogram? Try adjusting the minimum and maximum values of the range keyword in the call to tvim to increase the contrast'''

#Figure Format
plt.figure(1, facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('1: Reading in the data')

#Plotting
plotfits(221,'Intensity',imgdata1,'gray','Continuum Intensity')
plotfits(222,'Magnetogram',imgdata2,'gray','Gauss')
plotfits(223,'Dopplergram',imgdata3,'RdBu_r','Velocity (m/s)')

#Text in figure
plt.figtext(0.75, 0.25,tuttext1, wrap=True, ha='center', fontsize=12, bbox=dict(facecolor='white',edgecolor='black'))
plt.figtext(0.75, 0.20,tuttext1q1, wrap=True, ha='center', fontsize=12, bbox=dict(facecolor='white',edgecolor='black'))
plt.figtext(0.75, 0.17,tuttext1q2, wrap=True, ha='center', fontsize=12, bbox=dict(facecolor='white',edgecolor='black'))
plt.figtext(0.75, 0.102,tuttext1q3, wrap=True, ha='center', fontsize=12, bbox=dict(facecolor='white',edgecolor='black'))

#Open Figure Maximized
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()