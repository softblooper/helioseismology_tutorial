import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
import numpy as np
import matplotlib.animation as animation

#Functions.
#Allows plotting of slices.
def plotfitsslices(row,col,title,imgdata,color,cbarlabel):
    ax = plt.subplot2grid((2, 3),(row, col),rowspan=2)
    plt.title(title)
    plt.imshow(imgdata, cmap=color)
    ibar= plt.colorbar()
    ibar.set_label(cbarlabel)
    plt.gca().invert_yaxis()
    plt.grid(False)
    
#Allows animations and plots them.
def animate(row,col,title,imgdata,color):
    ax = plt.subplot2grid((2, 3),(row, col))
    plt.title(title)
    plt.gca().invert_yaxis()
    ims = []
    for i in range(512):
        im = plt.imshow(imgdata[i], cmap = color)
        ims.append([im])
    ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True,
        repeat_delay=1000)
    return ani

#Data Info
fitsnplot=np.array([['data1.fits','Data 1 Time Slices','gray','Velocity (m/s)'],
    ['data2.fits','Data 2 Time Slices','gray','Velocity (m/s)']])

imgdata4=fits.getdata(fitsnplot[0,0])
imgdata5=fits.getdata(fitsnplot[1,0])

#Figure information
fig=plt.figure(1, facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('4: Time Slices')

#Plotting Slices
plotfitsslices(0,0,fitsnplot[0,1],imgdata4[:,:,0],fitsnplot[0,2],fitsnplot[0,3])
plotfitsslices(0,2,fitsnplot[1,1],imgdata5[:,:,0],fitsnplot[1,2],fitsnplot[1,3])

#Plotting animations
ani1 = animate(0,1,'Data 1',imgdata4,'gray')
ani2 = animate(1,1,'Data 2',imgdata5,'gray')
    
#Made autoanimation. How to control playback?

#Open Figure Maximized
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()
