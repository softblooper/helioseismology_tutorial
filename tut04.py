import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
import numpy as np
import matplotlib.animation as animation

#Function for plotting slices. Allows uneven distribution of figure.  
def plotfitsslices(row,col,title,imgdata,color,cbarlabel):
    ax = plt.subplot2grid((2, 3),(row, col),rowspan=2)
    plt.title(title)
    plt.imshow(imgdata, cmap=color)
    ibar= plt.colorbar()
    ibar.set_label(cbarlabel)
    plt.gca().invert_yaxis()
    plt.grid(False)
    
#Make function work?
'''def animate(title,imgdata,color):
    ims = []
    for i in range(512):
        im = plt.imshow(imgdata[i], cmap = color)
        plt.gcf().canvas.set_window_title(title)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.gca().invert_yaxis()
        ims.append([im])
    #ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True,
        #repeat_delay=1000)
        

def plotfitsani(row,col,title,imgdata,color,cbarlabel):
    ax = plt.subplot2grid((2, 3),(row, col))
    plt.title(title)
    plt.imshow(imgdata, cmap=color)
    ibar= plt.colorbar()
    ibar.set_label(cbarlabel)
    plt.gca().invert_yaxis()'''


#Data Info
fitsnplot=np.array([['data1.fits','Data 1 Time Slices','gray','Velocity (m/s)'],
    ['data2.fits','Data 2 Time Slices','gray','Velocity (m/s)']])

#Data and Text
fitsimage4=fits.open(fitsnplot[0,0])
imgdata4=fitsimage4[0].data
fitsimage5=fits.open(fitsnplot[1,0])
imgdata5=fitsimage4[0].data

#Figure information
fig=plt.figure(1, facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('4: Time Slices')

#Plotting Slices
plotfitsslices(0,0,fitsnplot[0,1],imgdata4[:,:,0],fitsnplot[0,2],fitsnplot[0,3])
plotfitsslices(0,2,fitsnplot[1,1],imgdata4[:,:,0],fitsnplot[1,2],fitsnplot[1,3])

#Plotting animations
ax3 = plt.subplot2grid((2, 3),(0, 1))
plt.title('Data 1')
ims1 = []
for i in range(512):
    im1 = plt.imshow(imgdata4[i], cmap = 'gray')
    ims1.append([im1])
plt.gca().invert_yaxis()
ani1 = animation.ArtistAnimation(fig, ims1, interval=100, blit=True,
    repeat_delay=1000)
  
ax4 = plt.subplot2grid((2, 3),(1, 1))
plt.title('Data 2')   
ims2 = []
for i in range(512):
    im2 = plt.imshow(imgdata5[i], cmap = 'gray')
    ims2.append([im2])
plt.gca().invert_yaxis()
ani2 = animation.ArtistAnimation(fig, ims2, interval=100, blit=True,
    repeat_delay=1000)
    
#Made autoanimation. How to control playback?

#Open Figure Maximized
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()