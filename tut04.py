import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
import numpy as np
import matplotlib.animation as animation


def plotfits(subsize,title,imgdata,color,cbarlabel):
    plt.subplot(subsize)
    plt.title(title)
    plt.imshow(imgdata, cmap=color)
    ibar= plt.colorbar()
    ibar.set_label(cbarlabel)
    plt.gca().invert_yaxis()
    
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
        #repeat_delay=1000)'''


#Data Info
fitsnplot=np.array([['data1.fits','Data 1 Time Slices','gray','Velocity (m/s)'],['data2.fits','Data 2 Time Slices','gray','Velocity (m/s)']])

#Data and Text
fitsimage4=fits.open(fitsnplot[0,0])
imgdata4=fitsimage4[0].data
fitsimage5=fits.open(fitsnplot[1,0])
imgdata5=fitsimage4[0].data

plt.figure(1, facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('4: Time Slices')
plotfits(121,fitsnplot[0,1],imgdata4[:,:,0],fitsnplot[0,2],fitsnplot[0,3])
plotfits(122,fitsnplot[1,1],imgdata4[:,:,0],fitsnplot[1,2],fitsnplot[1,3])

plt.figure(2, facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('4: Time Slices')
#animate('Dat 1',imgdata4,'gray')
'''ims = []
for i in range(512):
    im = plt.imshow(imgdata4[i], cmap = 'gray')
    plt.gcf().canvas.set_window_title('data1')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().invert_yaxis()
    ims.append([im])'''
#Fix old animation. Maybe try to function it? Maybe it's wrong? Ask Dr. Jwicz

plt.show()