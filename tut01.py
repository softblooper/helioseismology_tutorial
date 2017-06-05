import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import matplotlib.image as mpimg
#Data Info
fitsnplot=np.array([['fd_Ic_6h_01d.fits','Intensity','gray','Continuum Intensity'],['fd_M_96m_01.fits','Magnetogram','gray','Gauss'],['fd_V_01h.fits','Dopplergram','RdBu_r','Velocity (m/s)']])

#Data and Text
fitsimage1=fits.open(fitsnplot[0,0])
imgdata1=fitsimage1[0].data
fitsimage2=fits.open(fitsnplot[1,0])
imgdata2=fitsimage2[0].data
fitsimage3=fits.open(fitsnplot[2,0])
imgdata3=fitsimage3[0].data
tuttext1='''Here we have two datacubes, as well as read in and display \nsample images from the MDI instrument taken at during the \nsame interval (2002 March 31 15:00 - 23:32 UT). \nThe data cubes wll be the focus of the remainder of the exrcises. \nHowever, for now, examine the sample Dopplergram, Magnetogram \nand Continuum-Intensity image from the MDI instrument. \nIt is probable that you have seen simlar looking images as the \nintensity image, and magnetogram. The Dopplergram shows \nthe line-of-sight velocity of the photosphere at each pixel. \nA positive velocity indicates motion away from the observer (redshift). \nNotice that the other side of the Dopplergram is blue. Also notice the patterns of red and blue structures which increase in contrast towards the limb. This is supergranulation.'''
tuttext1q1='''If the units of the Dopplegram are in m/s, and the Sun has a radius of 696 Mm, can you determine th solar rotaion period?'''
tuttext1q2='''Why do you think supergranulation is more visible near the limb in these Dopplergrams?'''
tuttext1q3='''Can you see evidence of supergranulation in the continuum image or magnetogram? Try adjusting the minimum and maximum values of the range keyword in the call to tvim to increase the contrast'''

#Figure Format
plt.figure(1, facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('1: Reading in the data')

#Plotting
plt.subplot(2,2,1)
plt.title(fitsnplot[0,1])
plt.imshow(imgdata1, cmap=fitsnplot[0,2])
ibar1 = plt.colorbar()
ibar1.set_label(fitsnplot[0,3])
plt.gca().invert_yaxis()

plt.subplot(2,2,2)
plt.title(fitsnplot[1,1])
plt.imshow(imgdata2, cmap=fitsnplot[1,2])
ibar2 = plt.colorbar()
ibar2.set_label(fitsnplot[1,3])
plt.gca().invert_yaxis()

plt.subplot(2,2,3)
plt.title(fitsnplot[2,1])
plt.imshow(imgdata3, cmap=fitsnplot[2,2])
ibar3 = plt.colorbar()
ibar3.set_label(fitsnplot[2,3])
plt.gca().invert_yaxis()

#Text in figure
plt.figtext(0.75, 0.25,tuttext1, wrap=True, ha='center', fontsize=12, bbox=dict(facecolor='white',edgecolor='black'))
plt.figtext(0.75, 0.20,tuttext1q1, wrap=True, ha='center', fontsize=12, bbox=dict(facecolor='white',edgecolor='black'))
plt.figtext(0.75, 0.17,tuttext1q2, wrap=True, ha='center', fontsize=12, bbox=dict(facecolor='white',edgecolor='black'))
plt.figtext(0.75, 0.102,tuttext1q3, wrap=True, ha='center', fontsize=12, bbox=dict(facecolor='white',edgecolor='black'))

#Open Figure Maximized
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()