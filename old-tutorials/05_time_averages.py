import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *
from astropy.io import fits

image_data1 = fits.getdata('data1.fits')
image_data2 = fits.getdata('data2.fits')
plt.gcf().canvas.set_window_title('Showing data1 & data2')
subplot(2,1,1)
plt.imshow(np.mean(image_data1, 0), cmap = 'gray')
plt.xlabel('x')
plt.ylabel('y')
divider=make_axes_locatable(plt.gca())
cax=divider.append_axes("right", size="10%", pad=0.05)
ibar1 = plt.colorbar(cax=cax)
ibar1.set_label('m/s')
plt.gca().invert_yaxis()
plt.subplot(2,1,2)
plt.imshow(np.mean(image_data2, 0), cmap = 'gray')
plt.xlabel('x')
plt.ylabel('y')
divider=make_axes_locatable(plt.gca())
cax=divider.append_axes("right", size="10%", pad=0.05)
ibar1 = plt.colorbar(cax=cax)
ibar1.set_label('m/s')
plt.gca().invert_yaxis()
plt.show()
