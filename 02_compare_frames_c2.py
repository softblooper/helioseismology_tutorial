import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *
from astropy.io import fits
import textwrap as tw

image_data1 = fits.getdata('data1.fits')
image_data2 = fits.getdata('data2.fits')
plt.gcf().canvas.set_window_title('Showing data1 & data2')
subplot(2,2,1)
plt.imshow(image_data1[0],cmap = 'gray')
plt.gca().invert_yaxis()
divider=make_axes_locatable(plt.gca())
cax=divider.append_axes("right", size="10%", pad=0.05)
plt.colorbar(cax=cax)
plt.subplot(2,2,3)
plt.imshow(image_data2[0],cmap = 'gray')
plt.gca().invert_yaxis()
divider=make_axes_locatable(plt.gca())
cax=divider.append_axes("right", size="10%", pad=0.05)
plt.colorbar(cax=cax)
txt_1 = '''
The two datasets (data1 and data2) are three-dimensional arrays
(read in from the original FITS format).
One of these simultaneous datacubes is from the MDI instrument
and the other is from the ground-based GONG network.'''
fig_txt = tw.fill(tw.dedent(txt_1.rstrip()), width=50)
plt.figtext(0.56, 0.7, fig_txt, horizontalalignment='left',
            fontsize=16, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1))
txt_2 = '''
Can you guess which is GONG and which is MDI? How? If you cannot,
more clues will accumulate along the way.'''
fig_txt = tw.fill(tw.dedent(txt_2.rstrip()), width=50)
plt.figtext(0.56, 0.2, fig_txt, horizontalalignment='left',
            fontsize=16, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1), fontweight='bold')

plt.figure(6, figsize = (8, 8), facecolor = 'white', edgecolor = 'k')
plt.gcf().canvas.set_window_title('data1 vs. data2')
X = np.linspace(-800, 1000, 1800, endpoint = True)
Z = np.linspace(0, 0, 1800, endpoint = True)
Y = X
subplot(1,2,1)
plt.plot(image_data1[0], image_data2[0],'ro', color='gray')
plt.plot(X, Y, color='black', linewidth = 1.5, linestyle = '-')
plt.plot(X, Z, color='black', linewidth = 1.0, linestyle = '--')
plt.plot(Z, Y, color='black', linewidth = 1.0, linestyle = '--')
plt.xlabel('data1 frame 0')
plt.ylabel('data2 frame 0')
plt.xlim(-600, 800)
plt.ylim(-600, 800)
txt_3 = '''
Both of the datasets are calibrated in units of m/s.
However the scatter plot does not have a slope of unity,
even though the Dopplergrams were taken at the same time.'''
fig_txt = tw.fill(tw.dedent(txt_3.rstrip()), width=50)
plt.figtext(0.53, 0.7, fig_txt, horizontalalignment='left',
            fontsize=14, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1))
txt_4 = '''
Can you understand why the scatter plot looks the way it does?'''
fig_txt = tw.fill(tw.dedent(txt_4.rstrip()), width=50)
plt.figtext(0.53, 0.58, fig_txt, horizontalalignment='left',
            fontsize=14, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1), fontweight='bold')
txt_5 = '''
The spatial area covered by the datacubes represent only a
small 256 by 256 pixel subarray of the full disk images you
displayed in exercise 1. The pixels here are each .002 solar
radii, or 1.39 Mm in size. The area is centered on a sunspot,
and the datacube is tracked along with a constant solar rotation
rate (known as the Carrington rate).'''
fig_txt = tw.fill(tw.dedent(txt_5.rstrip()), width=50)
plt.figtext(0.53, 0.25, fig_txt, horizontalalignment='left',
            fontsize=14, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1))
txt_6 = '''
Sunspots have outflows called Evershed flows.
Why does the Evershed outflow look the way it does here?'''
fig_txt = tw.fill(tw.dedent(txt_6.rstrip()), width=50)
plt.figtext(0.53, 0.1, fig_txt, horizontalalignment='left',
            fontsize=14, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1), fontweight='bold')
plt.show()
