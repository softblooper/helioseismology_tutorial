import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
import textwrap as tw

#Opens FITS files
image_data1 = fits.getdata('data1.fits')
image_data2 = fits.getdata('data2.fits')
plt.gcf().canvas.set_window_title('Showing data1 & data2')

#Main function, creates the subplots and provides the style for the graphs.
subplot(2,2,1)
plt.imshow((image_data1[0] + image_data1[2])/2, cmap = 'gray', vmin = -500.00e00, vmax = 800e00)
plt.title("Average of frame 0 & 2")
plt.xlabel('x')
plt.ylabel('y')
ibar1 = plt.colorbar()
plt.gca().invert_yaxis()
plt.subplot(2,2,2)
plt.imshow(image_data1[0]-image_data1[2],cmap = 'gray', vmin = -300.00e00, vmax = 500e00)
plt.title("Difference of frame 0 & 2")
ibar2 = plt.colorbar()
ibar2.set_label('m/s')
plt.gca().invert_yaxis()

#Add the text. Font, style, location, width of the text box and more is created here.
txt_1 = '''
The images in the datacubes here are taken at a rate of one per minute.
A quick way to "see" the oscillations is to subtract two frames taken
about two minutes apart.'''
fig_txt = tw.fill(tw.dedent(txt_1.rstrip()), width=50)
plt.figtext(0.1, 0.3, fig_txt, horizontalalignment='left',
            fontsize=16, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1))
txt_2 = '''
Here are displayed both the average and the difference of two frames taken two
minutes apart for data1. After running it as is and examining the results,
modify this or copy to another batch file to display the same for data2.'''
fig_txt = tw.fill(tw.dedent(txt_2.rstrip()), width=50)
plt.figtext(0.1, 0.1, fig_txt, horizontalalignment='left',
            fontsize=16, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1))
txt_3 = '''
Why does the Evershed outflow disappear in the difference image?'''
fig_txt = tw.fill(tw.dedent(txt_3.rstrip()), width=50)
plt.figtext(0.55, 0.35, fig_txt, horizontalalignment='left',
            fontsize=16, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1), fontweight='bold')
txt_4 ='''What is appearing in the difference image?'''
fig_txt = tw.fill(tw.dedent(txt_4.rstrip()), width=50)
plt.figtext(0.55, 0.24, fig_txt, horizontalalignment='left',
            fontsize=16, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1), fontweight='bold')
txt_5 ='''How do the average and difference images compare between datasets 1 and 2?'''
fig_txt = tw.fill(tw.dedent(txt_5.rstrip()), width=50)
plt.figtext(0.55, 0.1, fig_txt, horizontalalignment='left',
            fontsize=16, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='white',
                      ec="0.5", pad=0.5, alpha=1), fontweight='bold')

plt.show()
