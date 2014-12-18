#!/usr/local/bin/python
#------------------------------------------------------------------------------------------------
import numpy as np
import glob
import os
import astropy.io.fits as pf
from circAperture import makeCircularMask
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
#------------------------------------------------------------------------------------------------
'''
Code to align two images in x, and y for saturated PSF. 

I use the aery rings to find the best matching x,y values and shift the images. 
To do this, I use an annular ring (size can be adjusted appropriately), and an
initial guess of the stellar center (within +- 5px).

To run use the following command: 
python pyAlign.py

v1 12/14/2014 - Abhijith Rajan
'''
#------------------------------------------------------------------------------------------------
def main():
    fitslst = glob.glob('*skysub.fits')
    fitslst.sort()

    ref = pf.getdata(fitslst[0])
    xref,yref = np.loadtxt(fitslst[0]+'.coo')
    refmask = makeCircularMask(ref.shape,15,xref,yref) - makeCircularMask(ref.shape,8,xref,yref) 
    refmed = np.median(ref[refmask])
    pf.writeto('aligned_'+fitslst[0], ref)

    data = [ref]
    for fits in fitslst[1:]:
        targ = pf.getdata(fits)
        xtarg,ytarg = np.loadtxt(fits+'.coo')

        shim = shift(targ,(yref-ytarg,xref-xtarg),order=1)
        shval = alignIm(shim,ref,refmask) # first 
        shim = shift(shim,shval,order=1)
        targmed = np.median(shim[refmask])
        pf.writeto('aligned_'+fits, refmed*shim/targmed)
        data.append(shim)
    pf.writeto('Combo_%s.fits'%(fitslst[0].split('_')[0]),np.median(np.array(data),axis=0))
    if not os.path.exists('../Aligned'): os.mkdir('../Aligned')
    os.system('mv aligned*fits ../Aligned/')
    os.system('mv Combo_*fits ../Aligned/')


#------------------------------------------------------------------------------------------------
def alignIm(targ,ref,mask):
    xsh = np.arange(-5.,5.,1.)
    ysh = np.arange(-5.,5.,1.)
    shifts = []
    val = []
    for x in xsh:
        for y in ysh:
            shifts.append([x,y])
            val.append( np.sum( ( (shift(targ,(x,y),order=1) - ref)*mask) **2) )
    xsh1, ysh1 = shifts[np.argmin(val)][0], shifts[np.argmin(val)][1]

    xsh = np.arange(-.5+xsh1,.5+xsh1,0.1)
    ysh = np.arange(-.5+ysh1,.5+ysh1,0.1)
    shifts = []
    val = []
    for x in xsh:
        for y in ysh:
            shifts.append([x,y])
            val.append( np.sum( ((shift(targ,(x,y),order=1) - ref)*mask)**2) )
    print val[np.argmin(val)], shifts[np.argmin(val)]

    return shifts[np.argmin(val)]

#------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
