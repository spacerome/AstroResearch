from astropy.io import fits
import astropy as ap
import glob
import numpy as np
import math
import scipy.signal as signal
import scipy.ndimage as ndimage
import time

start_time = time.time()
obs = glob.glob('*_obs.fits')
tru = glob.glob('*_tru.fits')
crm = glob.glob('*_crm.fits')
#####
obs1 = sorted(obs, key=lambda item: (int(item.partition(' ')[0])
                                  if item[0].isdigit()
                                  else float('inf'), item))
tru1 = sorted(tru, key=lambda item: (int(item.partition(' ')[0])
                                      if item[0].isdigit()
                                      else float('inf'), item))
crm1 = sorted(crm, key=lambda item: (int(item.partition(' ')[0])
                                            if item[0].isdigit()
                                            else float('inf'), item))
#####
nums = len(obs1)
xsize = 2048
ysize = 4096
q = []
q1 = []
q2 = []
img1 = []
img2 = []
bb = []
cc = []
crays = 0
cray1 = 0
cray2 = 0
for i in range(nums):
    s1 = fits.open(obs1[i])
    sdata1 = s1[1].data
    photflam1 = s1[1].header['photflam']
    b1 = fits.open(tru1[i])
    bdata1 = b1[0].data
    exptime1 = b1[0].header['exptime']
    cr1 = fits.open(crm1[i])
    cdata1 = cr1[0].data
    npix=0
    crays = np.where(cdata1)  
    # RAW - SKY (Gets Cray Values)
    # Np Where cdata1 = 1 (Cray)
    # if np where = 1 then we input the cray value
    # extra condition: IF value is -1 then we multiply it by -1 
    s1.close()
    b1.close()
    cr1.close()
bavg = np.average(bb)
print(bavg)
fract = npix/crays
perc = fract * 100
print(perc)
q = np.array(q)
q1=np.array(q1)
print(crays)
print(npix)
np.savetxt('cray1.txt',q,delimiter=',')
print("--- %s seconds ---" % (time.time() - start_time))
