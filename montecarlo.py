import numpy as np
import math as m
import matplotlib.pyplot as plt
import numpy.random as npr
import scipy.stats as s
from scipy.stats import ksone
import statistics as ss
from astropy.stats import SigmaClip
# Percentage = 6.1539224180665295%
# Gain = 1.000175
# Read Noise = 4.97749985
# Value = 446.214
# Avg = 43.016094
cray = np.loadtxt('cray.txt', delimiter = ',') #gets cosmic ray values from txt
mcray = max(cray) #Maximum cray value
flux = 43.016094
gain = 1.000175
n = gain * flux
g2 = gain ** 2
ng = n / g2
rn = 4.97749985
dna = ng - rn
sigma = m.sqrt(dna) #gets sigma
#print(sigma)
#plt.ylim(0,mcray)
#perc = 6.1539224180665295
perc = 6.1539224180665295
prob = perc/100
plt.close('all')
#plt.figure(0)
#plt.title('Number of Cosmic Rays for GRB050509B')
#plt.xlim(0,1500)
#plt.autoscale(enable=True,axis='y')
#plt.xlabel('Cosmic Ray Counts (DN)')
#plt.ylabel('Frequency')
#plt.hist(cray,bins='auto')
a = []
b = [] #mean
c = [] #median
d = [] #sigclip
aa = []
bb = [] #mean
cc = [] #median
dd = [] #sigclip
lm = [] #test
lm1 = [] #test
# COSMIC RAY REPEATS 10K TIMES
# % is what exp is 2k = 6% 1k = 3%
# Basically assume it is linear.
for j in range (10000):
    for k in range(3): #5 pixels
        q = npr.normal(flux,sigma)
        p = npr.uniform()
        if p < prob:
            ll = npr.choice(cray,1)
            while ll < 0:
                ll = npr.choice(cray,1)
            ll.tolist()
            a.append(ll)
        else:
            a.append(q)
    ma = np.mean(a)
    ms = np.median(a)
    #msc, up, low = s.sigmaclip(a,low=sigma,high=sigma) #SIGMA CLIP
    #msc = sigma_clip(a,sigma=sigma)
    b.append(ma)
    c.append(ms)
##    msc.tolist() #turns numpy array to python list
##    for i in range(len(msc)):
##        d.append(msc[i])
    a = []
 #NO COSMIC RAY
for j in range (10000):
    for k in range(3): #5 pixels
        qq = npr.normal(flux,sigma)
        aa.append(qq)
    ma1 = np.mean(aa)
    ms1 = np.median(aa)
    #msc1, upp,loww = s.sigmaclip(aa,low=sigma,high=sigma) #SIGMA CLIP
    #lm1.append(aa)
    bb.append(ma1)
    cc.append(ms1)
    #msc1.tolist() #turns numpy array to python list
##    for i in range(len(msc1)):
##        dd.append(msc1[i])
    aa = []
b = np.array(b) # Figure 1 Mean (CR)
c = np.array(c) # Figure 2 Median (CR)
#d = np.array(d) # Figure 3 SigClip (CR)
lm = np.array(lm)
#lm1 = np.array(lm1)
bb = np.array(bb) # Figure 4 Mean (No CR)
cc = np.array(cc) # Figure 5 Median (No CR)
#dd = np.array(dd) # Figure 6 SigClip (No CR)
sigclip = SigmaClip(sigma=sigma)
sigclip1 = SigmaClip(sigma=sigma)
##b = np.hstack(b)
##c = np.hstack(c)
##d = np.hstack(d)
##bb = np.hstack(bb)
##cc = np.hstack(cc)
##dd = np.hstack(dd)
print("Figure 1 STANDARD DEVIATION:",np.std(b),"MEAN:",np.mean(b))
print("Figure 2 STANDARD DEVIATION:",np.std(c),"MEAN:",np.mean(c))
print("Figure 3 STANDARD DEVIATION:",np.std(sigclip(b)),"MEAN:",np.mean(sigclip(b)))
print("Figure 4 STANDARD DEVIATION:",np.std(sigclip1(c)),"MEAN:",np.mean(sigclip1(c)))
print("Figure 5 STANDARD DEVIATION:",np.std(bb),"MEAN:",np.mean(bb))
print("Figure 6 STANDARD DEVIATION:",np.std(cc),"MEAN:",np.mean(cc))
plt.figure(1)
plt.title('Mean (Cosmic Rays) for n = 5')
plt.autoscale(enable=True,axis='both')
plt.xlabel('Counts (DN)')
plt.ylabel('Frequency')
plt.hist(b,bins='auto')
plt.figure(2)
plt.title('Median (Cosmic Rays) for n = 5')
plt.autoscale(enable=True,axis='both')
plt.xlabel('Counts (DN)')
plt.ylabel('Frequency')
plt.hist(c,bins='auto')
plt.figure(3)
plt.title('Sigma Clipped Mean (Cosmic Rays) for n = 5')
plt.autoscale(enable=True,axis='both')
plt.xlabel('Cosmic Ray Counts (DN)')
plt.ylabel('Frequency')
plt.hist(sigclip(b),bins='auto')
plt.figure(4)
plt.title('Sigma Clipped Median (No Cosmic Rays) for n = 5')
plt.autoscale(enable=True,axis='both')
plt.xlabel('Counts (DN)')
plt.ylabel('Frequency')
plt.hist(sigclip1(c),bins='auto')
plt.figure(5)
plt.title('Mean (No Cosmic Rays) for n = 5')
plt.autoscale(enable=True,axis='both')
plt.xlabel('Counts (DN)')
plt.ylabel('Frequency')
plt.hist(bb,bins='auto')
plt.figure(6)
plt.title('Median (No Cosmic Rays) for n = 5')
plt.autoscale(enable=True,axis='both')
plt.xlabel('Counts (DN)')
plt.ylabel('Frequency')
plt.hist(cc,bins='auto')
plt.show()
