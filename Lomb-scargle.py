# This algorithm computes the Lomb-Scargle Periodogram.
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
from scipy.optimize import curve_fit  
wmin = 1.0e-04                 # minimum frequency 
wmax = 5.0e-02              # maximum frequency
freq = np.linspace(wmin,wmax,5000)
from astropy.stats import LombScargle
def func(x, a0, a1):
  return a0 + a1*x
#def func(x, a0, a1, a2, a3, a4, a5):
#  return a0 + a1*x + a2*x*x+ a3*x**3 + a4*x**4 + a5*x**5
#dy = 0.1
#t,A,B,C,D,E,F,dyy = np.loadtxt('methanol_peaks02.csv', unpack = True)
B,t = np.loadtxt('merged-file', unpack = True)

#ts = data[1:,i]
n = len(B)
m = len(t)
tsnew = np.empty(m)
dy = np.full(m,0.1)
p0 = [10, -0.001]
popt, pcov = curve_fit(func, t, B, p0)
#print(popt)
fitfunc2 = func(t, popt[0], popt[1])
for j in range(0,len(B)):
    tsnew[j] = B[j] - fitfunc2[j]
ls = LombScargle(t,tsnew,dy)
nu, power = ls.autopower()
power = LombScargle(t,tsnew,dy).power(freq)
fileout = open("merged-file","w")
print(1.0e0/freq[np.argmax(power)])
print(freq[np.argmax(power)])
pmax = (max(power))
#print(ls.false_alarm_probability(power.max()))
probabilities = [0.1, 0.05, 0.001]
#plevel = (ls.false_alarm_level(probabilities))
#plt.axes().set_aspect('1')
plt.clf()
plt.cla()
plt.xlim(wmin,wmax)
#plt.ylim(15,52)

plt.xlabel('Frequency (day$^{-1}$)',fontsize=18)
plt.ylabel('Power',fontsize=18)
#plt.plot(freq,power1,'--', color='b',label='OH 1667 MHz; P=216.1 days')
#plt.plot(freq,power3, '-.',color='r',label='OH 1665 MHz; P=215.9 days')
plt.plot(freq,power,'-', color='k',label='CH$_3$OH 6.7 GHz; P= 201 days')
#plt.axvline(freq[np.argmax(power2)], color='k', linestyle='--',linewidth=0.75)
#plt.axhline(plevel[2], color='k', linestyle='--',linewidth=0.75)
#plt.axhline(max(power3)/2, color='k', linestyle='--',linewidth=0.75)
#plt.text(0.006,y1,'Half max 1667 MHz')
#plt.text(0.006,y2,'Half max 1665 MHz')
#plt.legend(frameon=False)
#plt.savefig('LSpgrams.pdf')
#plt.savefig('LSpgrams.eps')
for i in range(n):
  fileout.write("%15.9e  %15.9e \n"% (freq[i],power[i]))
fileout.close()
plt.show()
plt.clf()
plt.cla()
plt.close()

