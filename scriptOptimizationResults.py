# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 09:01:28 2017

@author: smullall
"""
import matplotlib.pyplot as plt

x=(79.9458,77.4996,74.9807,84.3753,83.1567,83.1567,82.8922)
y=(0.42133,0.26333,0.19092,3.4299, 1.3817,1.1246,1.2772)
lppdv=(1.98634, 1.74209,1.61313,3.33451,2.51208, 2.56355,2.66403)
messes=( 0.887778,0.976018,1.56767,0.953623,0.854554,0.874098,5.4593)
dvmodval1=(1.32081,6,3.28743,6,6,6,4.97862)

plt.plot(y,x,'bo')

plt.xlabel('Inefficiency: FP rate')
plt.ylabel('Completeness: TP rate')
plt.xlim((0,1.5))
plt.ylim((74,85))


plt.title('based on r62273')
#plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r62273/steveOptimizationResults-zoom.png')

plt.figure(figsize=(8,9))

plt.subplot(311)
plt.plot(y,lppdv,'ro')
plt.xlabel('Inefficiency: FP rate')
plt.ylabel('lppdv thresh')

plt.subplot(312)
plt.plot(y,messes,'r^')
plt.xlabel('Inefficiency: FP rate')
plt.ylabel('ses to mes thresh')

plt.subplot(313)
plt.plot(y,dvmodval1,'rv')
plt.xlabel('Inefficiency: FP rate')
plt.ylabel('dvmodval1 thresh')

plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r62273/steveOptimization-thresholds1.png')


#%%

