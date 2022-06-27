#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

data15581 = np.loadtxt('15581_500x500_seg.dat')
data15588 = np.loadtxt('15588_500x500_seg.dat')
data15604 = np.loadtxt('15604_500x500_seg.dat')
data5479  = np.loadtxt('5479_500x500_seg.dat')
data5548  = np.loadtxt('5548_500x500_seg.dat')
data9005  = np.loadtxt('9005-2_500x500_seg.dat')

plt.figure(figsize = (10, 8), dpi = 300)
plt.rc('font', size = 18)
plt.plot(data15581[:,0], linewidth = 2.0)
plt.plot(data15588[:,0], linewidth = 2.0)
plt.plot(data15604[:,0], linewidth = 2.0)
plt.plot(data5479[:,0],  linewidth = 2.0)
plt.plot(data5548[:,0],  linewidth = 2.0)
plt.plot(data9005[:,0],  linewidth = 2.0)
plt.xlabel('Correlation length $r$')
plt.ylabel('Surface-surface function $F_{ss}(r)$')
plt.legend(['Sandstone 1', 'Sandstone 2', 'Sandstone 3', 'Carbonate 1', 'Carbonate 2', 'Carbonate 3'])
plt.savefig('../images/KIK-ss.png')

plt.figure(figsize = (10, 8), dpi = 300)
plt.rc('font', size = 18)
plt.plot(data15581[:,1], linewidth = 2.0)
plt.plot(data15588[:,1], linewidth = 2.0)
plt.plot(data15604[:,1], linewidth = 2.0)
plt.plot(data5479[:,1],  linewidth = 2.0)
plt.plot(data5548[:,1],  linewidth = 2.0)
plt.plot(data9005[:,1],  linewidth = 2.0)
plt.xlabel('Correlation length $r$')
plt.ylabel('Surface-void function $F_{sv}(r)$')
plt.legend(['Sandstone 1', 'Sandstone 2', 'Sandstone 3', 'Carbonate 1', 'Carbonate 2', 'Carbonate 3'])
plt.savefig('../images/KIK-sv.png')
