import sys
import numpy as np
from pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt

def stats(data):
    volini=data[0, 19]
    print("Volume:",(1-min(data[:, 19])/volini)*100);
    #hvini = 5.-data[0, 15]
    hvini = 25
    print("h ini",hvini)
    print("min z:",min(5.-data[:, 15]))
    print("Vent height:",(1-min(5.0-data[:, 15])/hvini)*100);
    
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['font.weight'] = 'bold'

# load data from simulations
ta45 = np.loadtxt('D:\\output\\expTa45\\fisiopacer.txt')
ta60 = np.loadtxt('D:\\output\\expTa60\\fisiopacer.txt')
ta75 = np.loadtxt('D:\\output\\expTa75\\fisiopacer.txt')
#0 CA->stats->contSave , 
#1 CA->params->dt, 2 CA->stats->avgVel, 3 CA->stats->tol, 4 CA->stats->err, 
#5 getActiveTensionNormalized(0, CA),
#6 getActiveTensionDiscrete( 0, CA),
#7 getV(31, CA), 8 getV(421, CA), 
#9 getV(195, CA), 10 getVdiscret(0, CA), 
#11 getPressurePercent(CA), 12 getPressureDiscrete(CA),
#13 CA->stats->min[0], 14 ->stats->min[1], 15CA->stats->min[2], 
#16 CA->stats->max[0], 17 CA->stats->max[1], 18 CA->stats->max[2])

print("45")
stats(ta45)
print("60")
stats(ta60)
print("75")
stats(ta75)

###
plt.plot(ta45[:, 0]/1000, 5.-ta45[:, 15], '-.', color='gray', label='$T_a=45kPa$', lw=2)
plt.plot(ta60[:, 0]/1000, 5.-ta60[:, 15], '-', color='black', label='$T_a=60kPa$', lw=2)
plt.plot(ta75[:, 0]/1000, 5.-ta75[:, 15], '--', color='gray', label='$T_a=75kPa$', lw=2)

sizFnt=20

plt.xlabel(u"Time (s)")
plt.ylabel(u"Ventricular height (mm)")
plt.grid(True, lw=0., color='0.75', zorder=0)
plt.legend(loc='best', numpoints=1)
plt.axis([0, 1.3, 19, 25])

#plt.xticks([], [])
#plt.yticks([], [])
plt.savefig('ventricularHeight.pdf', inches_bbox='tight')
plt.savefig('ventricularHeight.png', inches_bbox='tight')
#plt.show()
plt.clf()
volini=ta60[0, 19]
plt.plot(ta45[:, 0]/1000, 100*ta45[:, 19]/volini, '-.', color='gray', label='AT45', lw=2)
plt.plot(ta60[:, 0]/1000, 100*ta60[:, 19]/volini, '-', color='black', label='AT60', lw=2)
plt.plot(ta75[:, 0]/1000, 100*ta75[:, 19]/volini, '--', color='gray', label='AT75', lw=2)
plt.axis([0, 1.3, 97.3, 100.1])
plt.grid(True, lw=0., color='0.75', zorder=0)
plt.legend(loc='best', numpoints=1)

plt.xlabel(u"Time (s)")
plt.ylabel(u"Volume (%)")

plt.savefig('volume.pdf')
