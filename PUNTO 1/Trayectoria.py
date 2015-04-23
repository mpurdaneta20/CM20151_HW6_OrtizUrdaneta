# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
from matplotlib import*
from pylab import*
from matplotlib import figure



o1= np.loadtxt("orbitas-1yrs.txt", skiprows=1)
xs = o1[:,0]
ys = o1[:,1]
xt= o1[:,2]
yt= o1[:,3]
xl= o1[:,4]
yl= o1[:,5]
xa= o1[:,6]
ya= o1[:,7]

#Genera la grafica en 2D - Orbita en 1 año
plt.figure()
plt.plot (xs,ys,label=r"$Sol$")
plt.plot (xl,yl,label=r"$Luna$")
plt.plot(xt,yt,label=r"$Tierra$")
plt.plot(xa,ya,label=r"$Asteroide$")
plt.title(r'Trayectoria de los 4 Astros (t=1)', fontsize=20)
plt.xlabel(r'$Eje\ X$',fontsize=20)
plt.ylabel(r'$Eje\ y$',fontsize=20)
plt.grid()
plt.ylim([-1.2,1.2])
plt.xlim([-3,3])
plt.legend(loc=0, fontsize=10)
plt.savefig("orbitas-1yrs.png", format='png', bbox_inches='tight')

o1000 = np.loadtxt("orbitas-1000yrs.txt", skiprows=1)
xs = o1000[:,0]
ys = o1000[:,1]
xt= o1000[:,2]
yt= o1000[:,3]
xl= o1000[:,4]
yl= o1000[:,5]
xa= o1000[:,6]
ya= o1000[:,7]

#Genera la grafica en 2D - Orbita en 1000 años
plt.figure(figsize=(15, 10))
plt.plot (xs,ys,label=r"$Sol$")
plt.plot (xl,yl,label=r"$Luna$")
plt.plot(xt,yt,label=r"$Tierra$")
plt.plot(xa,ya,label=r"$Asteroide$")
plt.title(r'Trayectoria de los 4 Astros (t=2)', fontsize=20)
plt.xlabel(r'$Eje\ X$',fontsize=20)
plt.ylabel(r'$Eje\ y$',fontsize=20)
plt.grid()
plt.ylim([-10,10])
plt.xlim([-10,10])
plt.legend(loc=0, fontsize=10)
plt.savefig("orbitas-1000yrs.png", format='png', bbox_inches='tight')


