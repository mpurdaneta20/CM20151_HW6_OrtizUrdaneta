# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Solución en 1D 

# <codecell>

from matplotlib import*
from pylab import*
from matplotlib import figure
from numpy import*

# <codecell>

#Se tienen los siguientes datos los cuales se mantendran constantes a lo largo del metodo
p=1000  #(ro)  densidad
cp=4181 #calor especifico agua
D=0.1  # coeficiente de dispersion 
d_H=100 #(delta_H)  calor de reaccion
Ra=1  #radio maximo del tubo
ra=0.5 #Radio tubo interno
L=10 #longitud total del tubo. Rango de [0,L]
E=10 #modulo de Young
C0=1 #concentracion inicial
T0=20 #Temp inicial en Kelvin
lam=0.5  #Lambda
vmax=0.5 #velocidad máxima
Cin=0.8
Tin=10

# <codecell>

# Función r(c,T) 
k0 = 10 #constante cinética 
E = 10 
R = 9 #consante de los gases
c=1
T=1 

def r(c,T):
    rcT=k0*exp(-E/(R*T))*c**2
    return rcT

#Función de Velocidad superficial del fluido
def v(ra):
    vr=vmax*(1-(ra/Ra)**2)
    return vr

# <codecell>

#arreglos para z , y tiempo
n_points=100
z = linspace(0.0,L,n_points)
t= linspace(0.0,1.0,n_points)

#Se definen los deltas
delta_z =np.abs(z[1]-z[0])
delta_t = 0.0005

#relación de alpha con delta_x y delta_t para obtener un parámetro
sigma = 0.05
#delta_t = sigma*delta_z/vmax
alpha=(-vmax*(delta_t/delta_z))
betha=(delta_t/delta_z**2)*D
alpha2=(-vmax*(delta_t/delta_z))
betha2=(delta_t/delta_z**2)*(lam/(p*cp))
betha3=(d_H/(p*cp))

print alpha
print betha
print delta_t
print delta_z

# <codecell>

#Funcion incial
c1 = Cin*(zeros(n_points))
t1= Tin*(ones(n_points))

# <codecell>

#c0 y t0 inicial=concentracion inicial y temp. inicial
c1[:] = C0
c1[0] = Cin
c1[-1] = c1[-2]

t1[:] = T0
t1[0] = Tin
t1[-1] = t1[-2]

#create a new variable to hold the previous value
c_past = c1.copy()
t_past = t1.copy()

plot(z, c_past)
plot(z, c1)

# <codecell>

plot(z, t_past)
plot(z, t1)

# <codecell>

n_tpoints=100
for i in range(n_tpoints):  # loop over time
    c_past = c1.copy()
    t_past = t1.copy()
    for i in range(1,n_points-1): #loop over space
        c1[i] = c_past[i-1]*(alpha+betha)+c_past[i]*(1-alpha-(2*betha))+(c_past[i+1]*betha)-(delta_t*r(c_past[i],t_past[i]))
        t1[i] = (t_past[i]*(1-alpha2-(2*betha2)))+(t_past[i-1]*(betha2-alpha2))+(t_past[i+1]*(betha2))+(betha3*(r(c_past[i],t_past[i]))) 
#print u1

# <codecell>

plot(z,c1)
print c1

# <codecell>

plot (z,t1)
print t1

# <codecell>


# <codecell>


