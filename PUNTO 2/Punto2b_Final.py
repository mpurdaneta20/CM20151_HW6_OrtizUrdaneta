


from matplotlib import*
from pylab import*
from matplotlib import figure
import numpy as np
from scipy import meshgrid
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Constantes

p=1000  #(ro)  densidad
cp=4181 #calor especifico agua
D=10**(-4) # coeficiente de dispersion 
d_H=-50 #(delta_H)  calor de reaccion
Ra=1  #radio maximo del tubo
ra=0.8 #Radio tubo interno
L=10 #longitud total del tubo. Rango de [0,L]
E=10 #modulo de Young
C0=0 #concentracion inicial
T0=10#Temp inicial en Kelvin
lam=1000  #Lambda
vmax=0.001#velocidad máxima
h = 2   # h de función de intercambio termico en la pared
Tw= 80

Cin=1
Tin=20

#Arrelgos para y y z
n_points=5
r = np.linspace(0.0,Ra,n_points)
z = np.linspace(0.0,L,n_points)

# Función r(c,T) 
k0 = 10 #constante cinética 
E = 10 

c=1
T=1 

# Velocidad- radio
def func_r(c,T):
    rcT=k0*exp(-E/(8.3*T))*c**2
    return rcT


#Función de Velocidad superficial del fluido
def func_v(ra):
    vr=vmax*(1-(ra/Ra)**2)
    return vr

#Funciones de concentraciones y temperaturas iniciales, hay una temperatura ambiente constante en todo el fluido
U1_0 = ones((n_points, n_points))
U2_0 = (ones((n_points, n_points))*16)

#Funciones de concentración (U1)(zxr) y Temperatura (U2)(zxr)
U1 = ones((n_points, n_points))
U2 = ones((n_points, n_points))

#Funciones de concentracion y temperatura pasadas
U1_2 = np.ones((n_points,n_points)) ##
U2_2 = np.ones((n_points,n_points))

#print U1
#print U2


# In[141]:


U1=U1_0.copy()
U2=U2_0.copy()

#Condiciones de Frontera-
# Simetria Radial sin transf. de Masa r=0
U1[0,:] = U1[1,:]
U2[0,:] = U2[1,:]

#Simetria Radial sin transf. de Masa r=R
U1[-1,:] = U1[-2,:]

# Condición de Intercambio Termico en la pared  r=R
U2[-1,:] = (h/lam)*(Tw-U2[-1,:])+U2[-2,:]

# Concentracion y Temperatura constantes al inicio del tubo, en z=0
U1[:,0] = U1_0[:,0]
U2[:,0] = U2_0[:,0]

# Difusion cero a la salida, en z=L //C(L,r)=C(L-1,r) // T(L,r)=T(L-1,r)
U1[:,-1] = U1[:,-2]
U2[:,-1] = U2[:,-2]


# In[142]:

# Simetria Radial sin transf. de Masa r=0

U1_0[0,:] = U1_0[1,:]
U2_0[0,:] = U2_0[1,:]

#Simetria Radial sin transf. de Masa r=R
U1_0[-1,:] = U1_0[-2,:]

# Condición de Intercambio Termico en la pared  r=R
U2_0[-1,:] = (h/lam)*(Tw-U2_0[-1,:])+U2_0[-2,:]
U2_0[:,-1] = (h/lam)*(Tw-U2_0[:,-1])+U2_0[:,-2]

# Concentracion y Temperatura constantes al inicio del tubo, en z=0
U1_0[:,0] = U1_0[:,0]
U2_0[:,0] = U2_0[:,0]

# Difusion cero a la salida, en z=L //C(L,r)=C(L-1,r) // T(L,r)=T(L-1,r)

U1_0[:,-1] = U1_0[:,-2]
U2_0[:,-1] = U2_0[:,-2]
U2_0[-1,:] = U2_0[-2,:]
       


# In[143]:

#plot representación de condiciones iniciales
fig=plt.figure()
ax=fig.gca(projection='3d')
Z,R=np.meshgrid(z,r) 
#representa la concentración inicial c = U1
wire1=ax.plot_wireframe(Z,R,U1[:],cmap=cm.coolwarm)
#representa la Temperatura inicial T=U2
wire2=ax.plot_wireframe(Z,R,U2[:],cmap=cm.coolwarm)


# In[144]:

#constantes para simplificar

#Delta r y Delta z
delta_r = r[1]-r[0]
delta_z = z[1]-z[0]
delta_t = 0.05


alpha_1=((delta_t/delta_z))
alpha_2=(D*(delta_t/(delta_z**2)))
alpha_3=(D*(delta_t/(delta_r**2)))
alpha_4=(D*(delta_t/(ra*(delta_r)))) #verificar el ra

betha_1=((delta_t/delta_z))
betha_2=((delta_t*lam)/(p*cp*(delta_z)**2))
betha_3=((delta_t*lam)/(p*cp*(delta_r)**2))
betha_4=((delta_t*lam)/(p*cp*ra))#verificar este ra
betha_5=((d_H*delta_t)/(p*cp*(delta_r**2)))

print delta_r
print delta_z
print alpha_1
print alpha_2
print alpha_3
print alpha_4
print betha_1
print betha_2
print betha_3
print betha_4
print betha_5


# In[145]:

n_time = 100

for n in range(n_time+1): ##loop across number of time steps
    U1_2 = U1.copy() #cpast y cnow
    U2_2 = U2.copy()
    
    #Se despejó para t+1 haciendo la derivada -backwards-ver doc anexo con las ecuaciones
    
    
    U1[1:-1,1:-1] = alpha_1*(-1*func_v(ra))*(U1_2[1:-1,1:-1]-U1_2[1:-1,0:-2])+(alpha_2*((U1_2[1:-1,0:-2])-(2*(U1_2[1:-1,1:-1]))+(U1_2[1:-1,0:-2])))+    (alpha_3*((U1_2[2:,1:-1])-(2*(U1_2[1:-1,1:-1])))+(U1_2[0:-2,1:-1]))+(alpha_4*((U1_2[1:-1,1:-1])-(U1_2[0:-2,1:-1])))+(U1_2[1:-1,1:-1])-(delta_t)*func_r(U1_2[1:-1,1:-1],(U2_2[1:-1,1:-1]))
    
    U2[1:-1,1:-1] = betha_1*(-1*func_v(ra))*(U2_2[1:-1,1:-1]-U2_2[1:-1,0:-2])+(betha_2*((U2_2[1:-1,0:-2])-(2*(U2_2[1:-1,1:-1]))+(U2_2[1:-1,0:-2])))+    (betha_3*((U2_2[2:,1:-1])-(2*(U2_2[1:-1,1:-1])))+(U2_2[0:-2,1:-1]))+(betha_4*((U2_2[1:-1,1:-1])-(U2_2[0:-2,1:-1])))+(U2_2[1:-1,1:-1])-(betha_5)*(func_r(U1_2[1:-1,1:-1],(U2_2[1:-1,1:-1])))

    # Simetria Radial sin transf. de Masa r=0
U1[0,:] = U1[1,:]
U2[0,:] = U2[1,:]

#Simetria Radial sin transf. de Masa r=R
U1[-1,:] = U1[-2,:]

# Condición de Intercambio Termico en la pared  r=R
U2[-1,:] = (h/lam)*(Tw-U2[-1,:])+U2[-2,:]

# Concentracion y Temperatura constantes al inicio del tubo, en z=0
U1[:,0] = U1_0[:,0]
U2[:,0] = U2_0[:,0]

# Difusion cero a la salida, en z=L //C(L,r)=C(L-1,r) // T(L,r)=T(L-1,r)
U1[:,-1] = U1[:,-2]
U2[:,-1] = U2[:,-2]

print U1
print U2
   
    


# In[149]:

fig=plt.figure()
ax=fig.gca(projection='3d')
Z,R=np.meshgrid(z,r) 
wire1=ax.plot_wireframe(Z,R,U1[:],cmap=cm.coolwarm)
savefig('Concentracion', format='png')

fig=plt.figure()
ax=fig.gca(projection='3d')
Z,R=np.meshgrid(z,r) 
wire1=ax.plot_wireframe(Z,R,U2[:],cmap=cm.coolwarm)
savefig('Temperatura', format='png')


# In[148]:

fig=plt.figure()
ax=fig.gca(projection='3d')
Z,R=np.meshgrid(z,r) 
wire1=ax.plot_wireframe(Z,R,U1[:],cmap=cm.coolwarm)
wire2=ax.plot_wireframe(Z,R,U2[:],cmap=cm.coolwarm)
savefig('Temperatura y Concentracion', format='png')


# In[ ]:



