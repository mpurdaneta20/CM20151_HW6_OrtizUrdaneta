# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline

# <codecell>

rt=[0,0,2]
rs=[0,0,1]

vt=[1,1,1]
vs=[2,2,2]

# <codecell>

#Se obtiene la derivada de la posicion lo cual genera la velocidad
def Posicion_punto(V,salto): 
    V[0]=V[0]+salto[0]
    V[1]=V[1]+salto[1]
    V[2]=V[2]+salto[2]
    return V

def Norma_cubo(r_i,r_j):
    return np.sqrt((r_i[0]-r_j[0])**2+(r_i[1]-r_j[1])**2+(r_i[2]-r_j[2])**2)

def salto(K,paso):
    K[0]=K[0]*paso
    K[1]=K[1]*paso
    K[2]=K[2]*paso
    return K

#Se obtiene la derivada de la velocidad lo cual genera la aceleracion
def Velocidad_punto_x(r_i,r_j,m,salto): #coord1 y coord 2 pueden ser x, y o z.
    G=1.0
    a=[0,0,0]
    a_x=(G*(r_i[0]-r_j[0])*m)/Norma_cubo(r_i,r_j) + salto
    a[0]=a_x
    a_y=(G*(r_i[1]-r_j[1])*m)/Norma_cubo(r_i,r_j)+ salto
    a[1]=a_y
    a_z=(G*(r_i[2]-r_j[2])*m)/Norma_cubo(r_i,r_j)+ salto
    a[2]=a_z
    return a

# <codecell>

r_i=[0,0,0]
r_j=[1,1,1]
Norma_cubo(r_i,r_j)

# <codecell>

n_points=100

r_T = np.zeros((n_points,3))
r_S = np.zeros((n_points,3))

v_T =np.zeros((n_points,3))
v_S =np.zeros((n_points,3))

# <codecell>

#Constane que determina el delta t
h=1
#Definicion de RungeKutta
def RungeKuttaFourthOrderStep(r_T,r_S,m_T,m_S,v_T,v_S, N): #tiempo??
    t= N*h
    
    s=[0,0,0]
    
    K1_PosT=Posicion_punto(v_T,s)
    K1_VelT=Velocidad_punto(r_T,r_S,m_S,s)

    K1_PosS=Posicion_punto(v_S,s)
    K1_VelS=Velocidad_punto(r_S,r_T,m_T,s)    
    
    #K2
    s=salto(K1_PosT,h/2)
    K2_PosT=Posicion_punto(v_T,s)
    K2_VelT=Velocidad_punto(r_T,r_S,m_S,s)
    
    K2_PosS=Posicion_punto(v_S,s)
    K2_VelS=Velocidad_punto(r_S,r_T,m_T,s)    
    
    #K3
    s=salto(K2_PosT,h/2)
    K3_PosT=Posicion_punto(v_T,s)
    K3_VelT=Velocidad_punto(r_T,r_S,m_S,s)
    
    K3_PosS=Posicion_punto(v_S,s)
    K3_VelS=Velocidad_punto(r_S,r_T,m_T,s)
    
    #K4
    s=salto(K3_PosT,h/2)
    K4_PosT=Posicion_punto(v_T,s)
    K4_VelT=Velocidad_punto(r_T,r_S,m_S,s)
    
    K4_PosS=Posicion_punto(v_S,s)
    K4_VelS=Velocidad_punto(r_S,r_T,m_T,s)

    
    
    return K1_PosT

# <codecell>

r_T[0]=np.array([1,1,0])
r_S[0]=[0,0,0]

v_T[0]=[1,0,0]
v_S[0]=[0,0,0]

m_T=1
m_S=1


for i in range(1,2):
    N=2
    r_T[i] = RungeKuttaFourthOrderStep(r_T[i-1],r_S[i-1],m_T,m_S,v_T[i-1],v_S[i-1], N)
    print r_T

# <codecell>


