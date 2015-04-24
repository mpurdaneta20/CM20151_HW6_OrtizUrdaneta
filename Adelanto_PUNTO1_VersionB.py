# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline

# <codecell>

def func_velo(v): 
    return v

def func_acel(r_cubo,coord1,coord2,m): #coord1 y coord 2 pueden ser x, y o z.
    G=1
    return((G*m*coord1-coord2/np.abs(r1-r2)**3))
    #return (G*m/np.abs(((x_i-x_j)**2)+((y_i-y_j)**2)+((z_i-z_j)**2))**3)*(((x_i-x_j)**2)+((y_i-y_j)**2)+((z_i-z_j)))#aceleracion

def r_cubo(x_i,y_i,z_i,x_j,y_j,z_j):
    return np.abs((x_i-x_j)**2+(y_i-y_j)**2+(z_i-z_j)**2)

# <codecell>

N=0
n_points=100

x_t = zeros(n_points)
y_t = zeros(n_points)
z_t = zeros(n_points)

x_s = zeros(n_points)
y_s = zeros(n_points)
z_s = zeros(n_points)

Vx_t = zeros(n_points)
Vy_t = zeros(n_points)
Vz_t = zeros(n_points)

Vx_s = zeros(n_points)
Vy_s = zeros(n_points)
Vz_s = zeros(n_points)


# <codecell>

h=1
def RungeKuttaFourthOrderStep(Vx_t, Vy_t,Vz_t,Vx_s,Vy_s,Vz_s, x_t, y_t, z_t, x_s, y_s, z_s, N): #tiempo??
    t= N*h
    ms=2
    mt=1
    G=1
    
    K1_Tx=func_velo(Vx_t)
    K1_Ty=func_velo(Vy_t)
    K1_Tz=func_velo(Vz_t)
    
    r_cubo1=r_cubo(x_t,x_s,y_t,y_s,z_t,z_s)
    
    K1_TVx=func_acel(r_cubo1,x_t,x_s,ms)
    K1_TVy=func_acel(r_cubo1,y_t,y_s,ms)
    K1_TVz=func_acel(r_cubo1,z_t,z_s,ms)
    
    K1_Sx=func_velo(Vx_s)
    K1_Sy=func_velo(Vy_s)
    K1_Sz=func_velo(Vz_s)
    
    r_cubo2=r_cubo(x_s,x_t,y_s,y_t,z_s,z_t)
    
    K1_SVx=func_acel(r_cubo2,x_s,x_t,mt)
    K1_SVy=func_acel(r_cubo2,y_s,y_t,mt)
    K1_SVz=func_acel(r_cubo2,z_s,z_t,mt)
    
    #K2
    K2_Tx=func_velo(Vx_t+(K1_Tx*h/2))
    K2_Ty=func_velo(Vy_t+(K1_Ty*h/2))
    K2_Tz=func_velo(Vz_t+(K1_Tz*h/2))
    
    r_cubo1=(r_cubo(x_t,x_s,y_t,y_s,z_t,z_s))*(h/2)**(3/2)
    
    K2_TVx=func_acel(r_cubo1,x_t+(K1_Tx*h/2),x_s+(K1_Tx*h/2),ms)
    K2_TVy=func_acel(r_cubo1,y_t+(K1_Ty*h/2),y_s+(K1_Ty*h/2),ys,ms)
    K2_TVz=func_acel(r_cubo1,z_t+(K1_Tz*h/2),z_s+(K1_Tz*h/2),ms)
    
    K2_Sx=func_velo(Vx_s+(K1_Tx*h/2))
    K2_Sy=func_velo(Vy_s+(K1_Ty*h/2))
    K2_Sz=func_velo(Vz_s+(K1_Tz*h/2))
    
    r_cubo2=(r_cubo(x_t,x_s,y_t,y_s,z_t,z_s))*(h/2)**(3/2)
    
    K2_SVx=func_acel(r_cubo2,x_s+(K1_Tx*h/2),x_t+(K1_Tx*h/2),mt)
    K2_SVy=func_acel(r_cubo2,y_s+(K1_Ty*h/2),y_t+(K1_Ty*h/2),mt)
    K2_SVz=func_acel(r_cubo2,z_s+(K1_Tz*h/2),z_t+(K1_Tz*h/2),mt)
    
    #K3
    K3_Tx=func_velo(Vx_t+(K2_Tx*h/2))
    K3_Ty=func_velo(Vy_t+(K2_Ty*h/2))
    K3_Tz=func_velo(Vz_t+(K2_Tz*h/2))
    
    r_cubo1=(r_cubo(x_t,x_s,y_t,y_s,z_t,z_s))*(h/2)**(3/2)
    
    K3_TVx=func_acel(r_cubo1,x_t+(K2_Tx*h/2),x_s+(K2_Tx*h/2),ms)
    K3_TVy=func_acel(r_cubo1,y_t+(K2_Ty*h/2),y_s+(K2_Ty*h/2),ys,ms)
    K3_TVz=func_acel(r_cubo1,z_t+(K2_Tz*h/2),z_s+(K2_Tz*h/2),ms)
    
    K3_Sx=func_velo(Vx_s+(K2_Tx*h/2))
    K3_Sy=func_velo(Vy_s+(K2_Ty*h/2))
    K3_Sz=func_velo(Vz_s+(K2_Tz*h/2))
    
    r_cubo2=(r_cubo(x_t,x_s,y_t,y_s,z_t,z_s))*(h/2)**(3/2)
    
    K3_SVx=func_acel(r_cubo2,x_s+(K2_Tx*h/2),x_t+(K2_Tx*h/2),mt)
    K3_SVy=func_acel(r_cubo2,y_s+(K2_Ty*h/2),y_t+(K2_Ty*h/2),mt)
    K3_SVz=func_acel(r_cubo2,z_s+(K2_Tz*h/2),z_t+(K2_Tz*h/2),mt)
    
    #K4
    K4_Tx=func_velo(Vx_t+(K3_Tx*h/2))
    K4_Ty=func_velo(Vy_t+(K3_Ty*h/2))
    K4_Tz=func_velo(Vz_t+(K3_Tz*h/2))
    
    r_cubo1=(r_cubo(x_t,x_s,y_t,y_s,z_t,z_s))*(h/2)**(3/2)
    
    K4_TVx=func_acel(r_cubo1,x_t+(K3_Tx*h/2),x_s+(K3_Tx*h/2),ms)
    K4_TVy=func_acel(r_cubo1,y_t+(K3_Ty*h/2),y_s+(K3_Ty*h/2),ys,ms)
    K4_TVz=func_acel(r_cubo1,z_t+(K3_Tz*h/2),z_s+(K3_Tz*h/2),ms)
    
    K4_Sx=func_velo(Vx_s+(K3_Tx*h/2))
    K4_Sy=func_velo(Vy_s+(K3_Ty*h/2))
    K4_Sz=func_velo(Vz_s+(K3_Tz*h/2))
    
    r_cubo2=(r_cubo(x_t,x_s,y_t,y_s,z_t,z_s))*(h/2)**(3/2)
    
    K4_SVx=func_acel(r_cubo2,x_s+(K3_Tx*h/2),x_t+(K3_Tx*h/2),mt)
    K4_SVy=func_acel(r_cubo2,y_s+(K3_Ty*h/2),y_t+(K3_Ty*h/2),mt)
    K4_SVz=func_acel(r_cubo2,z_s+(K3_Tz*h/2),z_t+(K3_Tz*h/2),mt)
    
    
    #Actualizacion
    average_Tx= (1.0/6.0)*(k1_Tx + 2.0*k2_Tx + 2.0*k3_Tx + k4_Tx)+(x_t)
    average_Ty= (1.0/6.0)*(k1_Ty + 2.0*k2_Ty + 2.0*k3_Ty + k4_Ty)+(y_t)
    average_Tz= (1.0/6.0)*(k1_Tz + 2.0*k2_Tz + 2.0*k3_Tz + k4_Tz)+(z_t)
    
    average_VTx= (1.0/6.0)*(k1_TVx + 2.0*k2_TVx + 2.0*k3_TVx + k4_TVx)+(Vx_t)
    average_VTy= (1.0/6.0)*(k1_TVy + 2.0*k2_TVy + 2.0*k3_TVy + k4_TVy)+(Vy_t)
    average_VTz= (1.0/6.0)*(k1_TVz + 2.0*k2_TVz + 2.0*k3_TVz + k4_TVz)+(Vz_t)
    
    average_Sx= (h/6.0)*(k1_Sx + 2.0*k2_Sx + 2.0*k3_Sx + k4_Sx)+(s_x)
    average_Sy= (h/6.0)*(k1_Sy + 2.0*k2_Sy + 2.0*k3_Sy + k4_Sy)+(s_y)
    average_Sz= (h/6.0)*(k1_Sz + 2.0*k2_Sz + 2.0*k3_Sz + k4_Sz)+(s_z)
    
    average_VSx= (h/6.0)*(k1_SVx + 2.0*k2_SVx + 2.0*k3_SVx + k4_SVx)+(Vx_s)
    average_VSy= (h/6.0)*(k1_SVy + 2.0*k2_SVy + 2.0*k3_SVy + k4_SVy)+(Vy_s)
    average_VSz= (h/6.0)*(k1_SVz + 2.0*k2_SVz + 2.0*k3_SVz + k4_SVz)+(Vz_s)

    
    x_t=average_Tx
    y_t=average_Ty
    z_t=average_Tz
    x_s=average_Sx
    y_s=average_Sy
    z_s=average_Sz
    
    Vx_t=average_VTx
    Vy_t=average_VTy
    Vz_t=average_VTz
    Vx_s=average_VSx
    Vy_s=average_VSy
    Vz_s=average_VSz
    
    
    return Vx_t, Vy_t,Vz_t,Vx_s,Vy_s,Vz_s, x_t, y_t, z_t, x_s, y_s, z_s, t

# <codecell>

x_t[0]=1
y_t[0]=0
z_t[0]=0
x_s[0]=0
y_s[0]=0
z_s[0]=0

Vx_t[0]=0
Vy_t[0]=0
Vz_t[0]=0
Vx_s[0]=0
Vy_s[0]=0
Vz_s[0]=0
    
for i in range(1,n_points):
    N=i
    x_t[i],y_t[i],z_t[i],x_s[i],y_s[i], z_s[i], Vx_t[i], Vy_t[i], Vz_t[i], Vx_s[i], Vy_s[i], Vz_s[i] = RungeKuttaFourthOrderStep(Vx_t[i-1], Vy_t[i-1],Vz_t[i-1],Vx_s[i-1],Vy_s[i-1],Vz_s[i-1],x_t[i-1],y_t[i-1],z_t[i-1], x_s[i-1], y_s[i-1], z_s[i-1], N)

# <codecell>


