


import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
import sys
from scipy.integrate import tplquad
from mpl_toolkits.mplot3d import Axes3D


uu=sys.argv[1]
NR=int(input("Number of dipole on one row?"))
a=float(input("cube size in nm?"))
nombre=6*(NR**2)
#b=a/5

with open(uu, "w") as f:
    f.write("")

phi=np.zeros(nombre)
theta=np.zeros(nombre)
psi=np.zeros(nombre)
x=np.zeros(nombre)
y=np.zeros(nombre)
z=np.zeros(nombre)

#***************************************

for i in range(0,NR*NR):    
    phi[i]=np.pi/2
    theta[i]=np.pi/2
    psi[i]=0

for i in range(NR*NR,2*NR*NR):    
    phi[i]=3*np.pi/2
    theta[i]=np.pi/2
    psi[i]=0

for i in range(2*NR*NR,3*NR*NR):    
    phi[i]=0
    theta[i]=np.pi/2
    psi[i]=0

for i in range(3*NR*NR,4*NR*NR):    
    phi[i]=2*np.pi/2
    theta[i]=np.pi/2
    psi[i]=0

for i in range(4*NR*NR,5*NR*NR):    
    phi[i]=0
    theta[i]=0
    psi[i]=0

for i in range(5*NR*NR,6*NR*NR):    
    phi[i]=0
    theta[i]=np.pi
    psi[i]=0

kkk=0
for i in range(0,NR):
    for j in range(0,NR):
        x[kkk]=(a/NR)*(i)-(a/2)+(a/(2*NR))
        y[kkk]=(a/2)
        z[kkk]=(a/NR)*(j)-(a/2)+(a/(2*NR))
        kkk=kkk+1

for i in range(0,NR):
    for j in range(0,NR):
        x[kkk]=(a/NR)*(i)-(a/2)+(a/(2*NR))
        y[kkk]=-a/2
        z[kkk]=(a/NR)*(j)-(a/2)+(a/(2*NR))
        kkk=kkk+1

for i in range(0,NR):
    for j in range(0,NR):
        x[kkk]=a/2
        y[kkk]=(a/NR)*(i)-(a/2)+(a/(2*NR))
        z[kkk]=(a/NR)*(j)-(a/2)+(a/(2*NR))
        kkk=kkk+1

for i in range(0,NR):
    for j in range(0,NR):
        x[kkk]=-a/2
        y[kkk]=(a/NR)*(i)-(a/2)+(a/(2*NR))
        z[kkk]=(a/NR)*(j)-(a/2)+(a/(2*NR))
        kkk=kkk+1
        
for i in range(0,NR):
    for j in range(0,NR):
        x[kkk]=(a/NR)*(i)-(a/2)+(a/(2*NR))
        y[kkk]=(a/NR)*(j)-(a/2)+(a/(2*NR))
        z[kkk]=a/2
        kkk=kkk+1
        
for i in range(0,NR):
    for j in range(0,NR):
        x[kkk]=(a/NR)*(i)-(a/2)+(a/(2*NR))
        y[kkk]=(a/NR)*(j)-(a/2)+(a/(2*NR))
        z[kkk]=-a/2
        kkk=kkk+1

xmax=np.amax(x)
ymax=np.amax(y)
zmax=np.amax(z)
b=(xmax+ymax+zmax)/10      
print("nombre de dipole")
print(nombre)
for j in range(nombre):
    with open(uu, "a") as f:
        valeurs = "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" % (phi[j],theta[j],psi[j],x[j],y[j],z[j]) 
        f.write(valeurs)
  #  xiO=np.zeros(nombre)
  #  yiO=np.zeros(nombre)
   # ziO=np.zeros(nombre)

    #xiH1=np.zeros(nombre)
    #yiH1=np.zeros(nombre)
    #ziH1=np.zeros(nombre)

   # xiH2=np.zeros(nombre)
    #yiH2=np.zeros(nombre)
    #ziH2=np.zeros(nombre)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.quiver(x,y,z,b*np.cos(phi)*np.sin(theta),b*np.sin(phi)*np.sin(theta),b*np.cos(theta))
plt.savefig(uu)
plt.show()
        


