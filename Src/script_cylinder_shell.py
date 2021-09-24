
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import matplotlib as mpl
import math as mt
import csv
import sys
import random
from scipy.integrate import tplquad
from mpl_toolkits.mplot3d import Axes3D


uu=sys.argv[1]
with open(uu, "w") as f:
    f.write("")

a=float(input("cylinder radius in nm?"))
l=float(input("length of the cylinder in nm?"))
l2=float(input("length in the directiobn in nm?"))
nombre1=int(input("Number of dipoles on inner radius?"))
nombre2=int(input("Number of dipoles on the length ?"))
nombre3=int(input("Number of dipoles on the direction ?"))
ind=float(input("indice de desordre, sigma distribution ?"))

ntotal=0
for op in range(1,nombre3+1):
    rayon=(l2)*(op)/nombre3
    Nrayon=int((rayon*nombre1)/a)
    print(Nrayon)
    ntotal=ntotal+Nrayon

nombrefinal=nombre2*ntotal
print("total number of dipoles on the cylinder:%d" % nombrefinal)
cc=float(l/(nombre2))

phi=np.zeros(nombrefinal)
theta=np.zeros(nombrefinal)
psi=np.zeros(nombrefinal)
x=np.zeros(nombrefinal)
y=np.zeros(nombrefinal)
z=np.zeros(nombrefinal)


#***************************************

kk=0
for j in range(nombre2):
    for op in range(1,nombre3+1): 
        rayon=(l2)*(op)/nombre3
        Nrayon=int((rayon*nombre1)/a)
        for i in range(Nrayon):
            phi[kk]=(np.random.randn(1)*ind+(i*2*np.pi/Nrayon))[0]
            theta[kk]=(np.random.randn(1)*ind+np.pi/2)[0]
            psi[kk]=0
            x[kk]=rayon*np.cos(i*2*np.pi/Nrayon)
            y[kk]=rayon*np.sin(i*2*np.pi/Nrayon)
            z[kk]=-(l/2)+cc*j
            kk=kk+1

xmax=np.amax(x)
ymax=np.amax(y)
zmax=np.amax(z)
b=(xmax+ymax+zmax)/10 

for j in range(0,nombrefinal):
    with open(uu, "a") as f:
        valeurs = "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" % (phi[j],theta[j],psi[j],x[j],y[j],z[j]) 
        f.write(valeurs)
        
print("total number of dipoles on the cylinder:%d" % nombrefinal)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.quiver(x,y,z,b*np.cos(phi)*np.sin(theta),b*np.sin(phi)*np.sin(theta),b*np.cos(theta))
plt.savefig(uu)
plt.show()
