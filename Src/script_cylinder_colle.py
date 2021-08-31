
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import matplotlib as mpl
import math as mt
import csv
import sys
from scipy.integrate import tplquad
from mpl_toolkits.mplot3d import Axes3D


uu=sys.argv[1]
with open(uu, "w") as f:
    f.write("")

a=float(input("cylinder radius in nm?"))
l=float(input("length of the cylinder in nm?"))
nombre1=int(input("Number of dipoles per radius?"))
nombre2=int(input("Number of dipoles on the length ?"))
#nombre3=int(input("Number of tube ?"))
nombrefinal=(nombre1*nombre2)*3#*nombre3
cc=float(l/(nombre2))

phi=np.zeros(nombrefinal)
theta=np.zeros(nombrefinal)
psi=np.zeros(nombrefinal)
x=np.zeros(nombrefinal)
y=np.zeros(nombrefinal)
z=np.zeros(nombrefinal)


#***************************************
postube=0
kk=0
#for kui in range(nombre3):
for i in range(nombre1):
    for j in range(nombre2):
        phi[kk]=i*2*np.pi/nombre1
        theta[kk]=np.pi/2
        psi[kk]=0
        x[kk]=a*np.cos(phi[kk])
        y[kk]=a*np.sin(phi[kk])
        z[kk]=-(l/2)+cc*j
        kk=kk+1

for i in range(nombre1):
    for j in range(nombre2):
        phi[kk]=i*2*np.pi/nombre1
        theta[kk]=np.pi/2
        psi[kk]=0
        x[kk]=5+a*np.cos(phi[kk])
        y[kk]=5+a*np.sin(phi[kk])
        z[kk]=-(l/2)+cc*j
        kk=kk+1

for i in range(nombre1):
    for j in range(nombre2):
        phi[kk]=i*2*np.pi/nombre1
        theta[kk]=np.pi/2
        psi[kk]=0
        x[kk]=10+a*np.cos(phi[kk])
        y[kk]=10+a*np.sin(phi[kk])
        z[kk]=-(l/2)+cc*j
        kk=kk+1


xmax=np.amax(x)
ymax=np.amax(y)
zmax=np.amax(z)
b=1#(xmax+ymax+zmax)/10 

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
