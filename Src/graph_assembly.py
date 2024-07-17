import math, random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import csv
import sys


uuv=sys.argv[1]

nombre=int(input("number of dipoles?"))
data=[0]*nombre
value=[0]*nombre
dx=np.zeros(nombre)
dy=np.zeros(nombre)
dz=np.zeros(nombre)
x=np.zeros(nombre)
y=np.zeros(nombre)
z=np.zeros(nombre)

for i in range(0,nombre):
    with open('out_arrows_gnuplot','r') as f:
        data[i]=f.readlines()[i]
        value[i] = data[i].split()
        x[i]=float(value[i][0])
        y[i]=float(value[i][1])
        z[i]=float(value[i][2])
        dx[i]=float(value[i][3])
        dy[i]=float(value[i][4])
        dz[i]=float(value[i][5])

a_x = np.amax(x)
a_y= np.amax(y)
a_z= np.amax(z)
a_list=[a_x,a_y,a_z]
a=np.amax(a_list)
fig = plt.figure()
ax = fig.gca(projection='3d')
#ax._axis3don = False
ax.quiver(x,y,z,(a/3)*dx,(a/3)*dy,(a/3)*dz,color='b')
plt.savefig(uuv)
plt.show()
