
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

#Reading number of cores---------------------------
f=open("AFSSH.inp", "r")
i=0

lines=f.readlines()
for line in lines:
    i=i+1
    line=line.strip().split()
    if (i==3):
        cores=int(line[0])
        print(cores)
f.close()
#Averaging-----------------------------------------

init=np.loadtxt(os.path.join("1",sys.argv[1]))
data=init[:,1]
#state=init[:,2]
for j in range(2,cores+1):
    raw=np.loadtxt(os.path.join(str(j),sys.argv[1]))
    data=data+raw[:,1]
    #state=state+raw[:,2]

data=data/cores
#state=state/cores

x=init[:,0]
y=data
#z=state
#--------------------------------------------------
data = np.column_stack((x, y))
np.savetxt('output.txt', data, fmt='%.6f', header='X Y', comments='')
#plt.plot(init[:,0],1-data)
#plt.show()

