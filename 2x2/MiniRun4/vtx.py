import numpy as np
import h5py
import h5flow
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.colors import LogNorm
#from matplotlib.colors import Normalize

x_boundaries = np.array([-63.931, -3.069, 3.069, 63.931])
y_boundaries = np.array([-329.8543, -207.8543])
z_boundaries = np.array([1236.3163,  -2.6837, 2.6837, 1364.3163])

center = np.array([0, -268, 1300])

def isIn(vertex):
  x, y, z, t = vertex
  if (x_boundaries[0] <= x <= x_boundaries[3]) and (y_boundaries[0] <= y <= y_boundaries[1]) and (z_boundaries[0] <= z <= z_boundaries[3]):
    return True
  else:
    return False

x,y,z=[],[],[]
xIn,yIn,zIn=[],[],[]

run = "MiniRun4"
Nfiles = 1

for d in range(Nfiles):
  #fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5"
  fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun4/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW/MiniRun4_1E19_RHC.flow.{d:05d}.FLOW.h5"
  if not os.path.exists(fname): continue                ### check if file exists. if not, continue
  print(fname)
  if d%10==0: print(d*100./Nfiles, " percent")
  
  file = h5py.File(fname, 'r')

  print(file.keys())
  print(file['mc_truth'].keys())
  print(file['charge'].keys())
 
  Ints = file["mc_truth/interactions/data"]
  stack = file["mc_truth/stack/data"]
  trajs = file["mc_truth/trajectories/data"]
  segments = file["mc_truth/segments/data"]
  light = file["mc_truth/light/data"]
  #charge = file["charge/data"]
  f = h5flow.data.H5FlowDataManager(fname, 'r')
 
  """  
  print(Ints.dtype.names)
  print(stack.dtype.names)
  print(trajs.dtype.names)
  print(segments.dtype.names)
  """

  #print("10", Ints['vertex'][0:10])

  vtx = Ints['vertex']

  for vt in vtx:

    x.append(vt[0])
    y.append(vt[1])
    z.append(vt[2])

    if isIn(vt)==True:
      xIn.append(vt[0])
      yIn.append(vt[1])
      zIn.append(vt[2])


  #print(vtx, len(vtx))
print(len(x), len(xIn))
print(len(y), len(yIn))
print(len(z), len(zIn))

print(xIn)
print(yIn)
print(zIn)

plt.figure(1)
plt.scatter(z, x)
plt.xlim(-500+center[2],500+center[2])
plt.ylim(-500+center[0],500+center[0])
rectangle = patches.Rectangle((z_boundaries[0], x_boundaries[0]), width=128, height=126, fill=False, edgecolor='red')
#rectangle = patches.Rectangle((-1236.3163, -63.931), width=128, height=126, fill=False, edgecolor='red')
plt.gca().add_patch(rectangle)
plt.xlabel("z")
plt.ylabel("x")
plt.grid()
plt.title("Interaction vertices on xz plane")
plt.savefig("xz.png")

plt.figure(2)
plt.scatter(z, y)
plt.xlim(-500+center[2],500+center[2])
plt.ylim(-500+center[1],500+center[1])
rectangle = patches.Rectangle((z_boundaries[0], y_boundaries[0]), width=128, height=122, fill=False, edgecolor='red')
plt.gca().add_patch(rectangle)
plt.xlabel("z")
plt.ylabel("y")
plt.grid()
plt.title("Interaction vertices on yz plane")
plt.savefig("yz.png")

plt.figure(3)
plt.scatter(zIn, xIn)
plt.xlim(-100+center[2],100+center[2])
plt.ylim(-100+center[0],100+center[0])
rectangle = patches.Rectangle((z_boundaries[0], x_boundaries[0]), width=128, height=126, fill=False, edgecolor='red')
plt.gca().add_patch(rectangle)
plt.xlabel("z")
plt.ylabel("x")
plt.grid()
plt.title("Interaction vertices on xz plane")
plt.savefig("xzIn.png")

plt.figure(4)
plt.scatter(zIn, yIn)
plt.xlim(-100+center[2],100+center[2])
plt.ylim(-100+center[1],100+center[1])
#rectangle = patches.Rectangle((-64.3163, -19.8543), width=128, height=122, fill=False, edgecolor='red', linestyle='--', label='Rectangle')
rectangle = patches.Rectangle((z_boundaries[0], y_boundaries[0]), width=128, height=122, fill=False, edgecolor='red')
plt.gca().add_patch(rectangle)
plt.xlabel("z")
plt.ylabel("y")
plt.grid()
plt.title("Interaction vertices on yz plane")
plt.savefig("yzIn.png")
