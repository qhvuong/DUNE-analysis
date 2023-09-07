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
y_boundaries = np.array([-19.8543, 103.8543])
z_boundaries = np.array([-64.3163,  -2.6837, 2.6837, 64.3163])

def isIn(vertex):
  x, y, z, t = vertex
  if (x_boundaries[0] <= x <= x_boundaries[3]) and (y_boundaries[0] <= y <= y_boundaries[1]) and (z_boundaries[0] <= z <= z_boundaries[3]):
    return True
  else:
    return False

x,y,z=[],[],[]
xIn,yIn,zIn=[],[],[]

run = "MiniRun3"
Nfiles = 1

for d in range(Nfiles):
  fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5"
  #fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun4/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW/MiniRun4_1E19_RHC.flow.{d:05d}.FLOW.h5"
  if not os.path.exists(fname): continue                ### check if file exists. if not, continue
  print(fname)
  if d%10==0: print(d*100./Nfiles, " percent")
  
  file = h5py.File(fname, 'r')

  print(file.keys())
  print(file['mc_truth'].keys())
 
  Ints = file["mc_truth/interactions/data"]
  stack = file["mc_truth/stack/data"]
  trajs = file["mc_truth/trajectories/data"]
  #segments = file["mc_truth/segments/data"]
  f = h5flow.data.H5FlowDataManager(fname, 'r')
 
  """
  print(Ints.dtype.names)
  print(stack.dtype.names)
  print(trajs.dtype.names)
  print(segments.dtype.names)
  """

  print("10", Ints['vertex'][0:10])

  vtx = Ints['vertex']

  for vt in vtx:
    #vt = vtx[i]

    x.append(vt[0])
    y.append(vt[1])
    z.append(vt[2])

    #print(vt)

    #print(vt)
    if isIn(vt)==True:
      xIn.append(vt[0])
      yIn.append(vt[1])
      zIn.append(vt[2])


  #print(vtx, len(vtx))
print(len(x), len(xIn))
print(len(y), len(yIn))
print(len(z), len(zIn))
"""
plt.figure(1)
plt.scatter(z, x)
plt.xlim(-500,500)
plt.ylim(-500,500)
rectangle = patches.Rectangle((-64.3163, -63.931), width=128, height=126, fill=False, edgecolor='red')
plt.gca().add_patch(rectangle)
plt.xlabel("z")
plt.ylabel("x")
plt.grid()
plt.title("Interaction vertices on xz plane")
plt.savefig("xz.png")

plt.figure(2)
plt.scatter(z, y)
plt.xlim(-500,500)
plt.ylim(-500,500)
#rectangle = patches.Rectangle((-64.3163, -19.8543), width=128, height=122, fill=False, edgecolor='red', linestyle='--', label='Rectangle')
rectangle = patches.Rectangle((-64.3163, -19.8543), width=128, height=122, fill=False, edgecolor='red')
plt.gca().add_patch(rectangle)
plt.xlabel("z")
plt.ylabel("y")
plt.grid()
plt.title("Interaction vertices on yz plane")
plt.savefig("yz.png")

plt.figure(3)
plt.scatter(zIn, xIn)
plt.xlim(-100,100)
plt.ylim(-100,100)
rectangle = patches.Rectangle((-64.3163, -63.931), width=128, height=126, fill=False, edgecolor='red')
plt.gca().add_patch(rectangle)
plt.xlabel("z")
plt.ylabel("x")
plt.grid()
plt.title("Interaction vertices on xz plane")
plt.savefig("xzIn.png")

plt.figure(4)
plt.scatter(zIn, yIn)
plt.xlim(-100,100)
plt.ylim(-50,150)
#rectangle = patches.Rectangle((-64.3163, -19.8543), width=128, height=122, fill=False, edgecolor='red', linestyle='--', label='Rectangle')
rectangle = patches.Rectangle((-64.3163, -19.8543), width=128, height=122, fill=False, edgecolor='red')
plt.gca().add_patch(rectangle)
plt.xlabel("z")
plt.ylabel("y")
plt.grid()
plt.title("Interaction vertices on yz plane")
plt.savefig("yzIn.png")
"""
"""
  for vtx in Ints['vertex_id']:
    print(Ints[Ints['vertex_id']==vtx]['vertex'])
    #print(trajs[(trajs['vertex_id']==vtx) & (trajs['pdg_id']==pid)])
    if abs(Ints[Ints['vertex_id']==vtx]['vertex'][0][0])>63.931 or abs(Ints[Ints['vertex_id']==vtx]['vertex'][0][1]-42)>61.8543 or abs(Ints[Ints['vertex_id']==vtx]['vertex'][0][2])>64.3163: continue
    #print(fname)
    if Ints[Ints['vertex_id']==vtx]['target']!=18: continue
    print("inside", Ints[Ints['vertex_id']==vtx]['vertex'])
"""

