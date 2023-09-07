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

for n in range(Nfiles):
  #fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5"
  fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun4/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW/MiniRun4_1E19_RHC.flow.{n:05d}.FLOW.h5"
  if not os.path.exists(fname): continue                ### check if file exists. if not, continue
  print(fname)
  if n%10==0: print(n*100./Nfiles, " percent")
  
  file = h5py.File(fname, 'r')

  """
  print(file.keys())
  for key in file.keys():
    print("\n", key, file[key].keys())
    for k,v in file[key].items():
      print(k, v)
  """
 
  Ints = file["mc_truth/interactions/data"]
  stack = file["mc_truth/stack/data"]
  trajs = file["mc_truth/trajectories/data"]
  segments = file["mc_truth/segments/data"]
  wvfm = file["light/wvfm/data"]
  sipm_hits = file["light/sipm_hits/data"]
  sum_hits = file["light/sum_hits/data"]
  #charge = file["charge/data"]
  f = h5flow.data.H5FlowDataManager(fname, 'r')

  #print(wvfm.dtype.names)
  #print(sipm_hits.dtype.names)
  #print(sum_hits.dtype.names)

  print(wvfm['samples'].shape)
  #print(wvfm['samples'][0][0][9])


  ev = 1
  TPCs = 8
  dummies, filled = 0, 0
  for tpc in range(TPCs):
    mod = tpc//2
    #tpc = 0
    print(tpc, mod)
    
    for i in range(64):
      SAMPLES = len(wvfm['samples'][ev][tpc][i])
      if not [x for x in abs(wvfm['samples'][ev][tpc][i]) if x != 0]:  
        dummies+=1
        continue

      filled+=1
      BIT = min(x for x in abs(wvfm['samples'][ev][tpc][i]) if x != 0)
      plt.clf()
      fig = plt.figure(figsize=(10,4))
      plt.plot(np.linspace(0,SAMPLES-1,SAMPLES),wvfm['samples'][ev][tpc][i]/BIT, label=f'Opt. Chan. {i}')
      plt.title(f'Event {ev}, Module {tpc//2}, Optical Channel {i}', fontsize=16)
      plt.xlabel(r'Time Sample [0.016 $\mu$s]', fontsize=14)
      plt.ylabel('SiPM Channel Output', fontsize=14)
      plt.savefig(f"wvfm{ev}{tpc}{i:02d}.png")
      plt.close()
    
    print(dummies, filled) 

  
