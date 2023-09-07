import numpy as np
import h5py
#import h5flow
import glob
import matplotlib.pyplot as plt
import os

### 2x2 boundaries
x_boundaries = np.array([-63.931, -3.069, 3.069, 63.931])
y_boundaries = np.array([-19.8543, 103.8543])
z_boundaries = np.array([-64.3163,  -2.6837, 2.6837, 64.3163])

def isInside(vertex):
  x, y, z = vertex
  if (x_boundaries[0] <= x <= x_boundaries[3]) and (y_boundaries[0] <= y <= y_boundaries[1]) and (z_boundaries[0] <= z <= z_boundaries[3]):
    return True
  else:
    return False

Enu_CC, Enu_NC = [], []
 
### lists of kaons produced in the 2x2
Enu_KpCC, Elep_KpCC, EKpCC = [], [], []
Enu_KpNC, Elep_KpNC, EKpNC = [], [], []

Kpid = 321

#larnd_sim location
#file_paths = glob.glob('/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.larnd_v2/MiniRun3_1E19_RHC.larnd_v2.01008.LARNDSIM.h5')

#ndlar_flow location
#file_paths = glob.glob("/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.*.FLOW.h5")
#fnames = [f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5" for d in range(1024)]

fout_KpCC = open("KpCCEvs.txt", "w")
fout_KpNC = open("KpNCEvs.txt", "w")
fout = open("KpEvs.txt", "w")

f_count = 0
Nfiles = 1024

for d in range(Nfiles):
  fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5"
  if not os.path.exists(fname): continue 		### check if file exists. if not, continue
  f_count += 1
  if f_count%10==0: print(f_count*100.0/Nfiles, "percent")

  file = h5py.File(fname, 'r')
  Ints = file["mc_truth/interactions/data"]		### genie interactions
  stack = file["mc_truth/stack/data"]			### genie stack

  ### Get the K+ vertices
  Kp_stack = stack[(stack['part_pdg']==Kpid) & (stack['part_status']==1)]['vertexID']

  for i in Kp_stack:
    if abs(Ints[Ints['vertexID']==i]['vertex'][0][0])>63.931 or abs(Ints[Ints['vertexID']==i]['vertex'][0][1]-42)>61.8543 or abs(Ints[Ints['vertexID']==i]['vertex'][0][2])>64.3163: continue
    if Ints[Ints['vertexID']==i]['target']!=18: continue
  
    fout.write(f"{d}, {i}\n")
    if Ints[Ints['vertexID']==i]['isCC']==1: 
      print("CC", fname, i)
      fout_KpCC.write(f"{d}, {i}\n")
    elif Ints[Ints['vertexID']==i]['isCC']==0: 
      print("NC", fname, i)
      fout_KpNC.write(f"{d}, {i}\n")

  print("\n")

fout_KpCC.close()
fout_KpNC.close()
fout.close()
