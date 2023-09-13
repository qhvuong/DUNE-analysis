import numpy as np
import h5py
import pandas as pd
import h5flow
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

## 3D PLOTTING
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.gridspec as gridspec
from matplotlib import cm, colors
import matplotlib.patches as mpatches

#from matplotlib.colors import LogNorm
#from matplotlib.colors import Normalize

x_boundaries = np.array([-63.931, -3.069, 3.069, 63.931])
y_boundaries = np.array([-329.8543, -207.8543])
z_boundaries = np.array([1236.3163,  1298.6837, 1302.6837, 1364.3163])

center = np.array([0, -268, 1300])

def isIn(vertex):
  x, y, z = vertex
  if (x_boundaries[0] <= x <= x_boundaries[3]) and (y_boundaries[0] <= y <= y_boundaries[1]) and (z_boundaries[0] <= z <= z_boundaries[3]):
    return True
  else:
    return False

ep = 11
mup = 13
gamma = 22
M_mu = 105.6583745

Mc_E, Mc_p = [], []
Ion_E, brem_E, others_E = [], [], []

run = "MiniRun4"
Nfiles = 1025

for n in range(Nfiles):
  #fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun4/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.larnd/LARNDSIM/MiniRun4_1E19_RHC.larnd.{n:05d}.LARNDSIM.h5"
  fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun4/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW/MiniRun4_1E19_RHC.flow.{n:05d}.FLOW.h5"
  if not os.path.exists(fname): continue                ### check if file exists. if not, continue
  print(fname)
  if n%10==0: print(n*100./Nfiles, " percent")
  
  file = h5py.File(fname, 'r')

  Ints = file["mc_truth/interactions/data"]
  stack = file["mc_truth/stack/data"]
  trajs = file["mc_truth/trajectories/data"]
  segments = file["mc_truth/segments/data"]

  print(trajs.dtype.names)

  mu_trajIDs = trajs[abs(trajs['pdg_id'])==mup]['traj_id']
  print(len(mu_trajIDs))

  for t in mu_trajIDs:
    pids = trajs[trajs['parent_id']==t]['pdg_id']
    #print(pids)
    if ep not in pids: continue
    #if abs(trajs[trajs['traj_id']==t]['E_end']-M_mu)>1E-3: continue

    mup_xyz = np.array(trajs[(trajs['traj_id']==t)]['xyz_end'][0])
    if isIn(mup_xyz)==False: continue

    Mc_trajIDs = trajs[(trajs['parent_id']==t) & (abs(trajs['pdg_id'])==ep) & (trajs['start_process']==6)  & (trajs['start_subprocess']==201)]['traj_id']
    #Mc_trajIDs = trajs[(trajs['parent_id']==t) & (abs(trajs['pdg_id'])==ep) & (trajs['E_start']>10)]['traj_id']
    for Mc_t in Mc_trajIDs:
      print("Michel e", trajs[trajs['traj_id']==Mc_t], "\n")
      print("brem", trajs[(trajs['parent_id']==Mc_t) & (abs(trajs['pdg_id'])==gamma) & (trajs['start_process']==2) & (trajs['start_subprocess']==3)], "\n")
      print("ionization e", trajs[(trajs['parent_id']==Mc_t) & (abs(trajs['pdg_id'])==ep) & (trajs['start_process']==2) & (trajs['start_subprocess']==2)], "\n")
      print("others", trajs[(trajs['parent_id']==Mc_t) & (abs(trajs['pdg_id'])==ep)], "\n")
      Mc_E.extend(trajs[trajs['traj_id']==Mc_t]['E_start'])
      Ion_E.append(sum(trajs[(trajs['parent_id']==Mc_t) & (abs(trajs['pdg_id'])==ep) & (trajs['start_process']==2) & (trajs['start_subprocess']==2)]['E_start']))
      brem_E.append(sum(trajs[(trajs['parent_id']==Mc_t) & (abs(trajs['pdg_id'])==gamma) & (trajs['start_process']==2) & (trajs['start_subprocess']==3)]['E_start']))
      others_E.append(sum(trajs[(trajs['parent_id']==Mc_t) & (abs(trajs['pdg_id'])==ep) & (trajs['start_process']==2) & (trajs['start_subprocess']!=2)]['E_start']))
      #print(Mc_E, Ion_E)
      #print(trajs[trajs['traj_id']==Mc_t]['E_start'])
      #print(trajs[(trajs['parent_id']==Mc_t) & (abs(trajs['pdg_id'])==ep) & (trajs['start_process']==2)]['E_start'])

print(Mc_E, len(Mc_E))
print(Ion_E, len(Ion_E))
print(others_E, len(others_E))

plt.figure(1)
plt.hist2d(Mc_E, Ion_E, bins=(30,30), cmap='viridis_r', range=([0,60],[0,60]))
plt.title("Michel Energy vs Ionization Energy")
plt.xlabel("Michel electron Energy (MeV)")
plt.ylabel("Ionization Energy (MeV)")
plt.grid()
plt.colorbar()
plt.savefig("McVsIon_E.png")


plt.figure(4)
plt.hist2d(Mc_E, brem_E, bins=(30,30), cmap='viridis_r', range=([0,60],[0,60]))
plt.title("Michel Energy vs Bremsstrahlung Energy")
plt.xlabel("Michel electron Energy (MeV)")
plt.ylabel("Bremsstrahlung Energy (MeV)")
plt.grid()
plt.colorbar()
plt.savefig("McVsBrem_E.png")


plt.figure(2)
plt.hist2d(Mc_E, others_E, bins=(30,30), cmap='viridis_r', range=([0,60],[0,60]))
plt.title("Michel Energy vs Others (MeV)")
plt.xlabel("Michel electron Energy (MeV)")
plt.ylabel("Others")
plt.grid()
plt.colorbar()
plt.savefig("McVsOthers_E.png")


plt.figure(3)
plt.hist(Mc_E, bins=50)
plt.title("Michel Energy")
plt.xlabel("Michel electron Energy (MeV)")
plt.ylabel("counts")
plt.grid()
plt.savefig("Mc_E.png")
