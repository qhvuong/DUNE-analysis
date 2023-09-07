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
M_mu = 105.6583745

Mc_E, Mc_p = [], []
ep_E = []

run = "MiniRun4"
Nfiles = 1 

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

  mu_trajIDs = trajs[trajs['pdg_id']==mup]['traj_id']
  print(len(mu_trajIDs))

  for t in mu_trajIDs:
    pids = trajs[trajs['parent_id']==t]['pdg_id']
    #print(pids)
    if ep not in pids: continue
    if abs(trajs[trajs['traj_id']==t]['E_end']-M_mu)>1E-3: continue

    mup_xyz = np.array(trajs[(trajs['traj_id']==t)]['xyz_end'][0])
    #print(trajs[(trajs['traj_id']==t) & (trajs['pdg_id']==mup)]['xyz_end'][0])
    if isIn(mup_xyz)==False: continue

    #Mc_tIDs = trajs[(trajs['parent_id']==t) & (trajs['pdg_id']==ep) & (trajs['start_process']!=2)]['traj_id']
    Mc_trajIDs = trajs[(trajs['parent_id']==t) & (trajs['pdg_id']==ep) & (trajs['E_start']>10)]['traj_id']
    print(trajs[trajs['traj_id']==Mc_trajIDs])

"""
    for Mc_t in Mc_trajIDs:
      ep_xyz = np.array(trajs[(trajs['traj_id']==Mc_t)]['xyz_start'][0])
  
      if not np.array_equal(ep_xyz, mup_xyz): continue
      print(trajs[trajs['traj_id']==t])
      print(trajs[trajs['traj_id']==Mc_t])
      Mc_E.extend(trajs[(trajs['traj_id']==Mc_t)]['E_start'])

      Mc_pxyz = trajs[(trajs['traj_id']==Mc_t)]['pxyz_start'][0]
      Mc_p.append(np.sqrt(Mc_pxyz[0]**2 + Mc_pxyz[1]**2 + Mc_pxyz[2]**2))

      ep_trajIDs = trajs[(trajs['parent_id']==Mc_t) & (trajs['pdg_id']==ep)]['traj_id'] 
      print(len(ep_trajIDs)) 
      print(trajs[trajs['traj_id']==ep_trajIDs])

  print(len(Mc_E), len(Mc_p))
"""

"""
plt.figure(1)
label_EMichel = [f"entries: {len(Mc_E)}"]
plt.hist(Mc_E, bins=50, label=label_EMichel, range=(0,50))
plt.grid()
plt.title("Michel Electron Energy")
plt.xlabel("MeV")
plt.ylabel("entries")
plt.savefig("Mc_E_noG4ID.png")


plt.figure(2)
label_pMichel = [f"entries: {len(Mc_p)}"]
plt.hist(Mc_p, bins=50, label=label_pMichel, range=(0,50))
plt.grid()
plt.title("Michel Electron Momentum")
plt.xlabel("MeV")
plt.ylabel("entries")
plt.savefig("Mc_p_noG4ID.png")
"""
