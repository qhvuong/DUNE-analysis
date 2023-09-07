import numpy as np
import h5py
import h5flow
import glob
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.colors import LogNorm
#from matplotlib.colors import Normalize

x_boundaries = np.array([-63.931, -3.069, 3.069, 63.931])
y_boundaries = np.array([-19.8543, 103.8543])
z_boundaries = np.array([-64.3163,  -2.6837, 2.6837, 64.3163])

def isInside(vertex):
  x, y, z = vertex
  if (x_boundaries[0] <= x <= x_boundaries[3]) and (y_boundaries[0] <= y <= y_boundaries[1]) and (z_boundaries[0] <= z <= z_boundaries[3]):
    return True
  else:
    return False

mKp = 493.677 #MeV

Lambda0 = 3122
Sigmap = 3222
Sigma0 = 3212
Sigmam = 3112
Xi0 = 3322
Xim = 3312

pids = {321: "Kp", -13: "mup", 211: "pip", -211: "pim"}

list_dict = {}
for pid in pids:
  list_dict[pid] = {
    'dt_trajs': [],
    'dE': [],
    'dx': [],
    'dEdx': [],
    't0': []
  }

#ids = {111: "$\\pi^0$", 211: "$\\pi^+$", -211: "$\\pi^-$", 311: "$K^0$", 321: "$K^+$", -321: "$K^-$", 2212: "$p$", 2112: "$n$", 3122: "$\\Lambda$", 3222: "$\\Sigma^+$", 3212: "$\\Sigma^0$", 3112: "$\\Sigma^-$", 22: "$\\gamma$"}
#strange_ids = [321, 3122, 3222, 3212, 3112, 3322, 3312]

#larnd_sim location
#file_paths = glob.glob('/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.larnd_v2/MiniRun3_1E19_RHC.larnd_v2.01008.LARNDSIM.h5')

run = "MiniRun4"
Nfiles = 10

for d in range(Nfiles):
  #fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5"
  fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun4/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW/MiniRun4_1E19_RHC.flow.{d:05d}.FLOW.h5"
  if not os.path.exists(fname): continue                ### check if file exists. if not, continue
  print(fname)
  if d%10==0: print(d*100./Nfiles, " percent")
  
  file = h5py.File(fname, 'r')
  
  Ints = file["mc_truth/interactions/data"]
  stack = file["mc_truth/stack/data"]
  trajs = file["mc_truth/trajectories/data"]
  tracks = file["mc_truth/tracks/data"]
  f = h5flow.data.H5FlowDataManager(fname, 'r')

  for vtx in Ints['vertexID']:
    if abs(Ints[Ints['vertexID']==vtx]['vertex'][0][0])>63.931 or abs(Ints[Ints['vertexID']==vtx]['vertex'][0][1]-42)>61.8543 or abs(Ints[Ints['vertexID']==vtx]['vertex'][0][2])>64.3163: continue
    if Ints[Ints['vertexID']==vtx]['target']!=18: continue
    for pid in pids:
      if pid not in trajs[trajs['vertexID']==vtx]['pdgId']: continue
      list_dict[pid]['dt_trajs'].extend(trajs[(trajs['vertexID']==vtx) & (trajs['pdgId']==pid)]['t_end'] - trajs[(trajs['vertexID']==vtx) & (trajs['pdgId']==pid)]['t_start'])
      list_dict[pid]['dE'].extend(tracks[(tracks['vertexID']==vtx) & (tracks['pdgId']==pid)]['dE'])
      list_dict[pid]['dEdx'].extend(tracks[(tracks['vertexID']==vtx) & (tracks['pdgId']==pid)]['dEdx'])
      list_dict[pid]['dx'].extend(tracks[(tracks['vertexID']==vtx) & (tracks['pdgId']==pid)]['dx'])
      list_dict[pid]['t0'].extend(tracks[(tracks['vertexID']==vtx) & (tracks['pdgId']==pid)]['t0'])
      #list_dict[pid]['dt'].extend(tracks[(tracks['vertexID']==vtx) & (tracks['pdgId']==pid)]['t0_end'] - tracks[(tracks['vertexID']==vtx) & (tracks['pdgId']==pid)]['t0_start'])
  
  #print(len(list_dict[pid]['dt_trajs']), len(list_dict[pid]['dEdx']))


for pid in pids:
  t0 = [num - list_dict[pid]['t0'][0] for num in list_dict[pid]['t0']]

  if len(list_dict[pid]['dt_trajs'])==0: continue

  name = pids[pid]

  if pid==-13:
    r = (0,20)
  else:
    r = (0,8)

  plt.figure(pid+1)
  plt.hist(list_dict[pid]['dt_trajs'], bins=40, range=r)
  plt.title(f"t_end - t_start for {name} trajectories")
  plt.xlabel("dt (us)")
  plt.ylabel("counts")
  plt.grid()
  plt.savefig(f"{run}_{name}_dt.png")
  plt.close(pid+1) 
  """  
  plt.figure(pid+2)
  plt.hist2d(t0, list_dict[pid]['dEdx'], bins = (20,20), cmap='viridis_r', range=([0,max(t0)],[0,50]))
  #plt.hist2d(list_dict[pid]['t0'], list_dict[pid]['dEdx'], bins = (20,20), cmap='viridis_r', range=([min(list_dict[pid]['t0']),max(list_dict[pid]['t0'])],[0,100]))
  #plt.hist2d(list_dict[pid]['t0'], list_dict[pid]['dEdx'], bins = (20,20), cmap='viridis_r')
  plt.title(f"dEdx vs t0 for {name} segments")
  plt.xlabel("t0 (us)")
  plt.ylabel("dEdx (MeV/cm)")
  plt.grid()
  plt.colorbar()
  plt.savefig(f"{name}_dEdxVst0.png")
  plt.close(pid+2) 
  """
  plt.figure(pid+3)
  plt.hist2d(list_dict[pid]['dx'], list_dict[pid]['dE'], bins = (20,20), cmap='viridis_r', range=([0,15],[0,30]))
  #plt.hist2d(list_dict[pid]['dx'], list_dict[pid]['dE'], bins = (20,20), cmap='viridis_r')
  plt.title(f"dE vs dx for {name} segments")
  plt.xlabel("dx (cm)")
  plt.ylabel("dE (MeV)")
  plt.grid()
  plt.colorbar()
  plt.savefig(f"{run}_{name}_dEVsdx.png")
  plt.close(pid+3) 

"""
print(dvtxCCmup)
print(dvtxNCmup)

dxCCmup, dyCCmup, dzCCmup = [], [], []
dxNCmup, dyNCmup, dzNCmup = [], [], []
for vtx in dvtxCCmup:
  dxCCmup.append(vtx[0])
  dyCCmup.append(vtx[1])
  dzCCmup.append(vtx[2])
for vtx in dvtxNCmup:
  dxNCmup.append(vtx[0])
  dyNCmup.append(vtx[1])
  dzNCmup.append(vtx[2])

#print(dtKpCC)
nKpCC = len(dtKpCC)
nKpNC = len(dtKpNC)

nKpCCmup = len(dtKpCCmup)
nKpNCmup = len(dtKpNCmup)

nCCmup = len(dtCCmup)
nNCmup = len(dtNCmup)

plt.figure(1)
label_dtKp = [f"CC: {nKpCC}", f"NC: {nKpNC}"]
plt.hist([dtKpCC, dtKpNC], histtype='stepfilled', stacked=True, label=label_dtKp, bins=20)
plt.title("(tKp_end - tKp_start) for K+ contained in the 2x2")
plt.xlabel("ns")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("dtKp.png")

plt.figure(2)
label_dtKpmup = [f"CC: {nKpCCmup}", f"(NC: {nKpNCmup}"]
plt.hist([dtKpCCmup, dtKpNCmup], histtype='stepfilled', stacked=True, label=label_dtKpmup, bins=20)
plt.title("(tKp_end - tKp_start) for K+ contained in the 2x2 with mu+ as a product")
plt.xlabel("ns")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("dtKpmup.png")

plt.figure(3)
label_dtmup = [f"CC: {nCCmup}", f"NC: {nNCmup}"]
plt.hist([dtCCmup, dtNCmup], histtype='stepfilled', stacked=True, label=label_dtmup, bins=20)
plt.title("(tmu_start - tKp_end) for K+ contained in the 2x2 with mu+ as a product")
plt.xlabel("ns")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("dtmup.png")


plt.figure(4)
label_dxmup = [f"CC: {nCCmup}", f"NC: {nNCmup}"]
plt.hist([dxCCmup, dxNCmup], histtype='stepfilled', stacked=True, label=label_dxmup, bins=20)
plt.title("(xmu_start - xKp_end) for K+ contained in the 2x2 with mu+ as a product")
plt.xlabel("cm")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("dxmup.png")

plt.figure(5)
label_dymup = [f"CC: {nCCmup}", f"NC: {nNCmup}"]
plt.hist([dyCCmup, dyNCmup], histtype='stepfilled', stacked=True, label=label_dymup, bins=20)
plt.title("(ymu_start - yKp_end) for K+ contained in the 2x2 with mu+ as a product")
plt.xlabel("cm")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("dymup.png")

plt.figure(6)
label_dzmup = [f"CC: {nCCmup}", f"NC: {nNCmup}"]
plt.hist([dzCCmup, dzNCmup], histtype='stepfilled', stacked=True, label=label_dzmup, bins=20)
plt.title("(zmu_start - zKp_end) for K+ contained in the 2x2 with mu+ as a product")
plt.xlabel("cm")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("dzmup.png")
"""


