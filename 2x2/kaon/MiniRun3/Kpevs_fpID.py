import numpy as np
import h5py
import h5flow
import glob
import os
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
#from matplotlib.colors import Normalize

x_boundaries = np.array([-63.931, -3.069, 3.069, 63.931])
y_boundaries = np.array([-19.8543, 103.8543])
z_boundaries = np.array([-64.3163,  -2.6837, 2.6837, 64.3163])

tgtCC, tgtNC = [], []
vtxCC, vtxNC = [], []
Enu_CC, Enu_NC = [], []
CCQES, CCMEC, CCRES, CCDIS, CCCOH = [], [], [], [], []
NCQES, NCMEC, NCRES, NCDIS, NCCOH = [], [], [], [], []
  
xCC, yCC, zCC = [], [], []
xNC, yNC, zNC = [], [], []

Enu_KpCC, Elep_KpCC, EKpCC = [], [], []
Enu_KpNC, Elep_KpNC, EKpNC = [], [], []

xKpCC, yKpCC, zKpCC = [], [], []
xKpNC, yKpNC, zKpNC = [], [], []

mKp = 493.677 #MeV
Kpid = 321
mupid = -13
Lambda0 = 3122
Sigmap = 3222
Sigma0 = 3212
Sigmam = 3112
Xi0 = 3322
Xim = 3312

ids = {111: "$\\pi^0$", 211: "$\\pi^+$", -211: "$\\pi^-$", 311: "$K^0$", 321: "$K^+$", -321: "$K^-$", 2212: "$p$", 2112: "$n$", 3122: "$\\Lambda$", 3222: "$\\Sigma^+$", 3212: "$\\Sigma^0$", 3112: "$\\Sigma^-$", 22: "$\\gamma$"}

strange_ids = [321, 3122, 3222, 3212, 3112, 3322, 3312]

CCints = {}
NCints = {}

run = "MiniRun3"
Nfiles = 10

for d in range(Nfiles):
  fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5"
  if not os.path.exists(fname): continue                ### check if file exists. if not, continue
    
  file = h5py.File(fname, 'r')
  
  Ints = file["mc_truth/interactions/data"]
  stack = file["mc_truth/stack/data"]
  trajs = file["mc_truth/trajectories/data"]
  tracks = file["mc_truth/tracks/data"]
  f = h5flow.data.H5FlowDataManager(fname, 'r')
  
  Kpvtx = stack[(stack['part_status']==1) & (stack['part_pdg']==Kpid)]['vertexID']

  for vtx in Kpvtx:
    if abs(Ints[Ints['vertexID']==vtx]['vertex'][0][0])>63.931 or abs(Ints[Ints['vertexID']==vtx]['vertex'][0][1]-42)>61.8543 or abs(Ints[Ints['vertexID']==vtx]['vertex'][0][2])>64.3163: continue
    if Ints[Ints['vertexID']==vtx]['target']!=18: continue

    Kp_trajID = trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['trackID']
    Kp_fspIDs = trajs[trajs['parentID']==Kp_trajID]['pdgId']

    #pids = list(set(stack[(stack['vertexID']==vtx) & (stack['part_status']==1)]['part_pdg']))
 
    #print(strange_ids)
    #print(pids)
 
    fid_strange = np.intersect1d(strange_ids, Kp_fspIDs)
    print(fid_strange)
    if len(fid_strange)<=1:
      print(stack.dtype.names)
      print(stack[(stack['vertexID']==vtx) & (stack['part_status']==1)]['part_pdg']) 
      print(trajs.dtype.names)
      print("primary trajectories", trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)])
      #print(Ints[Ints['vertexID']==vtx]['isCC'], pids)   

"""
    for pid in ids.keys():
      name = ids[pid]
      #print(pid, name)
      if pid in pids:
        if pid==321:
          print(stack[(stack['vertexID']==vtx) & (stack['part_status']==1)]['part_pdg']) 
          print(Ints[Ints['vertexID']==vtx]['isCC'], pids)   
      #if abs(pid)==11 or abs(pid)==12 or abs(pid)==13 or abs(pid)==14 or abs(pid)==15 or abs(pid)==16: continue
        if Ints[Ints['vertexID']==vtx]['isCC']==1:
          CCints[name] = CCints.get(name,0) + 1
        elif Ints[Ints['vertexID']==vtx]['isCC']==0:
          NCints[name] = NCints.get(name,0) + 1
      else:
        if Ints[Ints['vertexID']==vtx]['isCC']==1:
          CCints[name] = CCints.get(name,0) + 0
        elif Ints[Ints['vertexID']==vtx]['isCC']==0:
          NCints[name] = NCints.get(name,0) + 0
  
  #print(count)

keys = list(CCints.keys())
CCval = list(CCints.values())
NCval = list(NCints.values())

scale = 1E5/count

print(CCval, NCval)

#CCval = [num*scale for num in CCval]
#NCval = [num*scale for num in NCval]

print("final", count, scale)
#print(CCval, NCval)

plt.figure(1)
plt.bar(keys, CCval, label="CC")
plt.bar(keys, NCval, bottom=CCval, label="NC")
plt.ylabel(f"Number of interactions with this product particle per {count} Ints")
plt.legend()
plt.grid()
plt.yscale('log')
plt.savefig("pid.png")

"""
