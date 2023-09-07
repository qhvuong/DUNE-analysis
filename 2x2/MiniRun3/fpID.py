import numpy as np
import h5py
import h5flow
import glob
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

ids = {111: "$\\pi^0$", 211: "$\\pi^+$", -211: "$\\pi^-$", 311: "$K^0$", 321: "$K^+$", -321: "$K^-$", 2212: "$p$", 2112: "$n$", 3122: "$\\Lambda$", 3222: "$\\Sigma^+$", 3212: "$\\Sigma^0$", 3112: "$\\Sigma^-$", 22: "$\\gamma$"}

print(ids.keys())

CCints = {}
NCints = {}

count = 0

#larnd_sim location
#file_paths = glob.glob('/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.larnd_v2/MiniRun3_1E19_RHC.larnd_v2.01008.LARNDSIM.h5')

#ndlar_flow location
file_paths = glob.glob("/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.000*.FLOW.h5")

f_count = 0

for fname in file_paths:
  f_count += 1
  if f_count%10==0: print(f_count*100.0/1022.0, "percent")

  file = h5py.File(fname, 'r')

  Ints = file["mc_truth/interactions/data"]
  stack = file["mc_truth/stack/data"]
  trajs = file["mc_truth/trajectories/data"]
  tracks = file["mc_truth/tracks/data"]
  f = h5flow.data.H5FlowDataManager(fname, 'r')
  
  Enu_CC.extend(Ints[Ints['isCC']==1]['Enu']/1E3)
  Enu_NC.extend(Ints[Ints['isCC']==0]['Enu']/1E3)

  #CC and NC vertices
  vtxCC.extend(Ints[Ints['isCC']==1]['vertex'])
  vtxNC.extend(Ints[Ints['isCC']==0]['vertex'])

  #print(len(Ints[Ints['isCC']==1]['vertexID']))
  #print(len(Ints[Ints['isCC']==0]['vertexID']))

  for i in Ints['vertexID']:
    if abs(Ints[Ints['vertexID']==i]['vertex'][0][0])>63.931 or abs(Ints[Ints['vertexID']==i]['vertex'][0][1]-42)>61.8543 or abs(Ints[Ints['vertexID']==i]['vertex'][0][2])>64.3163: continue
    count += 1
    pids = list(set(stack[(stack['vertexID']==i) & (stack['part_status']==1)]['part_pdg']))
    #print(stack[(stack['vertexID']==i) & (stack['part_status']==1)]['part_pdg']) 
    #print(Ints[Ints['vertexID']==i]['isCC'], pids)   
    for pid in ids.keys():
      name = ids[pid]
      #print(pid, name)
      if pid in pids:
        if pid==321:
          print(stack[(stack['vertexID']==i) & (stack['part_status']==1)]['part_pdg']) 
          print(Ints[Ints['vertexID']==i]['isCC'], pids)   
      #if abs(pid)==11 or abs(pid)==12 or abs(pid)==13 or abs(pid)==14 or abs(pid)==15 or abs(pid)==16: continue
        if Ints[Ints['vertexID']==i]['isCC']==1:
          CCints[name] = CCints.get(name,0) + 1
        elif Ints[Ints['vertexID']==i]['isCC']==0:
          NCints[name] = NCints.get(name,0) + 1
      else:
        if Ints[Ints['vertexID']==i]['isCC']==1:
          CCints[name] = CCints.get(name,0) + 0
        elif Ints[Ints['vertexID']==i]['isCC']==0:
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

