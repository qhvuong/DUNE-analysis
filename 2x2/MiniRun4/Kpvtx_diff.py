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


EKpCC_start, EKpCC_end= [], []
EKpNC_start, EKpNC_end= [], []
EKpCCmup_start, EKpCCmup_end= [], []
EKpNCmup_start, EKpNCmup_end= [], []

dtKpCC, dtKpCCmup, dtCCmup  = [], [], []
dtKpNC, dtKpNCmup, dtNCmup  = [], [], []

dvtxCCmup, dvtxNCmup = [], []

#larnd_sim location
#file_paths = glob.glob('/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.larnd_v2/MiniRun3_1E19_RHC.larnd_v2.01008.LARNDSIM.h5')

lcCC, lcNC, lc = 0,0,0
with open("KpCCEvs.txt", "r") as fin:
  for line in fin:
    lcCC += 1
with open("KpNCEvs.txt", "r") as fin:
  for line in fin:
    lcNC += 1
with open("KpEvs.txt", "r") as fin:
  for line in fin:
    lc += 1
print(lcCC, lcNC, lc)

lcount = 0
with open("KpNCEvs.txt", "r") as fin:
  for line in fin:
    lcount+=1
    if lcount>10: break
    d, vtx = map(int, line.strip().split(", "))
    fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5"
    if not os.path.exists(fname): continue                ### check if file exists. if not, continue

    file = h5py.File(fname, 'r')
    
    Ints = file["mc_truth/interactions/data"]
    stack = file["mc_truth/stack/data"]
    trajs = file["mc_truth/trajectories/data"]
    tracks = file["mc_truth/tracks/data"]
    f = h5flow.data.H5FlowDataManager(fname, 'r')
    
    if abs(Ints[Ints['vertexID']==vtx]['vertex'][0][0])>63.931 or abs(Ints[Ints['vertexID']==vtx]['vertex'][0][1]-42)>61.8543 or abs(Ints[Ints['vertexID']==vtx]['vertex'][0][2])>64.3163: continue
    #if isInside(trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_end'][0]) != True: continue

    ### trackID of Kp
    Kp_trajID = trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['trackID'] 
    Kp_fspID = trajs[trajs['parentID']==Kp_trajID]['pdgId']
    
    var = ['x_start', 'y_start', 'z_start', 't_start', '']

    print(trajs.dtype.names)
    print("primary trajectories", trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)])
    print(tracks.dtype.names)
    print("tracks", tracks[(tracks['vertexID']==vtx) & (tracks['trackID']==Kp_trajID)]['dE'], "\n") 
    #print(Kp_trajID, "\n")

    i = 1
    while Kpid in Kp_fspID: 
      i+=1
      print(f"{i}th trajectories", trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==Kpid)])
      print("tracks", tracks[(tracks['vertexID']==vtx) & (tracks['trackID']==Kp_trajID)]['dE'], "\n") 
      #Kp_trajID = trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==Kpid)]['trackID'] 
      Kp_trajID = trajs[(trajs['parentID']==Kp_trajID) & (trajs['pdgId']==Kpid)]['trackID'] 
      Kp_fspID = trajs[trajs['parentID']==Kp_trajID]['pdgId']
      #print(Kp_trajID, "\n")
      #print(f"{i+2}th trajectories", trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==Kpid)], "\n")
      
"""
    fs_trajID = trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==Kp_trajID)]['trackID'] 
    KpIDs = trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==Kp_trajID)]['pdgId']


    if any(pid in KpIDs for pid in strange_ids): continue
    print(fname)
    print("stack", stack[(stack['vertexID']==vtx) & (stack['part_status']==1) & (stack['part_pdg']==Kpid)])
    print("trajectories", trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)], "\n")
    print(tracks.dtype)
    print("tracks", tracks[(tracks['vertexID']==vtx) & (tracks['trackID']==Kp_trajID)], "\n") 
    print("dE", tracks[(tracks['vertexID']==vtx) & (tracks['trackID']==Kp_trajID)]['dE'], "\n") 
    print("dx", tracks[(tracks['vertexID']==vtx) & (tracks['trackID']==Kp_trajID)]['dx'], "\n") 
    print("dEdx", tracks[(tracks['vertexID']==vtx) & (tracks['trackID']==Kp_trajID)]['dEdx'], "\n") 
    print(trajs[(trajs['vertexID']==vtx) & (trajs['parentID']==Kp_trajID)], "\n")
    print(tracks.dtype)
    for fs in fs_trajID:
      print(tracks[(tracks['vertexID']==vtx) & (tracks['trackID']==fs)], "\n\n") 
   
    #print(Ints[Ints['vertexID']==i])  
    #print(stack[(stack['vertexID']==i) & (stack['part_status']==1) & (stack['part_pdg']==Kpid)])
    #print(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)])

    #if Ints[Ints['vertexID']==i]['isCC']==1:
      #EKpCC_start.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_start'])
      #EKpCC_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end'])
      #dEKpCC.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_start']) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end']))
      #dtKpCC.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_end']*1E3) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_start']*1E3))
      #dvtxKpCC.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_end']) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_start']))

      #if mupid in trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID)]['pdgId']:
        #EKpCCmup_start.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_start'])
        #EKpCCmup_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end'])
        #dEKpCCmup.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_start']) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end']))
        #dtKpCCmup.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_end']*1E3) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_start']*1E3))
        #dtCCmup.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==mupid)]['t_start']*1E3) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_end']*1E3))
        #dvtxCCmup.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==mupid)]['xyz_start']) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_end']))


    #elif Ints[Ints['vertexID']==i]['isCC']==0:
      #EKpNC_start.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_start'])
      #EKpNC_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end'])
      #dEKpNC.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_start']) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end']))
      #dtKpNC.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_end']*1E3) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_start']*1E3))
      #dvtxKpNC.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_end']) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_start']))

      #if mupid in trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID)]['pdgId']:
        #dtKpNCmup.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_end']*1E3) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_start']*1E3))
        #dtNCmup.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==mupid)]['t_start']*1E3) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_end']*1E3))
        #dvtxNCmup.extend((trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==mupid)]['xyz_start']) - (trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_end']))
"""

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


