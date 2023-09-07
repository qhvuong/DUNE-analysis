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
fnames = [f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5" for d in range(1024)]
f_count = 0

for fname in fnames:
  if not os.path.exists(fname): continue 		### check if file exists. if not, continue

  f_count += 1
  if f_count%10==0: print(f_count*100.0/len(fnames), "percent")

  file = h5py.File(fname, 'r')

  Ints = file["mc_truth/interactions/data"]		### genie interactions
  stack = file["mc_truth/stack/data"]			### genie stack
  trajs = file["mc_truth/trajectories/data"]		### GEANT4 trajectories
  tracks = file["mc_truth/tracks/data"]			### GEANT4 tracks

  ### show keys in different datasets:
  """
  print("keys in genie interactions: ", Ints.dtype.names, "\n")
  print("keys in genie stack: ", stack.dtype.names, "\n")
  print("keys in GEANT4 trajectories: ", trajs.dtype.names, "\n")
  print("keys in GEANT4 tracks: ", tracks.dtype.names, "\n")
  """
  
  ### Get the CC and NC vertices
  CCints = Ints[Ints['isCC']==1]['vertexID']
  NCints = Ints[Ints['isCC']==0]['vertexID']
  
  ### Get the K+ vertices
  Kp_stack = stack[(stack['part_pdg']==Kpid) & (stack['part_status']==1)]['vertexID']

  ### Find CC K+ and NC K+ vertices (since vertexID is globally unique and it can be used to match interactions)
  KpCC_vtx = np.intersect1d(CCints, Kp_stack)
  KpNC_vtx = np.intersect1d(NCints, Kp_stack)


  eKpNC, vtx_KpNC = [], []
  eKpCC, vtx_KpCC = [], []
  for i in Ints['vertexID']:
    ### Require event vertices to be within the 2x2 volume and the target is Aron:
    if abs(Ints[Ints['vertexID']==i]['vertex'][0][0])>63.931 or abs(Ints[Ints['vertexID']==i]['vertex'][0][1]-42)>61.8543 or abs(Ints[Ints['vertexID']==i]['vertex'][0][2])>64.3163: continue
    if Ints[Ints['vertexID']==i]['target']!=18: continue

    Enu_CC.extend(Ints[(Ints['vertexID']==i) & (Ints['isCC']==1)]['Enu']/1E3)
    Enu_NC.extend(Ints[(Ints['vertexID']==i) & (Ints['isCC']==0)]['Enu']/1E3)

    pids = list(set(stack[(stack['vertexID']==i) & (stack['part_status']==1)]['part_pdg']))	### list of pdg ids of particles produced from the neutrino-Ar interactions

    if Kpid in pids: 
      if Ints[Ints['vertexID']==i]['isCC']==1:
        Enu_KpCC.extend(Ints[Ints['vertexID']==i]['Enu']/1E3)
        Elep_KpCC.extend(Ints[Ints['vertexID']==i]['Elep']/1E3)
        eKpCC.append(stack[(stack['vertexID']==i) & (stack['part_status']==1) & (stack['part_pdg']==Kpid)]['part_4mom'][0])
      if Ints[Ints['vertexID']==i]['isCC']==0:
        Enu_KpNC.extend(Ints[Ints['vertexID']==i]['Enu']/1E3)
        Elep_KpNC.extend(Ints[Ints['vertexID']==i]['Elep']/1E3)
        eKpNC.append(stack[(stack['vertexID']==i) & (stack['part_status']==1) & (stack['part_pdg']==Kpid)]['part_4mom'][0])
    
      
  print(len(Enu_CC), len(Enu_NC), len(Enu_KpCC), len(Enu_KpNC))

  EKpCC.extend(arr[3]/1E3 for arr in eKpCC)
  EKpNC.extend(arr[3]/1E3 for arr in eKpNC)

nCC = len(Enu_CC)
nNC = len(Enu_NC)
nKpCC = len(EKpCC)
nKpNC = len(EKpNC)

nInts = nCC + nNC 	#total number of neutrino interactions within the active volume
nKp = nKpCC + nKpNC	#number of K+ events produced within the active volume


#scale=1E5/nInts		#scale to 1E5 neutrino interactions
scale = 1.0


nuCC = [x - y - z for x, y, z in zip(Enu_KpCC, Elep_KpCC, EKpCC)]
nuNC = [x - y - z for x, y, z in zip(Enu_KpNC, Elep_KpNC, EKpNC)]


plt.figure(1)
label_Enu = [f"CC: {nCC}", f"NC: {nNC}"]
#plt.hist([Enu_CC, Enu_NC], histtype='stepfilled', stacked=True, label=label_Enu, bins=50, range=(0,15), weights=[scale*np.ones_like(Enu_CC), scale*np.ones_like(Enu_NC)])
plt.hist([Enu_CC, Enu_NC], histtype='stepfilled', stacked=True, label=label_Enu, bins=50, range=(0,15))
plt.title("Neutrino energy (within the 2x2)")
plt.xlabel("GeV")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("Enu.png")


plt.figure(2)
label_EnuKp = [f"CC: {nKpCC}", f"NC: {nKpNC}"]
plt.hist([Enu_KpCC, Enu_KpNC], histtype='stepfilled', stacked=True, label=label_EnuKp, bins=50, range=(0,15))
plt.title("Neutrino energy for K+ events (within the 2x2)")
plt.xlabel("GeV")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("Enu_Kp.png")


plt.figure(3)
label_EKp = [f"CC: {nKpCC}", f"NC: {nKpNC}"]
plt.hist([EKpCC, EKpNC], histtype='stepfilled', stacked=True, label=label_EKp, bins=50, range=(0,5))
plt.title("K+ energy (within the 2x2)")
plt.xlabel("GeV")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("EKp.png")


plt.figure(4)
label_nu = [f"CC: {nKpCC}", f"NC: {nKpNC}"]
plt.hist([nuCC, nuNC], histtype='stepfilled', stacked=True, label=label_nu, bins=50, range=(0,5))
plt.title("Energy transfer Enu-Elep-EKp (within the 2x2)")
plt.xlabel("GeV")
plt.ylabel("counts")
plt.legend()
plt.grid()
plt.savefig("nu_Kp.png")


