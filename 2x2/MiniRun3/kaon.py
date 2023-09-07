import numpy as np
import h5py
import h5flow
import glob
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

def plot_track(ax, E_start, vtx_start, t_start, E_end, vtx_end, t_end, track_id):
  # Unpack coordinates
  x_start, y_start, z_start = vtx_start
  x_end, y_end, z_end = vtx_end
  
  # Plot the track as a line
  ax.plot([x_start, x_end], [y_start, y_end], [z_start, z_end], label=track_id)
  
  # Annotate the start and end points with time values
  #ax.text(x_start, y_start, z_start, f't1={t_start}', color='red')
  #ax.text(x_end, y_end, z_end, f't2={t_end}', color='red')


Kpid = 321
mupid = -13

EKp_end, EKpmup_end = [], []

tKp_start, tKp_end = [], []
tKpmup_start, tKpmup_end = [], []
tmup_start = []

vtxKp_start, vtxKp_end = [], []
vtxKpmup_start, vtxKpmup_end = [], []

dtKp, dtKpmup, dtmup = [], [], []

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
  
  #CC and NC vertices
  vtxCC = Ints[Ints['isCC']==1]['vertexID']
  vtxNC = Ints[Ints['isCC']==0]['vertexID']

  #kaon vertices:
  vtxKp = stack[(stack['part_status']==1) & (stack['part_pdg']==Kpid)]['vertexID']

  vtxKpCC = np.intersect1d(vtxCC, vtxKp) 
  vtxKpNC = np.intersect1d(vtxNC, vtxKp) 

  for i in vtxKpNC:
    if abs(Ints[Ints['vertexID']==i]['vertex'][0][0])>63.931 or abs(Ints[Ints['vertexID']==i]['vertex'][0][1]-42)>61.8543 or abs(Ints[Ints['vertexID']==i]['vertex'][0][2])>64.3163: continue
    Kp_trajID = trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['trackID'] 
    if isInside(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_end'][0]) != True: continue
    print(fname)
    print(Ints[Ints['vertexID']==i])  
    print(stack[(stack['vertexID']==i) & (stack['part_status']==1) & (stack['part_pdg']==Kpid)])
    print(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)])
    print("\n\n")

    EKp_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end']/1E3)
    tKp_start.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_start']*1E3)
    tKp_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_end']*1E3)
 
    vtxKp_start.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_start'])
    vtxKp_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_end'])


    if mupid not in trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID)]['pdgId']: continue
    EKpmup_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end']/1E3)
    tKpmup_start.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_start']*1E3)
    tKpmup_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['t_end']*1E3)
    tmup_start.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==mupid)]['t_start']*1E3)
    
    vtxKpmup_start.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_start'])
    vtxKpmup_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['xyz_end'])


dtKp = [a-b for a,b in zip(tKp_end, tKp_start)]
dtKpmup = [a-b for a,b in zip(tKpmup_end, tKpmup_start)]
dtmup = [a-b for a,b in zip(tmup_start, tKpmup_start)]
    
print(EKp_end)
print(EKpmup_end, "\n")

print(tKp_start)
print(tKp_end, "\n")

print(vtxKp_start)
print(vtxKp_end, "\n")

print(tKpmup_start)
print(tKpmup_end)
print(tmup_start, "\n")

print(vtxKpmup_start)
print(vtxKpmup_end, "\n")



"""
    Kp, Kpmup = [], []
    keys = ['E_start', 'xyz_start', 't_start', 'E_end', 'xyz_end', 't_end', 'pdgId']
    for key in keys:
      Kp.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)][key])
    

    if mupid not in trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID)]['pdgId']: continue
    EKpmup_end.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)]['E_end'])
    for key in keys:
      Kpmup.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==-1) & (trajs['pdgId']==Kpid)][key])

     
    for j in trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID)]['trackID']:
      #print(j, trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID) & (trajs['pdgId']==j)])
      product = []
      for key in keys:
        product.extend(trajs[(trajs['vertexID']==i) & (trajs['parentID']==Kp_trajID) & (trajs['trackID']==j)][key])
      print(product)
      products.append(product)
    print("\n\n")
   

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    plot_track(ax, *Kp)
    for product in products:
      plot_track(ax, *product)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.savefig(f"{i}_Kp_track.png")
"""

