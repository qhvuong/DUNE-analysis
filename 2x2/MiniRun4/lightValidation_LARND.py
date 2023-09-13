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
z_boundaries = np.array([1236.3163,  -2.6837, 2.6837, 1364.3163])

center = np.array([0, -268, 1300])

def isIn(vertex):
  x, y, z, t = vertex
  if (x_boundaries[0] <= x <= x_boundaries[3]) and (y_boundaries[0] <= y <= y_boundaries[1]) and (z_boundaries[0] <= z <= z_boundaries[3]):
    return True
  else:
    return False

run = "MiniRun4"
Nfiles = 10

for n in range(Nfiles):
  fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun4/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.larnd/LARNDSIM/MiniRun4_1E19_RHC.larnd.{n:05d}.LARNDSIM.h5"
  #fname = f"/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun4/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW/MiniRun4_1E19_RHC.flow.{n:05d}.FLOW.h5"
  if not os.path.exists(fname): continue                ### check if file exists. if not, continue
  print(fname)
  if n%10==0: print(n*100./Nfiles, " percent")
  
  file = h5py.File(fname, 'r')

   
  print(file.keys())
  for key in file.keys():
    if key=='_header': 
      print(key, file[key].attrs.items(), "\n")
    else:
      print(key, file[key].dtype.names, "\n")
  

  #Ints = file["mc_truth/interactions/data"]
  #stack = file["mc_truth/stack/data"]
  #trajs = file["mc_truth/trajectories/data"]
  #segments = file["mc_truth/segments/data"]
  
  mc_hdr = file['mc_hdr']
  mc_stack = file['mc_stack']

  mc_packets_assn = file['mc_packets_assn']
  messages = file['messages']
  packets = file['packets']
  segments = file['segments']
  trajectories = file['trajectories']
  vertices = file['vertices']

  light_dat = file['light_dat']
  light_trig = file['light_trig']
  light_wvfm = file["light_wvfm"]

  f = h5flow.data.H5FlowDataManager(fname, 'r')

  print(mc_hdr['event_id'].shape)
  print(light_trig['ts_s'].shape)
  print(light_wvfm.shape)

"""
  ## DEFINE COMMON VALUES
  SPILL_PERIOD = 1.2e7
  SAMPLES = len(light_wvfm[0][0])
  BIT = min(x for x in abs(light_wvfm[0][0]) if x != 0)


  ## LOAD NECESSARY DATASETS
  packet_type = np.array(packets['packet_type'])
  p_tstamp = np.array(packets['timestamp'])
  l_tsync = np.array(light_trig['ts_sync'])
  spillID = np.array(segments['event_id'])
  segmentID = np.array(segments['segment_id'])
  opt_chan = np.array(light_trig['op_channel'])
  io_group = np.array(packets['io_group'])
  io_channel = np.array(packets['io_channel'])  


  ## INSPECT PACKET TYPES
  print('There are '+str(len(l_tsync))+' light events\n')
  df = pd.DataFrame(packet_type)
  df1 = df.value_counts()
  print('Packet Type\nof Charge Trigger:')
  print('------------------')
  print('trig.|     \ntype |count')
  print('-----------')
  print(df1)


  ## PACKET TYPES 4 AND 6 CAN MESS UP THIS CORRECTION, SO EXCLUDE THEM
  tstamp_trig0 = p_tstamp[packet_type==0]
  tstamp_trig7 = p_tstamp[packet_type==7]


  ## IDENTIFY THE INDEX WHERE THE TURNOVER OCCURS
  charge_cutoff = np.where(tstamp_trig0 > 1.999**31)[0][-1]
  light_cutoff = np.where(tstamp_trig7 > 1.999**31)[0][-1]
  wvfm_cutoff = np.where(l_tsync > 1.999**31)[0][-1]
  
  
  ## ADD 2^31 TO ALL TIMESTAMPS FOLLOWING THE TURNOVER
  tstamp_real_trig0 = np.concatenate((tstamp_trig0[:(charge_cutoff+1)],((2**31)+tstamp_trig0[(charge_cutoff+1):])))
  tstamp_real_trig7 = np.concatenate((tstamp_trig7[:(light_cutoff+1)],((2**31)+tstamp_trig7[(light_cutoff+1):])))
  l_tsync_real = np.concatenate((l_tsync[:(wvfm_cutoff+1)],((2**31)+l_tsync[(wvfm_cutoff+1):])))
  
  
  ## DEFINE SPILLID (EVENTID) FOR PACKETS AND LIGHT
  light_spillIDs = (np.rint(l_tsync_real/SPILL_PERIOD)).astype(int)
  packet0_spillIDs = (np.rint(tstamp_real_trig0/SPILL_PERIOD)).astype(int)
  packet7_spillIDs = (np.rint(tstamp_real_trig7/SPILL_PERIOD)).astype(int)
  list_spillIDs = np.unique(light_spillIDs)
  
  ## CHECK LARPIX TRIGGER VS LIGHT TRIGGER  
  fig = plt.figure(figsize=(18,6))
  indices = np.arange(0,len(p_tstamp),1)
  indices_0 = indices[packet_type==0]
  indices_7 = indices[packet_type==7]
  plt.plot(tstamp_real_trig0,indices_0, "o", color='dodgerblue', label='larpix')
  plt.plot(tstamp_real_trig7,indices_7,".", color='tomato', label='light')
  plt.axvline(x=(2**31), label='LArPix Clock Rollover')
  plt.title('Larpix (Spill) Trigger vs. Light Trigger\n', fontsize=18)
  plt.xlabel(r'Timestamp [0.01$\mu$s]', fontsize=14)
  plt.ylabel('Packet Index', fontsize=16)
  #plt.xlim(175, 195)
  plt.legend(fontsize=16)
  plt.savefig("evIDvsIndex.png")


  ## TRIGGERS IN LIGHT AND CHARGE PER SPILL
  fig = plt.figure(figsize=(14,6))
  bins = np.linspace(min(packet7_spillIDs),max(packet7_spillIDs),392)
  bin_width = bins[2] - bins[1]
  counts, bins = np.histogram(np.array(light_spillIDs), bins=bins)
  plt.hist(bins[:-1], bins, weights=counts, color='tomato', label='Light: '+str(len(l_tsync))+' triggers')
  counts, bins = np.histogram(np.array(packet7_spillIDs), bins=bins)
  plt.hist(bins[:-1], bins, weights=counts, histtype="step", color='dodgerblue', label='Pacman: '+str(len(packet7_spillIDs))+' triggers')
  plt.title('Triggers Per Spill ('+str(len(list_spillIDs))+' Spills)\n', fontsize=16)
  plt.xlabel('Spill', fontsize=14)
  plt.ylabel('Triggers', fontsize=14)
  plt.ylim(0,max(counts)+2)
  plt.xlim(0,40)
  plt.grid(axis='y', color='0.85')
  plt.legend(loc='upper left', fontsize=14)
  plt.savefig("TrigsPerSpill.png")


  ## CHECK SHAPE AND MODULE
  print(np.shape(light_wvfm))
  print('opt_chan[0][0] = '+str(opt_chan[0][0])+': Ergo, light_wvfm[0] is readout from Module 1')
  
  
  ## PLOT SINGLE WAVEFORM
  fig = plt.figure(figsize=(10,4))
  plt.plot(np.linspace(0,SAMPLES-1,SAMPLES),light_wvfm[0][0]/BIT, label='Opt. Chan. 0')
  plt.title('Module 1, Event '+str(light_spillIDs[0])+', Optical Channel 1', fontsize=16)
  plt.xlabel(r'Time Sample [0.01 $\mu$s]', fontsize=14)
  plt.ylabel('SiPM Channel Output', fontsize=14)
  plt.savefig("wvfm0.png")


  ## SELECT A SPILL TO INSPECT
  SPILL = 0
  
  
  ## ASSIGN "SUM CHANNEL" POSITIONS (this would be one side of one TPC)
  SiPM_struct = np.array([0,0,0,0,0,0,
  1,1,1,1,1,1,
  2,2,2,2,2,2,
  3,3,3,3,3,3])
  
  ## SELECT DATASETS BELONGING TO YOUR SPILL
  spill_light = np.where(light_spillIDs == SPILL)[0]
  
  ## CREATE EMPTY DATASETS FOR EACH LIGHT ARRAY (one side of one TPC)
  l_mod1_1L = np.zeros((24,SAMPLES))
  l_mod1_1R = np.zeros((24,SAMPLES))
  l_mod1_2L = np.zeros((24,SAMPLES))
  l_mod1_2R = np.zeros((24,SAMPLES))
  
  l_mod2_3L = np.zeros((24,SAMPLES))
  l_mod2_3R = np.zeros((24,SAMPLES))
  l_mod2_4L = np.zeros((24,SAMPLES))
  l_mod2_4R = np.zeros((24,SAMPLES)) 
  
  l_mod3_5L = np.zeros((24,SAMPLES))
  l_mod3_5R = np.zeros((24,SAMPLES))
  l_mod3_6L = np.zeros((24,SAMPLES))
  l_mod3_6R = np.zeros((24,SAMPLES))
  
  l_mod4_7L = np.zeros((24,SAMPLES))
  l_mod4_7R = np.zeros((24,SAMPLES))
  l_mod4_8L = np.zeros((24,SAMPLES))
  l_mod4_8R = np.zeros((24,SAMPLES)) 
  
  ## SORT THE LIGHT DATA BY MODULE, TPC, and SIDE
  for j in spill_light:
    if (opt_chan[j][0]) == 0: 
      l_mod1_1L = np.add(l_mod1_1L,light_wvfm[j][0:24])
      l_mod1_1R = np.add(l_mod1_1R,light_wvfm[j][24:48])
      l_mod1_2R = np.add(l_mod1_2R,light_wvfm[j][48:72])
      l_mod1_2L = np.add(l_mod1_2L,light_wvfm[j][72:96])

    if opt_chan[j][0]==96:
      l_mod2_3L = np.add(l_mod2_3L,light_wvfm[j][0:24])
      l_mod2_3R = np.add(l_mod2_3R,light_wvfm[j][24:48])
      l_mod2_4R = np.add(l_mod2_4R,light_wvfm[j][48:72])
      l_mod2_4L = np.add(l_mod2_4L,light_wvfm[j][72:96])

    if opt_chan[j][0]==192:
      l_mod3_5L = np.add(l_mod3_5L,np.array(light_wvfm[j][0:24]))
      l_mod3_5R = np.add(l_mod3_5R,np.array(light_wvfm[j][24:48]))
      l_mod3_6R = np.add(l_mod3_6R,np.array(light_wvfm[j][48:72]))
      l_mod3_6L = np.add(l_mod3_6L,np.array(light_wvfm[j][72:96])) 

    if opt_chan[j][0] == 288:
      l_mod4_7L = np.add(l_mod4_7L,np.array(light_wvfm[j][0:24]))
      l_mod4_7R = np.add(l_mod4_7R,np.array(light_wvfm[j][24:48]))
      l_mod4_8R = np.add(l_mod4_8R,np.array(light_wvfm[j][48:72]))
      l_mod4_8L = np.add(l_mod4_8L,np.array(light_wvfm[j][72:96])) 


  def data_readout(io_first, io_second, spill):
  ## SET UP AN 18-PLOT DISPLAY    
    fig = plt.figure(figsize=(13.8,8),tight_layout=True)
    subfigs = fig.subfigures(1, 6, wspace=0, width_ratios=[0.8,1.5,0.8,0.8,1.5,0.8], height_ratios=[1])
    axs0 = subfigs[0].subplots(4, 1,sharey=True,gridspec_kw={'hspace': 0})
    axs1 = subfigs[1].subplots(1, 1)
    axs2 = subfigs[2].subplots(4, 1,sharey=True,gridspec_kw={'hspace': 0})
    axs3 = subfigs[3].subplots(4, 1,sharey=True,gridspec_kw={'hspace': 0})
    axs4 = subfigs[4].subplots(1, 1)
    axs5 = subfigs[5].subplots(4, 1,sharey=True,gridspec_kw={'hspace': 0})

  ## CREATE AN EMPTY ARRAY TO AVOID RE-PLOTTING TRACKS
    plotted_tracks = []

  ## SET UP LABELING AND COLOR SCHEME
    titles = ["mod. 2, io_group 3","mod. 1, io_group 1","mod. 2, io_group 4","mod. 1, io_group 2",
              "mod. 4, io_group 7","mod. 3, io_group 5","mod. 4, io_group 8","mod. 3, io_group 6"]
    colors = ['aqua','aqua','lightgreen','lightgreen','yellow','yellow','orangered','orangered']
    cmap = cm.jet
    
  ## ENFORCE GEOMETRY
    ios = [3,1,4,2,7,5,8,6]
    left_data = [l_mod2_3L,l_mod1_1L,l_mod2_4L,l_mod1_2L,l_mod4_7L,l_mod3_5L,l_mod4_8L,l_mod3_6L]
    right_data = [l_mod2_3R,l_mod1_1R,l_mod2_4R,l_mod1_2R,l_mod4_7R,l_mod3_5R,l_mod4_8R,l_mod3_6R]

  ## ENSURE THE TIMESTAMP TURNOVER ISN'T AN ISSUE
    packet_list = packets[packet_type==0][packet0_spillIDs==spill]
    mc_assoc = mc_packets_assn[packet_type==0][packet0_spillIDs==spill]

  ## MAP PACKETS TO TRACKS
    for ip,packet in enumerate(packet_list):
      #print(packet_list)
      print("ip", ip)
      print("packet", packet)
      track_ids = mc_assoc['track_ids'][ip]
      io_group = packet['io_group']
      print("track_ids", track_ids)
      print("io_group", io_group)
  
  ## GET THE POSITION OF CHARGE TRACKS AND SAVE TO THE CORRECT IO_GROUP
    for trackid in track_ids:
      print(trackid)
      if trackid >= 0 and trackid not in plotted_tracks:
        print(segments[trackid])
        plotted_tracks.append(trackid)
        if io_group==io_first:
          X = (segments[trackid]['x_start']*10,segments[trackid]['x_end']*10)
          Y = (segments[trackid]['y_start']*10,segments[trackid]['y_end']*10)
          Z = (segments[trackid]['z_start']*10,segments[trackid]['z_end']*10)
          print('1', Y, Z)
          axs1.plot(Z,Y,c=colors[ios.index(io_first)],alpha=1,lw=1.5)
        if io_group==io_second:
          X = (segments[trackid]['x_start']*10,segments[trackid]['x_end']*10)
          Y = (segments[trackid]['y_start']*10,segments[trackid]['y_end']*10)
          Z = (segments[trackid]['z_start']*10,segments[trackid]['z_end']*10)
          print('4', Y, Z)
          axs4.plot(Z,Y,c=colors[ios.index(io_second)],alpha=1,lw=1.5)
        else:
          print("else")
          #axs1.plot(0,0,c='navy',alpha=0.1)
          #axs4.plot(0,0,c='navy',alpha=0.1)
          pass

  ## LABEL THE LIGHT PLOTS                            
    axs0[0].set_title("Left:\nio_group "+str(io_first))
    axs2[0].set_title("Right:\nio_group "+str(io_first))
    axs3[0].set_title("Left:\nio_group "+str(io_second))
    axs5[0].set_title("Right:\nio_group "+str(io_second))
    axs0[3].set_xlabel(r"Samples [0.01 $\mu$s]")
    axs2[3].set_xlabel(r"Samples [0.01 $\mu$s]")
    axs3[3].set_xlabel(r"Samples [0.01 $\mu$s]")
    axs5[3].set_xlabel(r"Samples [0.01 $\mu$s]")
    fig.supylabel("Pulse Sum Over Light Collection Module",x=-0.07,y=0.53)

  ## SUM THE LIGHT DATA (IN PARTS)  
    all_sums=[]
    for i in range(4):
      if (i%2)==0:
        clr = 'greenyellow'
      else:
        clr = 'lightgreen'
      wvfm_scndL = [sum(w) for w in zip(*(left_data[ios.index(io_second)])[SiPM_struct==i]/BIT)]
      wvfm_scndR = [sum(w) for w in zip(*(right_data[ios.index(io_second)])[SiPM_struct==i]/BIT)]
      wvfm_frstL = [sum(w) for w in zip(*(left_data[ios.index(io_first)])[SiPM_struct==i]/BIT)]
      wvfm_frstR = [sum(w) for w in zip(*(right_data[ios.index(io_first)])[SiPM_struct==i]/BIT)]

      all_sums.extend(wvfm_scndL+wvfm_scndR+wvfm_frstL+wvfm_frstR)
        
  ## SET UNIVERSAL AXIS LIMITS
      y_min = (min(all_sums)-500)
      y_max = (max(all_sums))

  ## PLOT LIGHT WAVEFORMS
      axs0[i].plot(np.linspace(0,SAMPLES-1,SAMPLES),wvfm_frstL,color='k')
      axs0[i].set_facecolor(clr) 
      axs0[i].set_box_aspect(1)
      axs0[i].label_outer()
      axs0[i].set_ylim(y_min,y_max)
      
      axs2[i].plot(np.linspace(0,SAMPLES-1,SAMPLES),wvfm_frstR,color='k')
      axs2[i].set_facecolor(clr)
      axs2[i].label_outer()
      axs2[i].set_box_aspect(1)
      axs2[i].set_ylim(y_min,y_max)
      axs2[i].yaxis.set_ticklabels([])
      
      axs3[i].plot(np.linspace(0,SAMPLES-1,SAMPLES),wvfm_scndL,color='k')
      axs3[i].set_facecolor(clr)
      axs3[i].label_outer()
      axs3[i].set_box_aspect(1)
      axs3[i].set_ylim(y_min,y_max)
      axs3[i].yaxis.set_ticklabels([])
      
      axs5[i].plot(np.linspace(0,SAMPLES-1,SAMPLES),wvfm_scndR,color='k')
      axs5[i].set_facecolor(clr)
      axs5[i].label_outer()
      axs5[i].set_box_aspect(1)
      axs5[i].set_ylim(y_min,y_max)
      axs5[i].yaxis.set_ticklabels([])

  ## COLOR THE CHARGE PLOTS
    #tpc_rectL = plt.Rectangle((12350,-3350), 650, 1300, linewidth=0.75, edgecolor='b', facecolor=cmap(0),zorder=-1)
    #tpc_rectR = plt.Rectangle((13000,-3350), 650, 1300, linewidth=0.75, edgecolor='b', facecolor=cmap(0),zorder=-1)
    #tpc_rectL = plt.Rectangle(((12665 - 640/2), (-2680 - 1340/2)), 640, 1340, linewidth=0.75, edgecolor='b', facecolor=cmap(0),zorder=-1)
    #tpc_rectR = plt.Rectangle(((12665 + 640/2), (-2680 - 1340/2)), 640, 1340, linewidth=0.75, edgecolor='b', facecolor=cmap(0),zorder=-1)
    tpc_rectL = plt.Rectangle((1234.5,-333), 64, 130, linewidth=0.75, edgecolor='b', facecolor=cmap(0),zorder=-1)
    tpc_rectR = plt.Rectangle((1301.5,-333), 64, 130, linewidth=0.75, edgecolor='b', facecolor=cmap(0),zorder=-1)

  ## LABEL THE CHARGE PLOTS    
    axs1.add_patch(tpc_rectL)
    axs1.set_aspect("equal")
    axs1.set_xlabel("z [cm]")
    axs1.set_ylim(-334,-202)
    axs1.set_xlim(1233.5, 1299.5)
    axs1.set_title(titles[ios.index(io_first)])
    axs1.yaxis.set_ticklabels([]) 

    axs4.add_patch(tpc_rectR) 
    axs4.set_xlabel("z [cm]")
    axs4.set_ylim(-334,-202)
    axs4.set_xlim(1300.5, 1366.5)
    axs4.set_aspect("equal")
    axs4.set_title(titles[ios.index(io_second)])
    axs4.yaxis.set_ticklabels([]) 

  ## SAVE FIG
    fig.savefig(f"LightVsCharge_{io_first}{io_second}{spill}.png")


  data_readout(3,1,SPILL)
  #data_readout(4,2,SPILL)
  #data_readout(7,5,SPILL)
  #data_readout(8,6,SPILL)

"""
