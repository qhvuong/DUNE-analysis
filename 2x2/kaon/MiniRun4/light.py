## STANDARD IMPORTS
#import uproot, h5py
import h5py
import numpy as np
import glob
import pandas as pd

## 3D PLOTTING
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.gridspec as gridspec
from matplotlib import cm, colors
import matplotlib.patches as mpatches

#larnd_sim location
fnames = glob.glob('/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.larnd_v2/MiniRun3_1E19_RHC.larnd_v2.00000.LARNDSIM.h5')

for fname in fnames:
  h5 = h5py.File(fname,'r')

  ## INSPECT FILE STRUCTURE
  
  print('\n-------------------File Contents--------------------')
  print('File Name: '+str(fname)+'\n')
  print(f'Available branches: {[t for t in h5.keys()]}')
  print('\n')

  header = h5['_header']
  print(f'Stored Values in _header:\n')
  for key, val in header.attrs.items():
    print(f'    %s: %s \n' % (key, val))
    
  configs = h5['configs']
  print(f'Stored Values in configs:\n')
  for key, val in configs.attrs.items():
    print(f'    %s: %s \n' % (key, val))
    
  print(f'Available keys in configs: {[t for t in configs.dtype.names]}')
  print('\n')

  genie_hdr = h5['genie_hdr']
  print(f'Available keys in genie_hdr: {[t for t in genie_hdr.dtype.names]}')
  print('\n')
  
  genie_stack = h5['genie_stack']
  print(f'Available keys in genie_stack: {[t for t in genie_stack.dtype.names]}')
  print('\n')

  mc_assn = h5['mc_packets_assn']
  print(f'Available keys in mc_packets_assn: {[t for t in mc_assn.dtype.names]}')
  print('\n')
  
  messages = h5['messages']
  print(f'Available keys in messages: {[t for t in messages.dtype.names]}')
  print('\n')
  
  packets = h5['packets']
  print(f'Available keys in packets: {[t for t in packets.dtype.names]}')
  print('\n')
  
  tracks = h5['tracks']
  print(f'Available keys in tracks: {[t for t in tracks.dtype.names]}')
  print('\n')
  
  traject = h5['trajectories']
  print(f'Available keys in trajectories: {[t for t in traject.dtype.names]}')
  print('\n')
  
  vertices = h5['vertices']
  print(f'Available keys in vertices: {[t for t in vertices.dtype.names]}')
  print('\n')
  
  light_dat = h5['light_dat']
  print(f'Available keys in light_dat: {[t for t in light_dat.dtype.names]}')
  print('\n')
  
  light_trig = h5['light_trig']
  print(f'Available keys in light_trig: {[t for t in light_trig.dtype.names]}')
  print('\n')
  
  light_wvfm = h5['light_wvfm']
  print(f'Stored Values in light_wvfm:\n')
  for key, val in light_wvfm.attrs.items():
    print(f'    %s: %s \n' % (key, val))
 
  print(f'Available keys in light_wvfm: NONE, it is a single nested dataset')
  print('\n')

  print(light_wvfm)
  #print(light_wvfm[0][0])
  #print(light_wvfm[1])

  ## DEFINE COMMON VALUES
  
  SPILL_PERIOD = 1.2e7
  
  SAMPLES = len(light_wvfm[0][0])
  
  BIT = min(x for x in abs(light_wvfm[0][0]) if x != 0)

  ## LOAD NECESSARY DATASETS
  
  pack_type = np.array(packets['packet_type'])
  p_tstamp = np.array(packets['timestamp'])
  l_tsync = np.array(light_trig['ts_sync'])
  spillID = np.array(tracks['eventID'])
  trackID = np.array(tracks['trackID'])
  opt_chan = np.array(light_trig['op_channel'])
  io_group = np.array(packets['io_group'])
  io_channel = np.array(packets['io_channel'])
  
  ## INSPECT PACKET TYPES
  
  print('There are '+str(len(l_tsync))+' light events\n')

  print(pack_type[0], type(packets['packet_type']))
 
  """ 
  df = pd.DataFrame(pack_type)
  df1 = df.value_counts()
  print('Packet Type\nof Charge Trigger:')
  print('------------------')
  print('trig.|     \ntype |count')
  print('-----------')
  print(df1)
  
  ## SEE THAT TIMESTAMP TURNOVER MUST BE ACCOUNTED FOR 
  ## WHEN LINKING DATA-LIKE AND TRUTH-LIKE BRANCHES
  
  print(max(spillID))
  print(max(p_tstamp)/SPILL_PERIOD)
  print(max(l_tsync)/SPILL_PERIOD)
  
  ## PACKET TYPES 4 AND 6 CAN MESS UP THIS CORRECTION, SO EXCLUDE THEM
  
  tstamp_trig0 = p_tstamp[pack_type==0]
  tstamp_trig7 = p_tstamp[pack_type==7]
  
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
  
  fig = plt.figure(figsize=(18,6))
  
  indices = np.arange(0,len(p_tstamp),1)
  indices_0 = indices[pack_type==0]
  indices_7 = indices[pack_type==7]
  
  plt.plot(tstamp_real_trig0,indices_0, "o", color='dodgerblue', label='larpix')
  plt.plot(tstamp_real_trig7,indices_7,".", color='tomato', label='light')
  plt.axvline(x=(2**31), label='LArPix Clock Rollover')
  
  plt.title('Larpix (Spill) Trigger vs. Light Trigger\n', fontsize=18)
  plt.xlabel(r'Timestamp [0.01$\mu$s]', fontsize=14)
  plt.ylabel('Packet Index', fontsize=16)
  #plt.xlim(175, 195)
  plt.legend(fontsize=16)
  plt.savefig("SpillVsLightTrigger.png")

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
  plt.savefig("TriggersPerSpill.png")
  
  ## CHECK SHAPE AND MODULE
  
  print(np.shape(light_wvfm))
  print('opt_chan[0][0] = '+str(opt_chan[0][0])+': Ergo, light_wvfm[0] is readout from Module 1')
  
  ## PLOT SINGLE WAVEFORM
  
  fig = plt.figure(figsize=(10,4))
  plt.plot(np.linspace(0,SAMPLES-1,SAMPLES),light_wvfm[0][0]/BIT, label='Opt. Chan. 0')
  plt.title('Module 1, Event '+str(light_spillIDs[0])+', Optical Channel 1', fontsize=16)
  plt.xlabel(r'Time Sample [0.01 $\mu$s]', fontsize=14)
  plt.ylabel('SiPM Channel Output', fontsize=14)
  plt.savefig("SingleWaveform.png")
  
  ## A PEAK
  
  fig = plt.figure(figsize=(10,4))
  plt.plot(np.linspace(0,SAMPLES-1,SAMPLES),light_wvfm[19][0]/BIT, label='Opt. Chan. 0')
  plt.title('Module 3, Timestamp '+str(light_spillIDs[19])+', Optical Channel 196', fontsize=16)
  plt.xlabel(r'Time Sample [0.01 $\mu$s]', fontsize=14)
  plt.ylabel('SiPM Channel Output', fontsize=14)
  plt.savefig("APeak.png")
  """

