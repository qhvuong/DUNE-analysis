import numpy as np
import h5py
import h5flow
from h5flow.data import dereference
import matplotlib.pyplot as plt
import plotly.graph_objects as go

x_boundaries = np.array([-63.931, -3.069, 3.069, 63.931])
y_boundaries = np.array([-19.8543, 103.8543])
z_boundaries = np.array([-64.3163,  -2.6837, 2.6837, 64.3163])

def draw_cathode_planes(x_boundaries, y_boundaries, z_boundaries, **kwargs):
    
    traces = []
    for i_z in range(int(len(z_boundaries)/2)):
        for i_x in range(int(len(x_boundaries)/2)):
            z, y = np.meshgrid(np.linspace(z_boundaries[i_z * 2], z_boundaries[i_z * 2 + 1], 2), np.linspace(y_boundaries.min(), y_boundaries.max(),2))
            x = (x_boundaries[i_x * 2] + x_boundaries[i_x * 2 + 1]) * 0.5 * np.ones(z.shape)

            traces.append(go.Surface(x=x, y=y, z=z, **kwargs))

    return traces

def draw_anode_planes(x_boundaries, y_boundaries, z_boundaries, **kwargs):
    
    traces = []
    for i_z in range(int(len(z_boundaries)/2)):
        for i_x in range(int(len(x_boundaries))):           
            z, y = np.meshgrid(np.linspace(z_boundaries[i_z * 2], z_boundaries[i_z * 2 + 1], 2), np.linspace(y_boundaries.min(), y_boundaries.max(),2))
            x = x_boundaries[i_x] * np.ones(z.shape)

            traces.append(go.Surface(x=x, y=y, z=z, **kwargs))
    
    return traces

def vtx_in_LArActive(vtx, x_boundaries, y_boundaries, z_boundaries):
    vtx_in_Active = False
    if vtx[0] >= np.min(x_boundaries) and vtx[0] <= np.max(x_boundaries):
        if vtx[1] >= np.min(y_boundaries) and vtx[1] <= np.max(y_boundaries):
            if vtx[2] >= np.min(z_boundaries) and vtx[2] <= np.max(z_boundaries):
                vtx_in_Active = True
    
    return vtx_in_Active

def plot_segs(segs, **kwargs):
    
    def to_list(axis):
        return np.column_stack([
            segs[f'{axis}_start'],
            segs[f'{axis}_end'],
            np.full(len(segs), None)
        ]).flatten().tolist()
        
    x, y, z = (to_list(axis) for axis in 'xyz')

    trace = go.Scatter3d(x=x, y=y, z=z, **kwargs)
    
    return trace



fname = '/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.00008.FLOW.h5'

f = h5flow.data.H5FlowDataManager(fname, 'r')
file = h5py.File(fname, 'r')

evs = [71,]
if len(evs)!=0:
  for i_evt in evs:

    Ints_ev = f["charge/raw_events", "mc_truth/interactions", i_evt]['eventID'][0,0]	#associated eventID 
    
    PromptHits_ev = f["charge/events", "charge/calib_prompt_hits", i_evt]
    
    Segs_PromptHits = f["charge/raw_events","charge/calib_prompt_hits","charge/packets", "mc_truth/tracks", i_evt]
    
    Ints_PromptHits = f["charge/raw_events","charge/calib_prompt_hits","charge/packets", "mc_truth/tracks", "mc_truth/interactions", i_evt]
    print(Ints_PromptHits.shape)
    
    Ints_PHits_ev = Ints_PromptHits.data[0,:,0,:,0]
    
    PHits_Segs = f["charge/raw_events", "mc_truth/interactions", "mc_truth/trajectories", "mc_truth/tracks", "charge/packets", "charge/calib_prompt_hits", i_evt]
    print(PHits_Segs.shape)
    
    #Ints_PHits_ev = PHits_Segs.data[0,0,:,:,:,0]
    
    
    fig = go.Figure()
    
    plot_outActive_vertices = False
    
    ##########################
    # Draw the cathodes
    ##########################
    fig.add_traces(draw_cathode_planes(
    x_boundaries, y_boundaries, z_boundaries, 
    showscale=False,
    opacity=0.3,
    colorscale='Greys',
    ))
    
    ##########################
    # Draw the anodes
    ##########################
    fig.add_traces(draw_anode_planes(
    x_boundaries, y_boundaries, z_boundaries, 
    showscale=False,
    opacity=0.1,
    colorscale='ice',
    ))
    
    ##########################
    # Draw the prompt hits
    ##########################
    PHits_traces = go.Scatter3d(
    x=PromptHits_ev.data['x'].flatten(), y=PromptHits_ev.data['y'].flatten(), z=PromptHits_ev.data['z'].flatten(),
    marker_color=Segs_PromptHits[0,:,0,0]['trackID'],
    name='prompt hits',
    mode='markers',
    visible='legendonly',
    marker_size=3,
    marker_symbol='square',
    showlegend=True,
    opacity=0.7,
    customdata=Segs_PromptHits[0,:,0,0]['pdgId'].data,
    hovertemplate='<b>x:%{x:.3f}</b><br>y:%{y:.3f}<br>z:%{z:.3f}<br>pdg:%{customdata:d}',
    )
    fig.add_traces(PHits_traces)
    
    ##########################
    # Draw the vertices that are in the active LAr
    ##########################
    # Many hits will trace back to same vertices. Let's get rid some repetition here.
    vtxID_uniq_ev, vtx_idx = np.unique(Ints_PHits_ev['vertexID'], return_index=True)
    unique_Ints = np.take(Ints_PHits_ev,vtx_idx)
    Vtx_traces = []
    for i_vtx in range(len(unique_Ints)):
      if plot_outActive_vertices or vtx_in_LArActive(unique_Ints[i_vtx]['vertex'], x_boundaries, y_boundaries, z_boundaries):
        Vtx_traces.append(go.Scatter3d(
	x=[unique_Ints[i_vtx]['vertex'][0]], y=[unique_Ints[i_vtx]['vertex'][1]], z=[unique_Ints[i_vtx]['vertex'][2]],
	marker_color='red',
	name=f'vertex{i_vtx}',
	mode='markers',
	visible='legendonly',
	marker_size=5,
	marker_symbol='circle',
	showlegend=True,
	opacity=0.7
	))
    fig.add_traces(Vtx_traces)
    if len(Vtx_traces) == 0:
      print("No vertices in the active LAr volume!")

    ##########################
    # Draw associated segments
    ##########################
    Segs_traces = plot_segs(Segs_PromptHits[0,:,0,0], 
    mode="lines",
    name="edep segments",
    visible='legendonly',
    line_color='red',
    showlegend=True
    )
    fig.add_traces(Segs_traces)
    
    fig.update_layout(
    width=1024, height=768,
    legend_orientation="h",
    scene = dict(xaxis_title='x [cm]',
    yaxis_title='y [cm]',
    zaxis_title='z [cm]')    
    )
    
    fig.write_html(f"/dune/app/users/qvuong/data/kaon/{Ints_ev}.html")




