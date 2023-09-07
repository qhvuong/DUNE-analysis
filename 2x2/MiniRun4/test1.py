import h5flow
import numpy as np
import matplotlib.pyplot as plt

flow_out=h5flow.data.H5FlowDataManager("MiniRun3D/MiniRun3_1E19_RHC.flow_v6.00000.FLOW.h5","r")
particleIndices54=np.where(flow_out["mc_truth/tracks/data"]["trackID"]==54)[0]
trajx=[]
trajy=[]
trajz=[]

for i in particleIndices54:
    trajx.append(flow_out["mc_truth/tracks/data"][i]["x_start"])
    trajy.append(flow_out["mc_truth/tracks/data"][i]["y_start"])
    trajz.append(flow_out["mc_truth/tracks/data"][i]["z_start"])

trajz=np.array(trajz)
trajx=np.array(trajx)
trajy=np.array(trajy)
print(trajx[0],trajz[0])
plt.plot(trajx[0],trajz[0],"-")
plt.savefig("trajStartParticle54.png")
plt.close()
hitsFromEvent0=flow_out["charge/events/","charge/calib_prompt_hits",0]
hitsFromParticle54=[]


for i in particleIndices54:
    hitsFromParticle54.append(flow_out["mc_truth/tracks","charge/packets","charge/calib_prompt_hits",i][0])
x=[]
y=[]
z=[]

for hits in hitsFromParticle54:
    x.append(hits["x"])
    y.append(hits["y"])
    z.append(hits["z"])
x=np.array(x)
y=np.array(y)
z=np.array(z)
print(x[0],z[0])
plt.plot(x[0],z[0],"o")
plt.savefig("hitsFromParticle54.png")
