# Validation script of truth information for popular hadrons
import h5flow
import numpy as np
import matplotlib.pyplot as plt
# Open file
flow_out=h5flow.data.H5FlowDataManager("/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/MiniRun3_1E19_RHC.flow_v6.00000.FLOW.h5","r")
trajectories=flow_out["mc_truth/trajectories/data"]
segments=flow_out["mc_truth/tracks/data"]
# Variables for segment studies
length=[]
dEdx=[]
# Loop over segments
for j in segments:
    # Only include hadrons in the fiducial volume
    if (abs(j["x_start"])>62 or abs(j["y_start"]-42)>62 or abs(j["z_start"])>62): continue
    if (abs(j["pdgId"])!=211 and abs(j["pdgId"])!=2212 and abs(j["pdgId"])!=321): continue
    # Calculate segment length
    diffX=abs(j["x_end"]-j["x_start"])
    diffY=abs(j["y_end"]-j["y_start"])
    diffZ=abs(j["z_end"]-j["z_start"])
    # Also save dEdx
    dEdx.append(j["dEdx"])
    dist=np.sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ)
    length.append(dist)
# plot dE/dx from segments
plt.hist(dEdx,bins=100,range=[0,10])
plt.title(r"dE/dx for $\pi^{+/-}$, p, $K^{+/-}$ in the Active LAr")
plt.ylabel("Number of Segments")
plt.xlabel("dEdx [MeV/cm]")
plt.savefig("dEdxSegmentsHadrons.png")
plt.close()
# Set variables we are going to save for the trajectory study
energyFromTrackless=[]
lengthFromTrackless=[]
lengthSegFromTraj=[]
energyLeft=[]
xPointsFromTrackless=[]
zPointsFromTrackless=[]
xFromTrackless=[]
zFromTrackless=[]
xFromTracklessMicroBooNE=[]
zFromTracklessMicroBooNE=[]
a=0
# Loop over all trajectories
while a<len(trajectories):
    # Fill the trajectory information
    index=a
    j=flow_out["mc_truth/trajectories/data"][index]
    
    a=a+1
    # Only accept hadrons in the volume
    if (abs(j["xyz_start"][0])>62 or abs(j["xyz_start"][1]-42)>62 or abs(j["xyz_start"][2])>62): continue
    if (abs(j["pdgId"])!=211 and abs(j["pdgId"])!=2212 and abs(j["pdgId"])!=321): continue
    startingEnergy=j["E_start"]
    mass=938.27
    if (abs(j["pdgId"])==321):
        mass=493.677
    if (abs(j["pdgId"])==211):
        mass=139.57
    startingEnergy=startingEnergy-mass
    # Get segments from trajectory
    segmentsFromTraj=flow_out["mc_truth/trajectories","mc_truth/tracks",index][0]
    indicesFromTraj=np.where(flow_out["mc_truth/tracks/data"]["trackID"]==j["trackID"])[0]
    # Print info on trajectories with no segments
    if (len(indicesFromTraj)==0):
        print("TRAJECTORY FOUND WITH NO SEGMENTS")
        print("Number of flow objects associated: ",len(segmentsFromTraj), "Was it masked: ",segmentsFromTraj)
        print("Start Postion: ",j["xyz_start"][0],j["xyz_start"][1],j["xyz_start"][2])
        print("Particle ID ",j["trackID"],index)
        print("Parent ID: ",j["parentID"]," Parent PDG: ",flow_out["mc_truth/trajectories/data"][j["parentID"]]["pdgId"])
        daughter=np.where(flow_out["mc_truth/trajectories/data"]["parentID"]==j["trackID"])
        print("Daughter indices ",daughter)
        energy=j["E_start"]
        print("Energy of segmentless particle [MeV]: ",energy)
        print("PDG code of segmentless particle: ",j["pdgId"])
        print("End Process ",j["end_process"]," ",j["end_subprocess"])
        #print("Length of segmentless ",np.sqrt(dx*dx+dy*dy+dz*dz))
        #print("") 
        dx=j["xyz_start"][0]-j["xyz_end"][0]
        dy=j["xyz_start"][1]-j["xyz_end"][1]
        dz=j["xyz_start"][2]-j["xyz_end"][2]
        print("Length of segmentless ",np.sqrt(dx*dx+dy*dy+dz*dz))
        print("")
        if (np.sqrt(dx*dx+dy*dy+dz*dz)>0.3):
            xFromTrackless.append(j["xyz_start"][0])
            zFromTrackless.append(j["xyz_start"][2])
            xPointsFromTrackless.append([j["xyz_start"][0],j["xyz_end"][0]])
            zPointsFromTrackless.append([j["xyz_start"][2],j["xyz_end"][2]])
        if (np.sqrt(dx*dx+dy*dy+dz*dz)>0.03):
            xFromTracklessMicroBooNE.append(j["xyz_start"][0])
            zFromTracklessMicroBooNE.append(j["xyz_start"][2])
        lengthFromTrackless.append(np.sqrt(dx*dx+dy*dy+dz*dz))
        energyFromTrackless.append(energy)
        continue
    
    for i in segmentsFromTraj: #calculate energy loss after each segment
        remainingEnergy=startingEnergy-i["dE"]
        if (abs(i["x_start"])>62 or abs(i["y_start"]-42)>62 or abs(i["z_start"])>62): continue
   

        diffX=abs(i["x_end"]-i["x_start"])
        diffY=abs(i["y_end"]-i["y_start"])
        diffZ=abs(i["z_end"]-i["z_start"])

        dist=np.sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ)
        lengthSegFromTraj.append(dist) # save segment length and the remaining energy
        energyLeft.append(remainingEnergy)
# Print segment length vs energy left 
plt.scatter(lengthSegFromTraj,energyLeft)
plt.ylabel("Kinetic Energy Remaining After Segment [MeV]")
plt.title(r"Segment lengths for $\pi^{+/-}$, p, $K^{+/-}$ in the Active LAr")
plt.xlabel("Segment length [cm]")
plt.savefig("eRemainingPerSegment.png")

plt.xlim(0,1)
plt.savefig("eRemainingPerSegmentZoom.png")

plt.close()
print(lengthSegFromTraj)
print(energyLeft)
lengthSegFromTraj=np.ma.getdata(lengthSegFromTraj)
energyLeft=np.ma.getdata(energyLeft)
print(len(energyLeft),len(lengthSegFromTraj))
plt.hist2d(lengthSegFromTraj,energyLeft,bins=200,range=[[0,1],[0,500]])
plt.ylabel("Kinetic Energy Remaining After Segment [MeV]")
plt.title(r"Segment lengths for $\pi^{+/-}$, p, $K^{+/-}$ in the Active LAr")
plt.xlabel("Segment length [cm]")
plt.savefig("eRemainingPerSegment2DHist.png")
plt.close()

plt.scatter(xFromTrackless,zFromTrackless)
plt.ylabel("Z [cm]")
plt.xlabel("X [cm]")
plt.xlim(-63.93,63.93)
plt.ylim(-63.93,63.93)
plt.axhline(-1.4)
plt.axhline(1.4)
plt.axvline(-3.069)
plt.axvline(3.069)
plt.axvline(-33.5)
plt.axvline(33.5)
plt.title("Location of Hadron Trajectories>0.3 cm witout Segments")
plt.savefig("scatteringMissingSegments.png")
plt.close()



plt.scatter(xFromTracklessMicroBooNE,zFromTracklessMicroBooNE)
plt.xlim(-63.93,63.93)
plt.ylim(-63.93,63.93)
plt.axhline(-1.4)
plt.axhline(1.4)
plt.axvline(-3.069)
plt.axvline(3.069)
plt.axvline(-33.5)
plt.axvline(33.5)
plt.ylabel("Z [cm]")
plt.xlabel("X [cm]")
plt.title("Location of Hadron Trajectories>0.03 cm witout Segments")
plt.savefig("scatteringMissingSegments03.png")
plt.close()

i=0
while i<len(xPointsFromTrackless):
    x=xPointsFromTrackless[i]
    z=zPointsFromTrackless[i]
    plt.plot(x,z,"-")
    i=i+1
plt.xlim(-63.93,63.93)
plt.ylim(-63.93,63.93)
plt.axhline(-1.4,color="k")
plt.axhline(1.4,color="k")
plt.axvline(-3.069,color="k")
plt.axvline(3.069,color="k")
plt.axvline(-33.5,color="k")
plt.axvline(33.5,color="k")
plt.ylabel("Z [cm]")
plt.xlabel("X [cm]")
plt.title("Path of Hadron Trajectories>0.3 cm witout Segments")
plt.savefig("trajectoriesWithMissingSegments.png")
plt.close()



# Plot the trajectory length and energy for trajectories without any segments
plt.scatter(energyFromTrackless,lengthFromTrackless)
plt.ylabel("Length [cm]")
plt.xlabel("Energy [MeV]")
plt.title("Energy and Length of Trajectories with no Segments")
plt.plot()
plt.savefig("trajEnergyWithNoSegments.png")
plt.close()



# Plot segment length for segments hadrons in the fiducial volume

plt.hist(length, bins=100,range=[0,2])
plt.xlabel("Segment length [cm]")
plt.ylabel("Number of Segments")
plt.title(r"Segment lengths for $\pi^{+/-}$, p, $K^{+/-}$ in the Active LAr")
plt.savefig("segmentLengthInAr.png")
plt.close()


# I could not get this code to work, I do not have time to debug but the code above worked?
plt.hist(lengthSegFromTraj, bins=100,range=[0,2])
plt.xlabel("Segment length [cm]")
plt.ylabel("Number of Segments")
plt.title(r"Segment lengths from Trajs for $\pi^{+/-}$, p, $K^{+/-}$ in the Active LAr")
plt.savefig("segmentLengthFromTraj.png")
plt.close()

