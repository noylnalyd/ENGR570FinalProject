#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% Importing modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import os
import sys

#%% reading in arguments
nProc = int(sys.argv[1])
nSamples = int(sys.argv[2])
nParam = int(sys.argv[3])

#%% defining constants
aFile = "a.txt"
bFile = "b.txt"
cFile = "c.txt"
pFile = "p.txt"

fA = pd.DataFrame()
fB = pd.DataFrame()
fAB1 = pd.DataFrame()
fAB2 = pd.DataFrame()
fAB3 = pd.DataFrame()
fAB4 = pd.DataFrame()
fBA1 = pd.DataFrame()
fBA2 = pd.DataFrame()
fBA3 = pd.DataFrame()
fBA4 = pd.DataFrame()

#%% reading in input data

aIn = pd.read_csv(aFile, sep=",", names=["TsrmIndoors (K)","TsrmOutdoors (K)","TToRecovery (sec)","TMachine (K)"])
bIn = pd.read_csv(bFile, sep=",", names=["TsrmIndoors (K)","TsrmOutdoors (K)","TToRecovery (sec)","TMachine (K)"])
cIn = pd.read_csv(cFile, sep=",", names=["TsrmIndoors (K)","TsrmOutdoors (K)","TToRecovery (sec)","TMachine (K)"])
pIn = pd.read_csv(pFile, sep=",", names=["Start Index","Stop Index (inclusive)"])

#%% reading in output data

for proc in range(nProc):
    fAFile="runs/"+str(proc+1)+"/pID"+str(proc+1)+".A_out.txt"
    tempA=pd.read_csv(fAFile, names=["Time to Stable Temp"])
    fA=pd.concat([fA,tempA],ignore_index=True)
    
    fBFile="runs/"+str(proc+1)+"/pID"+str(proc+1)+".B_out.txt"
    tempB=pd.read_csv(fBFile, names=["Time to Stable Temp"])
    fB=pd.concat([fB,tempB],ignore_index=True)
    
    fAB1File="runs/"+str(proc+1)+"/pID"+str(proc+1)+".A_B_0_out.txt"
    tempAB1=pd.read_csv(fAB1File, names=["Time to Stable Temp"])
    fAB1=pd.concat([fAB1,tempAB1],ignore_index=True)
    
    fAB2File="runs/"+str(proc+1)+"/pID"+str(proc+1)+".A_B_1_out.txt"
    tempAB2=pd.read_csv(fAB2File, names=["Time to Stable Temp"])
    fAB2=pd.concat([fAB2,tempAB2],ignore_index=True)
    
    fAB3File="runs/"+str(proc+1)+"/pID"+str(proc+1)+".A_B_2_out.txt"
    tempAB3=pd.read_csv(fAB3File, names=["Time to Stable Temp"])
    fAB3=pd.concat([fAB3,tempAB3],ignore_index=True)
    
    fAB4File="runs/"+str(proc+1)+"/pID"+str(proc+1)+".A_B_3_out.txt"
    tempAB4=pd.read_csv(fAB4File, names=["Time to Stable Temp"])
    fAB4=pd.concat([fAB4,tempAB4],ignore_index=True)
    
    fBA1File="runs/"+str(proc+1)+"/pID"+str(proc+1)+".B_A_0_out.txt"
    tempBA1=pd.read_csv(fBA1File, names=["Time to Stable Temp"])
    fBA1=pd.concat([fBA1,tempBA1],ignore_index=True)
    
    fBA2File="runs/"+str(proc+1)+"/pID"+str(proc+1)+".B_A_1_out.txt"
    tempBA2=pd.read_csv(fBA2File, names=["Time to Stable Temp"])
    fBA2=pd.concat([fBA2,tempBA2],ignore_index=True)
    
    fBA3File="runs/"+str(proc+1)+"/pID"+str(proc+1)+".B_A_2_out.txt"
    tempBA3=pd.read_csv(fBA3File, names=["Time to Stable Temp"])
    fBA3=pd.concat([fBA3,tempBA3],ignore_index=True)
    
    fBA4File="runs/"+str(proc+1)+"/pID"+str(proc+1)+".B_A_3_out.txt"
    tempBA4=pd.read_csv(fBA4File, names=["Time to Stable Temp"])
    fBA4=pd.concat([fBA4,tempBA4],ignore_index=True)
    
    
fC=pd.concat([fA,fB],ignore_index=True)

#%% optionally plotting input spaces
plt.figure()
cIn["TsrmIndoors (K)"].hist(bins=20)
plt.title("Indoor Temperature Prior Samples")
plt.ylabel("Temp (K)")
plt.xlabel("Frequency")
plt.figure()
cIn["TsrmOutdoors (K)"].hist(bins=20)
plt.title("Outdoor Temperature Prior Samples")
plt.ylabel("Temp (K)")
plt.xlabel("Frequency")
plt.figure()
cIn["TToRecovery (sec)"].hist(bins=20)
plt.title("Ambulance Arrival Time Prior Samples")
plt.ylabel("Time (s)")
plt.xlabel("Frequency")
plt.figure()
cIn["TMachine (K)"].hist(bins=20)
plt.title("Temperature of Saline-Blood Mixture Prior Samples")
plt.ylabel("Temp (K)")
plt.xlabel("Frequency")

#%% plotting output distribution given parameter uncertainty

plt.figure()
fC["Time to Stable Temp"].hist(bins=20)
plt.title("Solver Outputs: Time to Stable Temperature")
plt.ylabel("Time (s)")
plt.xlabel("Frequency")

#%% Calculating stuff of interest

meanC = fC["Time to Stable Temp"].mean()
sT=np.zeros((nParam,))
s=np.zeros((nParam,))
for i in range(nParam):
    num=0; numT=0; denom=0
    if i==0:
        fAB = fAB1; fBA=fBA1
    elif i==1:
        fAB = fAB2; fBA=fBA2
    elif i==2:
        fAB = fAB3; fBA=fBA3
    else:
        fAB = fAB4; fBA=fBA4
    for j in range(nSamples):
        num=num+(fA.loc[j]["Time to Stable Temp"]*fBA.loc[j]["Time to Stable Temp"] \
                 - fA.loc[j]["Time to Stable Temp"]*fB.loc[j]["Time to Stable Temp"])
        numT=numT+(fA.loc[j]["Time to Stable Temp"]-fAB.loc[j]["Time to Stable Temp"])**2
    for j in range(2*nSamples):
        denom=denom+(fC.loc[j]["Time to Stable Temp"]**2)
    num=num/nSamples
    numT=numT/(2*nSamples)
    denom=denom/(2*nSamples)
    denom=denom-(meanC)**2
    sT[i]=numT/denom
    s[i]=num/denom
    
print(s)
print(sT)