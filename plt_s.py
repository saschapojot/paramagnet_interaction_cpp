import numpy as np
import matplotlib.pyplot as plt
import re
import glob
from pathlib import Path
import sys
#this script plots s and confidence interval

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()
groupNum=int(sys.argv[1])
rowNum=int(sys.argv[2])
inDir="./group"+str(groupNum)+"data/row"+str(rowNum)+"/"
in_sFileDir=inDir+"/sAll"+"/"
in_sFileNames=[]

TVals=[]

for file in glob.glob(in_sFileDir+"/*.txt"):
    in_sFileNames.append(file)
    matchT=re.search(r"T([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)s",file)
    TVals.append(float(matchT.group(1)))


T_inds=np.argsort(TVals)

TSortedVals=[TVals[ind] for ind in T_inds]


in_sFileNamesSorted=[in_sFileNames[ind] for ind in T_inds]

sValsAll=[]
hfLengthAll=[]

for file in in_sFileNamesSorted:
    fptr=open(file,"r")
    contents=fptr.readlines()
    for line in contents:
        matchs=re.search(r"s=(-?\d+(\.\d+)?)",line)
        if matchs:
            sTmp=float(matchs.group(1))
            sValsAll.append(sTmp)
        matchhfLength=re.search(r"hfLength=(-?\d+(\.\d+)?([eE][-+]?\d+)?)",line)
        if matchhfLength:
            hfTmp=float(matchhfLength.group(1))
            hfLengthAll.append(hfTmp)


fullLengthAll=[2*val for val in hfLengthAll]
plt.figure()
plt.errorbar(TSortedVals,sValsAll,yerr=fullLengthAll,fmt="o",ecolor="red",color="black",markersize=2)
plt.xlabel("$T$")
plt.ylabel("$s$")

xTicks=np.linspace(np.min(TSortedVals),np.max(TSortedVals),10)
plt.xticks(xTicks)
outDir=inDir+"Observables"
Path(outDir).mkdir(exist_ok=True,parents=True)
fileName="s.png"

plt.savefig(outDir+"/"+fileName)