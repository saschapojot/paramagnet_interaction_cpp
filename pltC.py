import numpy as np
import matplotlib.pyplot as plt
import re
import glob
from pathlib import Path
#this script plots C and confidence interval


part=1


inCFilesDir="./part"+str(part)+"CAll"
inCFileNames=[]
TVals=[]

for file in glob.glob(inCFilesDir+"/*.txt"):
    inCFileNames.append(file)
    matchT=re.search(r"T(\d+(\.\d+)?)C",file)
    TVals.append(float(matchT.group(1)))



T_inds=np.argsort(TVals)

TSortedVals=[TVals[ind] for ind in T_inds]

inCFileNamesSorted=[inCFileNames[ind] for ind in T_inds]

CValsAll=[]
hfLengthAll=[]

for file in inCFileNamesSorted:
    fptr=open(file,"r")
    contents=fptr.readlines()
    for line in contents:
        matchC=re.search(r"C=(-?\d+(\.\d+)?)",line)
        if matchC:
            CTmp=float(matchC.group(1))
            CValsAll.append(CTmp)
        matchhfLength=re.search(r"halfLength=(-?\d+(\.\d+)?([eE][-+]?\d+)?)",line)
        if matchhfLength:
            hfTmp=float(matchhfLength.group(1))
            hfLengthAll.append(hfTmp)



fullLengthAll=[2*val for val in hfLengthAll]
print(len(TSortedVals))
print(len(CValsAll))
print(len(fullLengthAll))
plt.figure()
plt.errorbar(TSortedVals,CValsAll,yerr=fullLengthAll,fmt="o",ecolor="red",color="black",markersize=2)
plt.xlabel("$T$")
plt.ylabel("$C$")
xTicks=np.linspace(np.min(TSortedVals),np.max(TSortedVals),10)
plt.xticks(xTicks)
outDir="./part"+str(part)+"Observables"
Path(outDir).mkdir(exist_ok=True,parents=True)
fileName="C.png"

plt.savefig(outDir+"/"+fileName)