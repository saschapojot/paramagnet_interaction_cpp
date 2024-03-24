import numpy as np
import matplotlib.pyplot as plt
import re
import glob
from pathlib import Path
#this script plots chi and confidence interval

groupNum=0
rowNum=4

inDir="./group"+str(groupNum)+"data/row"+str(rowNum)+"/"
inChiFilesDir=inDir+"chiAll"

inChiFileNames=[]
TVals=[]
for file in glob.glob(inChiFilesDir+"/*.txt"):
    inChiFileNames.append(file)
    matchT=re.search(r"T([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)chi",file)
    TVals.append(float(matchT.group(1)))

T_inds=np.argsort(TVals)

TSortedVals=[TVals[ind] for ind in T_inds]

inChiFileNamesSorted=[inChiFileNames[ind] for ind in T_inds]

chiValsAll=[]
hfLengthAll=[]

for file in inChiFileNamesSorted:
    fptr=open(file,"r")
    contents=fptr.readlines()
    for line in contents:
        # print(line)
        matchChi=re.search(r"chi=(-?\d+(\.\d+)?)",line)
        if matchChi:
            chiTmp=float(matchChi.group(1))
            chiValsAll.append(chiTmp)
        matchhfLength=re.search(r"hfLength=(-?\d+(\.\d+)?([eE][-+]?\d+)?)",line)
        if matchhfLength:
            hfTmp=float(matchhfLength.group(1))
            hfLengthAll.append(hfTmp)




fullLengthAll=[2*val for val in hfLengthAll]
# for i in range(0,len(TSortedVals)):
#     print("T="+str(TSortedVals[i])+", chi="+str(chiValsAll[i]))
plt.figure()
# plt.scatter(TSortedVals,chiValsAll,color="black")
plt.errorbar(TSortedVals,chiValsAll,yerr=fullLengthAll,fmt="o",ecolor="red",color="black",markersize=2)
# print(fullLengthAll)
plt.xlabel("$T$")
plt.ylabel("$\chi$")
xTicks=np.linspace(np.min(TSortedVals),np.max(TSortedVals),10)
plt.xticks(xTicks)
outDir=inDir+"Observables"
Path(outDir).mkdir(exist_ok=True,parents=True)
fileName="chi.png"

plt.savefig(outDir+"/"+fileName)