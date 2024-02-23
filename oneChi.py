import numpy as np
from datetime import datetime
import glob
import re
import xml.etree.ElementTree as ET
from copy import deepcopy
import sys
#This script computes chi for one temperature

if len(sys.argv)!=2:
    raise RuntimeError("Wrong number of arguments")


#parse summary.txt

summaryPath=sys.argv[1]+"/summary.txt"
summaryFileFptr=open(summaryPath,"r")
contents=summaryFileFptr.readlines()


lag=0
ferro=0
reachEq=0

for line in contents:
    matchEq=re.search(r"equilibrium reached:\s*(\d)",line)
    if matchEq:
        reachEq=int(matchEq.group(1))
    matchFerro=re.search(r"ferro:\s*(\d)",line)
    if matchFerro:
        ferro=int(matchFerro.group(1))
    matchLag=re.search(r"lag\s*=\s*(-?\d+)",line)
    if matchLag:
        lag=int(matchLag.group(1))


if reachEq!=1:
    raise RuntimeError("No equilibrium for file "+sys.argv[1])




xmlFilesPath=str(sys.argv[1])+"/sAll/"



inXMLFileNamesAll=[]
inXMLFileNamesBeforeEq=[]
inXMLFileNamesAfterEq=[]

for file in glob.glob(xmlFilesPath+"/*"):
    inXMLFileNamesAll.append(file)

for name in inXMLFileNamesAll:
    matchEq=re.search(r"AfterEq",name)
    if matchEq:
        inXMLFileNamesAfterEq.append(name)
    else:
        inXMLFileNamesBeforeEq.append(name)



startValsBeforeEq=[]
startValsAfterEq=[]

for fileName in inXMLFileNamesBeforeEq:
    matchStart=re.search(r"loopStart(-?\d+(\.\d+)?)loopEnd",fileName)
    if matchStart:
        startValsBeforeEq.append(matchStart.group(1))



for fileName in inXMLFileNamesAfterEq:
    matchStart=re.search(r"loopStart(-?\d+(\.\d+)?)loopEnd",fileName)
    if matchStart:
        startValsAfterEq.append(matchStart.group(1))

def str2int(valList):
    ret = [int(strTmp) for strTmp in valList]
    return ret
startValsBeforeEq=str2int(startValsBeforeEq)

startValsAfterEq=str2int(startValsAfterEq)

start_BeforeEq_inds=np.argsort(startValsBeforeEq)

start_AfterEq_inds=np.argsort(startValsAfterEq)
#sort files by starting value

inXMLFileNamesBeforeEq=[inXMLFileNamesBeforeEq[ind] for ind in start_BeforeEq_inds]

inXMLFileNamesAfterEq=[inXMLFileNamesAfterEq[ind] for ind in start_AfterEq_inds]

#ensure the file number before eq is a multiple of 3

if len(inXMLFileNamesBeforeEq)%3==0:
    xmlFileToBeParsed=deepcopy(inXMLFileNamesBeforeEq)
elif len(inXMLFileNamesBeforeEq)%3==1:
    xmlFileToBeParsed=deepcopy(inXMLFileNamesBeforeEq[1:])
else:
    xmlFileToBeParsed=deepcopy(inXMLFileNamesBeforeEq[2:])

xmlFileToBeParsed=xmlFileToBeParsed[int(len(xmlFileToBeParsed)/3):]

xmlFileToBeParsed+=inXMLFileNamesAfterEq

def parse1File(fileName):
    """

    :param fileName: xml file name
    :return: all vectors in the xml file
    """
    tree=ET.parse(fileName)
    root = tree.getroot()
    first_level_items = root.find('vecvec').findall('item')
    vectors=[]
    for item in first_level_items:
        oneVec=[float(value.text) for value in item.findall('item')]
        vectors.append(oneVec)
    return np.array(vectors)

sVecValsCombined=parse1File(xmlFileToBeParsed[0])
# print(sVecValsCombined)
# print(sVecValsCombined.shape)
for file in xmlFileToBeParsed[1:]:
    sVecsNext=parse1File(file)
    sVecValsCombined=np.r_[sVecValsCombined,sVecsNext]





if ferro==1 and lag<=0:
    sVecValsSelected=deepcopy(sVecValsCombined)
if ferro==0 and lag>0:
    sVecValsSelected=sVecValsCombined[0::lag]

matchT=re.search(r"T(\d+(\.\d+)?)",sys.argv[1])
if matchT:
    T=float(matchT.group(1))
# print("T="+str(T))
sMean=np.mean(sVecValsSelected,axis=1)#mean for each mc configuration
sAbsMean=np.abs(sMean)#absolute value of mean for each mc configuration
sAvgAbsMean=np.mean(sAbsMean)#mean of absolute value of mean for each mc configuration
sSquared=sAbsMean**2
meanS2=np.mean(sSquared)

chiVal=(meanS2-sAvgAbsMean**2)/T

print("chi="+str(chiVal))