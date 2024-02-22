import xml.etree.ElementTree as ET
import numpy as np
# from datetime import datetime
import statsmodels.api as sm
from scipy import stats
import glob
import sys
import re
from copy import deepcopy
import warnings

#This script checks if a vector reaches equilibrium

sigWrongArg="wrong number of arguments"
sigEq="equilibrium"
sigContinue="continue"
sigStop="stop"

if (len(sys.argv)!=2):
    print(sigWrongArg)
    exit()

xmlFilesPath=str(sys.argv[1])

#fetch files in the directory
inXMLFileNames=[]
startVals=[]
for file in glob.glob(xmlFilesPath+"/*"):
    inXMLFileNames.append(file)
    matchStart=re.search(r"loopStart(-?\d+(\.\d+)?)loopEnd",file)
    if matchStart:
        startVals.append(matchStart.group(1))


def str2int(valList):
    ret = [int(strTmp) for strTmp in valList]
    return ret


startVals = str2int(startVals)

start_inds = np.argsort(startVals)

#sort files by starting value
inXMLFileNames=[inXMLFileNames[ind] for ind in start_inds]


#ensure the file number is a multiple of 3
if len(inXMLFileNames)%3==0:
    xmlFileToBeParsed=deepcopy(inXMLFileNames)
elif len(inXMLFileNames)%3==1:
    xmlFileToBeParsed=deepcopy(inXMLFileNames[1:])
else:
    xmlFileToBeParsed=deepcopy(inXMLFileNames[2:])

xmlFileToBeParsed=xmlFileToBeParsed[int(len(xmlFileToBeParsed)/3):]

# print("xml file number: "+str(len(xmlFileToBeParsed)))
def parse1File(fileName):
    """

    :param fileName: xml file name
    :return: the values in the vector
    """

    tree=ET.parse(fileName)
    root = tree.getroot()
    vec=root.find("vec")
    vec_items=vec.findall('item')
    vecValsAll=[float(item.text) for item in vec_items]
    # vecValsAll=np.array(vecValsAll)
    return vecValsAll

#combine all vectors
vecValsCombined=parse1File(xmlFileToBeParsed[0])

for file in xmlFileToBeParsed[1:]:
    vecValsCombined+=parse1File(file)



vecValsCombined=np.array(vecValsCombined)
# print(len(vecValsCombined))
#ferromagnetic case: all values equal
#check if the whole vector has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        vecAutc=sm.tsa.acf(vecValsCombined)
    except Warning as w:
        print(sigStop+" ferro")
        exit()



halfLength=int(len(vecValsCombined)/2)

part0=vecValsCombined[:halfLength]
part1=vecValsCombined[halfLength:]

ferro0=False
ferro1=False
#check if the part0 has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        autc0=sm.tsa.acf(part0)
    except Warning as w:
        ferro0=True
#check if the part1 has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        autc1=sm.tsa.acf(part1)
    except Warning as w:
        ferro1=True

if ferro0 and ferro1:
    print(sigStop+" ferro")
    exit()
elif ferro0==True and ferro1==False:
    # print("f0 True f1 False")
    print(sigContinue)
    exit()
elif ferro0==False and ferro1==True:
    # print("f0 False f1 True")
    print(sigContinue)
    exit()



#computation of auto-correlation



acfOfVec=sm.tsa.acf(vecValsCombined)

# print("total elem number = "+str(len(vecValsCombined)))
eps=(1e-3)*5
pThreshHold=0.05
lagVal=0
if np.min(acfOfVec)>eps:
    # print("high correlation")
    print(sigContinue)
    exit()
else:
    lagVal=np.where(acfOfVec<=eps)[0][0]
    selectedFromPart0=part0[::lagVal]
    selectedFromPart1=part1[::lagVal]
    D,p=stats.ks_2samp(part0,part1)
    if p<pThreshHold:
        print(sigContinue)
        # print(p)
        exit()
    else:
        print(sigEq+" "+str(lagVal))
        exit()



