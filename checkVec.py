import xml.etree.ElementTree as ET
import numpy as np
from datetime import datetime
import statsmodels.api as sm
from scipy import stats
import glob
import sys
import re
from copy import deepcopy
import warnings

#This script checks if a vector reaches equilibrium

sigWrongArg="1"
sigFerro="2"
sigContinue="3"
sigStop="0"

if (len(sys.argv)!=2):
    raise RuntimeError(sigWrongArg)

xmlFilesPath=str(sys.argv[1])

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

inXMLFileNames=[inXMLFileNames[ind] for ind in start_inds]

if len(inXMLFileNames)%3==0:
    xmlFileToBeParsed=deepcopy(inXMLFileNames)
elif len(inXMLFileNames)%3==1:
    xmlFileToBeParsed=deepcopy(inXMLFileNames[1:])
else:
    xmlFileToBeParsed=deepcopy(inXMLFileNames[2:])

xmlFileToBeParsed=xmlFileToBeParsed[int(len(xmlFileToBeParsed)/3):]


def parse1File(fileName):
    tree=ET.parse(fileName)
    root = tree.getroot()
    vec=root.find("vec")
    vec_items=vec.findall('item')
    vecValsAll=[float(item.text) for item in vec_items]
    # vecValsAll=np.array(vecValsAll)
    return vecValsAll


vecValsCombined=parse1File(xmlFileToBeParsed[0])

for file in xmlFileToBeParsed[1:]:
    vecValsCombined+=parse1File(file)


#ferromagnetic case: all values equal
vecValsCombined=np.array(vecValsCombined)


with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        vecAutc=sm.tsa.acf(vecValsCombined)
    except Warning as w:
        print(sigFerro)



halfLength=int(len(vecValsCombined)/2)

part0=vecValsCombined[:halfLength]
part1=vecValsCombined[halfLength:]

with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        autc0=sm.tsa.acf(part0)
    except Warning as w:
        ferro0=True

with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        autc1=sm.tsa.acf(part1)
    except Warning as w:
        ferro1=True

if ferro0 or ferro1:
    print(sigFerro)




