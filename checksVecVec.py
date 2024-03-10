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
# import matplotlib.pyplot as plt


#This script checks if a vector<vector> reaches equilibrium
#applicable to s values

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

xmlFileToBeParsed=xmlFileToBeParsed[int(len(xmlFileToBeParsed)/3*2):]


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


def sAbsAvg(sVecVec):
    """

    :param sVecVec: an array of s values from mc steps
    :return: abs of avg of each row
    """
    sMeanTmp=np.mean(sVecVec,axis=1)
    sAbs=np.abs(sMeanTmp)
    return sAbs

vecValsCombined=sAbsAvg(parse1File(xmlFileToBeParsed[0]))

for file in xmlFileToBeParsed[1:]:
    sAbsNext=sAbsAvg(parse1File(file))
    vecValsCombined=np.r_[vecValsCombined,sAbsNext]
# print(len(vecValsCombined))
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
def Jackknife(vec):
    """

    :param vec:
    :return: the mean and half length  of 0.95 confidence interval computed by Jackkknife
    """

    psMean=np.mean(vec)

    psVar=np.var(vec,ddof=1)

    n=len(vec)

    hfLen=1.96*np.sqrt(psVar/n)
    return psMean,hfLen


acfOfVec=sm.tsa.acf(vecValsCombined)
# plt.figure()
# plt.plot(acfOfVec,color="black")
# plt.savefig("sAutc.png")
# plt.close()
eps=(1e-2)*5
pThreshHold=0.05
lagVal=0
if np.min(np.abs(acfOfVec))>eps:
    # print("high correlation")
    print(sigContinue)
    exit()
else:
    lagVal=np.where(np.abs(acfOfVec)<=eps)[0][0]
    selectedFromPart0=part0[::lagVal]
    selectedFromPart1=part1[::lagVal]

    # plt.subplot(1,2,1)
    # plt.hist(selectedFromPart0,bins=100)
    # plt.subplot(1,2,2)
    # plt.hist(selectedFromPart1,bins=100)
    # print(np.mean(selectedFromPart0), np.sqrt(np.var(selectedFromPart0))/np.sqrt(len(selectedFromPart0/60)))
    # print(np.mean(selectedFromPart1), np.sqrt(np.var(selectedFromPart1))/np.sqrt(len(selectedFromPart1/60)))
    # plt.savefig("histS.png")

    # D,p=stats.ks_2samp(part0,part1)
    mean0,hf0=Jackknife(part0)
    mean1,hf1=Jackknife(part1)
    if np.abs(mean0-mean1)<=hf0 or np.abs(mean0-mean1)<=hf1:
        print(sigEq+" "+str(lagVal))
        exit()
    else:
        print(sigContinue)
        exit()



