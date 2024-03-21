import numpy as np
import matplotlib.pyplot as plt
import re
import glob
from pathlib import Path
import xml.etree.ElementTree as ET
from datetime import datetime
#this script plots unfolded bands and confidence interval


groupNum=0
rowNum=0

inPartDir="./group"+str(groupNum)+"data/row"+str(rowNum)+"/"
inTDirs=[]
TVals=[]

for TDir in glob.glob(inPartDir+"T*"):
    inTDirs.append(TDir)
    matchT=re.search(r"T(\d+(\.\d+)?)",TDir)
    TVals.append(float(matchT.group(1)))


def xmlFileNames(TDir):
    """

    :param TDir: directory corresponding to one T
    :return: path of xml files containing information of unfolded bands
    """
    # print(TDir)

    for fileName in glob.glob(TDir+"/*.xml"):
        # print(fileName)
        matchHfLength=re.search(r"EHfLength",fileName)
        if matchHfLength:
            hfLengthFile=fileName

        matchEMean=re.search(r"EMean",fileName)
        if matchEMean:
            EMeanFile=fileName

        matchMarker=re.search(r"marker",fileName)
        if matchMarker:
            markerFile=fileName


    return hfLengthFile,EMeanFile,markerFile


def xml2Vec(TDir):
    """

    :param TDir: directory corresponding to one T
    :return: flattened data for half length, mean of E, error
    """
    hfLengthFile,EMeanFile,markerFile=xmlFileNames(TDir)

    #parse hflength
    tree0=ET.parse(hfLengthFile)
    root0=tree0.getroot()
    flHfLength=root0.find("flattenedEHfLength")
    hf_items=flHfLength.findall("item")
    hfLengthFlattened=[float(item.text) for item in hf_items]

    #parse EMean
    tree1=ET.parse(EMeanFile)
    root1=tree1.getroot()
    flEMean=root1.find("flatttenedEMean")
    mean_items=flEMean.findall("item")
    EMeanFlattened=[float(item.text) for item in mean_items]

    #parse marker size
    tree2=ET.parse(markerFile)
    root2=tree2.getroot()
    flMarkerSz=root2.find("flattenedSize")
    sz_items=flMarkerSz.findall("item")
    mkSzFlattened=[float(item.text) for item in sz_items]

    return hfLengthFlattened,EMeanFlattened,mkSzFlattened


L=10
M=20
vecLength=4*L

def vec2Array(flattenedVec):
    """

    :param flattenedVec:
    :return: deserialized data
    """

    retMat=np.zeros((M,L,vecLength),dtype=float)
    lenForEach_m=L*vecLength
    for m in range(0,M):
        for a in range(0,L):
            startingPos=m*lenForEach_m+a*vecLength
            for j in range(0,vecLength):
                retMat[m,a,j]=flattenedVec[startingPos+j]

    return retMat
#0..6
# ind=3
# print(inTDirs)
for ind in range(0,len(inTDirs)):
    TD=inTDirs[ind]
    print(TD)
    hfLengthFlattened,EMeanFlattened,mkSzFlattened=xml2Vec(TD)

    matEMean=vec2Array(EMeanFlattened)
    matHfLength=vec2Array(hfLengthFlattened)
    matMkSize=vec2Array(mkSzFlattened)

    NPrim=M*L
    kPrimValsAll=[2*j/NPrim for j in range(0,NPrim)]

    matEMean2Plt=np.zeros((M*L,vecLength),dtype=float)
    for m in range(0,M):
        for a in range(0,L):
            rowPos=m+a*M
            matEMean2Plt[rowPos,:]=matEMean[m,a,:]

    matFullLength2Plt=np.zeros((M*L,vecLength),dtype=float)
    for m in range(0,M):
        for a in range(0,L):
            rowPos=m+a*M
            matFullLength2Plt[rowPos,:]=2*matHfLength[m,a,:]

    matMkSize2Plt=np.zeros((M*L,vecLength),dtype=float)
    for m in range(0,M):
        for a in range(0,L):
            rowPos=m+a*M
            matMkSize2Plt[rowPos,:]=matMkSize[m,a,:]

    plt.figure()

    rowNum,colNum=matEMean2Plt.shape
    print(rowNum)
    plt.figure()
    for c in range(0,colNum):
        EValsTmp=matEMean2Plt[:,c]
        mkSizeTmp=matMkSize2Plt[:,c]
        lengthTmp=matFullLength2Plt[:,c]
        lengthTmpNew=[]
        # plt.scatter(kPrimValsAll,EValsTmp,s=mkSizeTmp,color="red")
        # plt.errorbar(kPrimValsAll,EValsTmp,yerr=lengthTmp,ls="None",color="red")
        for i in range(0,len(lengthTmp)):
            val=0 if mkSizeTmp[i]<1e-3 else lengthTmp[i]
            lengthTmpNew.append(val)
        # print(lengthTmpNew)
        plt.errorbar(kPrimValsAll,EValsTmp,yerr=lengthTmpNew,ls="None",color="red",elinewidth=0.1)
        plt.scatter(kPrimValsAll,EValsTmp,s=mkSizeTmp,color="blue")


    plt.xlabel("$k/\pi$")
    plt.ylabel("$\epsilon$")
    plt.title("$T=$"+str(TVals[ind]))
    plt.savefig(TD+"/band.png")
    plt.close()

