import xml.etree.ElementTree as ET
import numpy as np
import glob
import sys
import re
import statsmodels.api as sm
import matplotlib.pyplot as plt
# from copy import deepcopy

#This script loads EAll and sAll files to compute EAvg and sAvg, chi

L = 10# length of a supercell
M = 20#number of supercells

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()

pathPart=str(sys.argv[1])

inTFileNames=[]
TVals=[]
for file in glob.glob(pathPart+"/*"):
    inTFileNames.append(file)
    matchT=re.search(r"T(\d+(\.\d+)?)",file)
    TVals.append(float(matchT.group(1)))

T_inds=np.argsort(TVals)

TSortedVals=[TVals[ind] for ind in T_inds]
inTFileNamesSorted=[inTFileNames[ind] for ind in T_inds]


def parseSummaryFerro(summaryFile):
    """

    :param summaryFile:
    :return: whether the system is labeled as ferro in the summary file
    """
    fptr=open(summaryFile,"r")
    contents=fptr.readlines()

    for line in contents:
        match=re.search(r"ferro:\s*(\d+)",line)
        if match:
            ferro=int(match.group(1))
    return ferro


def EAndSFilesSelected(oneTFile):
    """

    :param oneTFile:
    :return: E files and s files to be parsed
    """
    smrFile=oneTFile+"/summary.txt"
    ferro=parseSummaryFerro(smrFile)
    fileNumSelected=0
    if ferro==1:
        fileNumSelected=1

    else:
        fileNumSelected=5
    EAllDir=oneTFile+"/EAll/*"
    sAllDir=oneTFile+"/sAll/*"

    inEAllFileNames=[]
    startEAllVals=[]

    for file in glob.glob(EAllDir):
        inEAllFileNames.append(file)
        matchEStart=re.search(r"loopStart(-?\d+(\.\d+)?)loopEnd",file)
        if matchEStart:
            startEAllVals.append(int(matchEStart.group(1)))

    start_E_inds=np.argsort(startEAllVals)
    sortedEAllFileNames=[inEAllFileNames[ind] for ind in start_E_inds]
    # print(sortedEAllFileNames)

    inSAllFileNames=[]
    startSAllVals=[]
    for file in glob.glob(sAllDir):
        inSAllFileNames.append(file)
        matchSAllStart=re.search(r"loopStart(-?\d+(\.\d+)?)loopEnd",file)
        if matchSAllStart:
            startSAllVals.append(int(matchSAllStart.group(1)))

    start_sAll_inds=np.argsort(startSAllVals)
    sortedSAllFilenames=[inSAllFileNames[ind] for ind in start_sAll_inds]

    # print(sortedSAllFilenames)
    retEAllFileNames=sortedEAllFileNames[-fileNumSelected:]
    retSAllFileNames=sortedSAllFilenames[-fileNumSelected:]
    return ferro, retEAllFileNames,retSAllFileNames






# EAllFilesSelected,sAllFilesSelected=EAndSFilesSelected(inTFileNamesSorted[20])

def parseEFile(EAllFileName):
    """

    :param EAllFileName: xml file containing E
    :return: values of E in this file
    """

    tree=ET.parse(EAllFileName)
    root = tree.getroot()
    vec=root.find("vec")
    vec_items=vec.findall('item')
    vecValsAll=[float(item.text) for item in vec_items]
    # vecValsAll=np.array(vecValsAll)
    return vecValsAll


def parseSFile(sAllFileName):
    """

    :param sAllFileName: xml file containing s vectors
    :return: all vectors in this xml file
    """
    tree=ET.parse(sAllFileName)
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


def combineValues(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: combined values of E and s from each file
    """
    ferro, ELastFileNames,sLastFileNames=EAndSFilesSelected(oneTFile)
    EVecValsCombined=parseEFile(ELastFileNames[0])
    for file in ELastFileNames[1:]:
        EVecValsCombined+=parseEFile(file)
    sMeanAbsVecCombined=sAbsAvg(parseSFile(sLastFileNames[0]))
    for file in  sLastFileNames[1:]:
        sMeanAbsNext=sAbsAvg(parseSFile(file))
        sMeanAbsVecCombined=np.r_[sMeanAbsVecCombined,sMeanAbsNext]

    return ferro,EVecValsCombined,sMeanAbsVecCombined



def lagVal(oneTFile):
    ferro,EVecValsCombined,sMeanAbsVecCombined=combineValues(oneTFile)
    print(oneTFile)
    if ferro==1:
        lag=-1
        return ferro,EVecValsCombined,sMeanAbsVecCombined,lag,0.05
    else:
        eps=0.05
        acfOfEVec=np.abs(sm.tsa.acf(EVecValsCombined))
        acfOfSVec=np.abs(sm.tsa.acf(sMeanAbsVecCombined))
        # print("ferro="+str(ferro))
        # print(np.min(acfOfEVec))
        # print(np.min(acfOfSVec))
        # print("================")
        minAcfEAll=np.min(acfOfEVec)
        minAcfSAll=np.min(acfOfSVec)
        # print(minAcfEAll)
        # print(minAcfSAll)
        if minAcfEAll>eps or minAcfSAll>eps:
            eps=np.max([minAcfEAll,minAcfSAll])
        # print(len(EVecValsCombined))
        lagEVec=np.where(acfOfEVec<=eps)[0][0]
        lagSVec=np.where(acfOfSVec<=eps)[0][0]
        lag=np.max([lagEVec,lagEVec])
        # return ferro,EVecValsCombined,sMeanAbsVecCombined,lag

        return ferro,EVecValsCombined,sMeanAbsVecCombined,lag,eps




def diagnosticsAndObservables(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: diagnostic plots and observable values
    """

    #diagnostics
    ferro,EVecValsCombined,sMeanAbsVecCombined,lag,eps=lagVal(oneTFile)
    TTmpMatch=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    if TTmpMatch:
        TTmp=float(TTmpMatch.group(1))
    if ferro==1:
        plt.figure()
        EPerSupercell=np.array(EVecValsCombined)/M
        plt.hist(EPerSupercell)
        sdE=np.sqrt(np.var(EPerSupercell)/len(EPerSupercell))
        meanE=np.mean(EPerSupercell)
        plt.title("ferromagnetic, T="+str(TTmp)+", mean="+str(meanE)+", sd="+str(sdE),fontsize=10)
        plt.xlabel("$E$")
        plt.savefig(oneTFile+"/EHist.png")
        plt.close()

        plt.figure()
        plt.hist(sMeanAbsVecCombined)
        meanS=np.mean(sMeanAbsVecCombined)
        sdS=np.sqrt(np.var(sMeanAbsVecCombined)/len(sMeanAbsVecCombined))
        plt.title("ferromagnetic, T="+str(TTmp)+", mean="+str(meanS)+", sd="+str(sdS),fontsize=10)
        plt.xlabel("$|s|$")
        plt.savefig(oneTFile+"/sHist.png")
        plt.close()
    else:
        halfLength=int(len(EVecValsCombined)/2)
        EVecPart0=EVecValsCombined[:halfLength]
        EvecPart1=EVecValsCombined[halfLength:]

        sVecPart0=sMeanAbsVecCombined[:halfLength]
        sVecPart1=sMeanAbsVecCombined[halfLength:]




diagnosticsAndObservables(inTFileNamesSorted[10])