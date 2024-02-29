import xml.etree.ElementTree as ET
import numpy as np
import glob
import sys
import re
import statsmodels.api as sm
import matplotlib.pyplot as plt
# from copy import deepcopy
from pathlib import Path
from multiprocessing import Pool
from datetime import datetime
#This script loads EAll and sAll files under one temperature to compute EAvg and sAvg, chi

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
    fileNumSelected=0#files' numbers to be parsed
    if ferro==1:
        fileNumSelected=1

    else:
        fileNumSelected=15
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
        # print(eps)
        lagEVec=np.where(acfOfEVec<=eps)[0][0]
        lagSVec=np.where(acfOfSVec<=eps)[0][0]
        lag=np.max([lagEVec,lagEVec])
        # return ferro,EVecValsCombined,sMeanAbsVecCombined,lag

        return ferro,EVecValsCombined,sMeanAbsVecCombined,lag,eps

def pseudoValueForChi(SVec,i,T):
    """

    :param SVec: a vector containing S values
    :param i: the index of the element to be deleted
    :param T: temperature
    :return: pseudovalue for chi with ith element deleted
    """
    SVec=np.array(SVec)
    SVecDeleted=[SVec[j] for j in range(0,len(SVec)) if j!=i]
    SVecDeleted=np.array(SVecDeleted)
    chi_i=np.mean(SVecDeleted**2)-(np.mean(SVecDeleted))**2
    chi_i*=1/T
    return chi_i

def JackknifeForChi(SVec,T):
    """

    :param SVec: a vector containing S values
    :param T: temperature
    :return: mean and half confidence interval for chi
    """
    chi_deletedAll=np.array([pseudoValueForChi(SVec,i,T) for i in range(0,len(SVec))])
    chi_ps_mean=np.mean(chi_deletedAll)
    chi_ps_var=np.sum((chi_deletedAll-chi_ps_mean)**2)/(len(SVec)-1)

    halfInterval=1.960*np.sqrt(chi_ps_var/len(SVec))

    return chi_ps_mean, halfInterval




def diagnosticsAndObservables(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: diagnostic plots and observable values
    """
    tOneFileStart=datetime.now()
    #diagnostics
    ferro,EVecValsCombined,sMeanAbsVecCombined,lag,eps=lagVal(oneTFile)

    TTmpMatch=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    matchPartNum=re.search(r"part(\d+)",oneTFile)
    if matchPartNum:
        partNum=int(matchPartNum.group(1))
    EHistAllDir="./part"+str(partNum)+"EHistAll"
    sHistAllDir="./part"+str(partNum)+"sHistAll"
    EBlkMeanDir="./part"+str(partNum)+"EBlkMean"
    Path(EHistAllDir).mkdir(parents=True, exist_ok=True)
    Path(sHistAllDir).mkdir(parents=True, exist_ok=True)
    Path(EBlkMeanDir).mkdir(parents=True,exist_ok=True)
    if TTmpMatch:
        TTmp=float(TTmpMatch.group(1))
    if ferro==1:
        nbins=100
        numElem=1000
        EVecValsCombined=EVecValsCombined[-numElem:]
        sMeanAbsVecCombined=sMeanAbsVecCombined[-numElem:]
        plt.figure()
        EPerSupercell=np.array(EVecValsCombined)/M
        plt.hist(EPerSupercell,bins=nbins)
        sdE=np.sqrt(np.var(EPerSupercell)/len(EPerSupercell))
        meanE=np.mean(EPerSupercell)
        plt.title("ferromagnetic, T="+str(TTmp)+", mean="+str(meanE)+", sd="+str(sdE),fontsize=10)
        plt.xlabel("$E$")
        EHistOut="T"+str(TTmp)+"EHist.png"
        plt.savefig(oneTFile+"/"+EHistOut)
        plt.savefig(EHistAllDir+"/"+EHistOut)
        plt.close()

        plt.figure()
        plt.hist(sMeanAbsVecCombined,bins=nbins)
        meanS=np.mean(sMeanAbsVecCombined)
        sdS=np.sqrt(np.var(sMeanAbsVecCombined)/len(sMeanAbsVecCombined))
        plt.title("ferromagnetic, T="+str(TTmp)+", mean="+str(meanS)+", sd="+str(sdS),fontsize=10)
        plt.xlabel("$|s|$")
        sHistOut="T"+str(TTmp)+"sHist.png"
        plt.savefig(oneTFile+"/"+sHistOut)
        plt.savefig(sHistAllDir+"/"+sHistOut)
        plt.close()

        #block mean
        def meanPerBlock(length):
            blockNum=int(len(EPerSupercell)/length)
            EMeanBlock=[]
            for blkNum in range(0,blockNum):
                blkE=EPerSupercell[blkNum*length:(blkNum+1)*length]
                EMeanBlock.append(np.mean(blkE))
            return EMeanBlock
        fig=plt.figure(figsize=(20,20))
        fig.tight_layout(pad=5.0)
        lengthVals=[20,50,100,300]
        for i in range(0,len(lengthVals)):
            # print("entering loop "+str(i))
            l=lengthVals[i]
            EMeanBlk=meanPerBlock(l)
            ax=fig.add_subplot(2,2,i+1)
            (n,_,_)=ax.hist(EMeanBlk,bins=100,color="aqua")
            xPosTextBlk=(np.max(EMeanBlk)-np.min(EMeanBlk))*1/7+np.min(EMeanBlk)
            yPosTextBlk=np.max(n)*3/4

            meanTmp=np.mean(EMeanBlk)
            meanTmp=np.round(meanTmp,3)
            sdTmp=np.sqrt(np.var(EMeanBlk))
            sdTmp=np.round(sdTmp,3)
            ax.set_title("L="+str(l))
            ax.text(xPosTextBlk,yPosTextBlk,"mean="+str(meanTmp)+", sd="+str(sdTmp))
        fig.suptitle("T="+str(TTmp)+", ferromagnetic")
        plt.savefig(oneTFile+"/T"+str(TTmp)+"EBlk.png")
        plt.savefig(EBlkMeanDir+"/T"+str(TTmp)+"EBlk.png")
        plt.close()
        #observables
        chi_ps,hfInterval=JackknifeForChi(sMeanAbsVecCombined,TTmp)
        chiOutFileName="T"+str(TTmp)+"chi.txt"
        contents=["chi="+str(chi_ps)+"\n","hfLength="+str(hfInterval)+"\n"+"lag="+str(lag)]
        fptr1=open(oneTFile+"/"+chiOutFileName,"w+")
        fptr1.writelines(contents)
        fptr1.close()

        chiAllDir="./part"+str(partNum)+"chiAll"
        Path(chiAllDir).mkdir(exist_ok=True,parents=True)
        fptr2=open(chiAllDir+"/"+chiOutFileName,"w+")
        fptr2.writelines(contents)
        fptr2.close()
    # else:
    #
    #     EPerSupercell=np.array(EVecValsCombined)/M
    #
    #
    #
    #     #block mean
    #     def meanPerBlock(length):
    #         blockNum=int(len(EPerSupercell)/length)
    #         EMeanBlock=[]
    #         for blkNum in range(0,blockNum):
    #             blkE=EPerSupercell[blkNum*length:(blkNum+1)*length]
    #             EMeanBlock.append(np.mean(blkE))
    #         return EMeanBlock
    #     fig=plt.figure(figsize=(20,20))
    #     fig.tight_layout(pad=5.0)
    #     lengthVals=[20,50,100,300]
    #     for i in range(0,len(lengthVals)):
    #         # print("entering loop "+str(i))
    #         l=lengthVals[i]
    #         EMeanBlk=meanPerBlock(l)
    #         ax=fig.add_subplot(2,2,i+1)
    #         (n,_,_)=ax.hist(EMeanBlk,bins=100,color="aqua")
    #         xPosTextBlk=(np.max(EMeanBlk)-np.min(EMeanBlk))*1/7+np.min(EMeanBlk)
    #         yPosTextBlk=np.max(n)*3/4
    #
    #         meanTmp=np.mean(EMeanBlk)
    #         meanTmp=np.round(meanTmp,3)
    #         sdTmp=np.sqrt(np.var(EMeanBlk))
    #         sdTmp=np.round(sdTmp,3)
    #         ax.set_title("L="+str(l))
    #         ax.text(xPosTextBlk,yPosTextBlk,"mean="+str(meanTmp)+", sd="+str(sdTmp))
    #     fig.suptitle("T="+str(TTmp)+", corr="+str(np.round(eps,3)))
    #     plt.savefig(oneTFile+"/T"+str(TTmp)+"EBlk.png")
    #     plt.savefig(EBlkMeanDir+"/T"+str(TTmp)+"EBlk.png")
    #     plt.close()








        halfLength=int(len(EPerSupercell)/2)
        #
        # #histogram of distribution of epsilon
        EVecPart0=EPerSupercell[:halfLength]
        EvecPart1=EPerSupercell[halfLength:]

        sVecPart0=sMeanAbsVecCombined[:halfLength]
        sVecPart1=sMeanAbsVecCombined[halfLength:]
        #
        # ESelectedFromPart0=EVecPart0[::lag]
        # ESelectedFromPart1=EvecPart1[::lag]
        #
        #
        # #diagnostics of E
        # meanEPart0=np.mean(ESelectedFromPart0)
        # varEPart0=np.var(ESelectedFromPart0)
        # sdEPart0=np.sqrt(varEPart0/len(ESelectedFromPart0))
        #
        # meanEPart1=np.mean(ESelectedFromPart1)
        # varEPart1=np.var(ESelectedFromPart1)
        # sdEPart1=np.sqrt(varEPart1/len(ESelectedFromPart1))
        #
        # nbins=100
        #
        # #histogram of E's part0 and E's part1
        # fig=plt.figure()
        # axE0=fig.add_subplot(1,2,1)
        # (n0,_,_)= axE0.hist(ESelectedFromPart0,bins=nbins)
        # meanEPart0=np.round(meanEPart0,4)
        # sdEPart0=np.round(sdEPart0,4)
        # axE0.set_title("part0, T="+str(np.round(TTmp,3)))
        # axE0.set_xlabel("$\epsilon$")
        # axE0.set_ylabel("#")
        # xPosE0Text=(np.max(ESelectedFromPart0)-np.min(ESelectedFromPart1))*1/2+np.min(ESelectedFromPart0)
        # yPosE0Text=np.max(n0)*2/3
        # axE0.text(xPosE0Text,yPosE0Text,"mean="+str(meanEPart0)+"\nsd="+str(sdEPart0)+"\nlag="+str(lag)+"\ncorr="+str(np.round(eps,4)))
        #
        # axE1=fig.add_subplot(1,2,2)
        # (n1,_,_)=axE1.hist(ESelectedFromPart1,bins=nbins)
        # meanEPart1=np.round(meanEPart1,4)
        # sdEPart1=np.round(sdEPart1,4)
        # axE1.set_title("part1, T="+str(np.round(TTmp,3)))
        # axE1.set_xlabel("$\epsilon$")
        # axE1.set_ylabel("#")
        # xPosE1Text=(np.max(ESelectedFromPart1)-np.min(ESelectedFromPart1))*1/2+np.min(ESelectedFromPart1)
        # yPosE1Text=np.max(n1)*2/3
        # axE1.text(xPosE1Text,yPosE1Text,"mean="+str(meanEPart1)+"\nsd="+str(sdEPart1)+"\nlag="+str(lag)+"\ncorr="+str(np.round(eps,4)))
        # EHistOut="T"+str(TTmp)+"EHist.png"
        #
        # plt.savefig(oneTFile+"/"+EHistOut)
        # plt.savefig(EHistAllDir+"/"+EHistOut)
        # plt.close()
        #
        # #diagnostics of s
        sSelectedFromPart0=sVecPart0[::lag]
        sSelectedFromPart1=sVecPart1[::lag]

        meanSPart0=np.mean(sSelectedFromPart0)
        varSPart0=np.var(sSelectedFromPart0)
        sdSPart0=np.sqrt(varSPart0/len(sSelectedFromPart0))

        meanSPart1=np.mean(sSelectedFromPart1)
        varSPart1=np.var(sSelectedFromPart1)
        sdSPart1=np.sqrt(varSPart1/len(sSelectedFromPart1))

        #histogram of s's part0 and s's part1
        # fig=plt.figure()
        # axS0=fig.add_subplot(1,2,1)
        # (n0,_,_)=axS0.hist(sSelectedFromPart0,bins=nbins)
        # meanSPart0=np.round(meanSPart0,4)
        # sdSPart0=np.round(sdSPart0,4)
        # axS0.set_title("part0, T="+str(np.round(TTmp,3)))
        # axS0.set_xlabel("$|s|$")
        # axS0.set_ylabel("#")
        # xPosS0Text=(np.max(sSelectedFromPart0)-np.min(sSelectedFromPart0))*1/2+np.min(sSelectedFromPart0)
        # yPosS0Text=np.max(n0)*2/3
        # axS0.text(xPosS0Text,yPosS0Text,"mean="+str(meanSPart0)+"\nsd="+str(sdSPart0)+"\nlag="+str(lag)+"\ncorr="+str(np.round(eps,4)))
        #
        # axS1=fig.add_subplot(1,2,2)
        # (n1,_,_)=axS1.hist(sSelectedFromPart1,bins=nbins)
        # meanSPart1=np.round(meanSPart1,4)
        # sdSPart1=np.round(sdSPart1,4)
        # axS1.set_title("part1, T="+str(np.round(TTmp,3)))
        # axS1.set_xlabel("$|s|$")
        # axS1.set_ylabel("#")
        # xPosS1Text=(np.max(sSelectedFromPart1)-np.min(sSelectedFromPart1))*1/2+np.min(sSelectedFromPart1)
        # yPosS1Text=np.max(n1)*2/3
        # axS1.text(xPosS1Text,yPosS1Text,"mean="+str(meanSPart1)+"\nsd="+str(sdSPart1)+"\nlag="+str(lag)+"\ncorr="+str(np.round(eps,4)))
        # sHistOut="T"+str(TTmp)+"sHist.png"
        # plt.savefig(oneTFile+"/"+sHistOut)
        # plt.savefig(sHistAllDir+"/"+sHistOut)
        # plt.close()

        #observavles

        #chi
        SVec=np.r_[sSelectedFromPart0,sSelectedFromPart1]
    #     # print(SVec)
    #     # print(len(SVec))
    #
    #
        chi_ps,hfInterval=JackknifeForChi(SVec,TTmp)
    #     # print(chi_ps)
    #     # print(hfInterval)
    #
        chiOutFileName="T"+str(TTmp)+"chi.txt"
        contents=["chi="+str(chi_ps)+"\n","hfLength="+str(hfInterval)+"\n"+"lag="+str(lag)]
        fptr1=open(oneTFile+"/"+chiOutFileName,"w+")
        fptr1.writelines(contents)
        fptr1.close()

        chiAllDir="./part"+str(partNum)+"chiAll"
        Path(chiAllDir).mkdir(exist_ok=True,parents=True)
        fptr2=open(chiAllDir+"/"+chiOutFileName,"w+")
        fptr2.writelines(contents)
        fptr2.close()
    tOneFileEnd=datetime.now()
    print("one file time: ",tOneFileEnd-tOneFileStart)







# diagnosticsAndObservables(inTFileNamesSorted[43])
tStart=datetime.now()
# procNum=48
# #parallel
# # pool0=Pool(procNum)
# #
# # ret=pool0.map(diagnosticsAndObservables,inTFileNamesSorted)
# #serial
for file in inTFileNamesSorted:
    diagnosticsAndObservables(file)
tEnd=datetime.now()
print("total time: ",tEnd-tStart)

