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

# if (len(sys.argv)!=2):
#     print("wrong number of arguments")
#     exit()
sweepNumInOneFlush=3000
groupNum=0
rowNum=0
pathGroupRow="./group"+str(groupNum)+"data/row"+str(rowNum)+"/"
inTFileNames=[]
TVals=[]
for file in glob.glob(pathGroupRow+"/T*"):
    # print(file)
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
    ferro=False
    sweepNumBeforeEq=-1
    lastFileNum=-1
    lag=-1
    for line in contents:
        match=re.search(r"ferro:\s*(\d+)",line)
        if match:
            ferro=int(match.group(1))
        matchSweepNum=re.search(r"total sweep number\s*:\s*(\d+)",line)
        if matchSweepNum:
            sweepNumBeforeEq=int(matchSweepNum.group(1))
        matchLag=re.search(r"lag=\s*(\d+)",line)
        if matchLag:
            lag=int(matchLag.group(1))
        matchLastFileNum=re.search(r"lastFileNum=\s*(\d+)",line)
        if matchLastFileNum:
            lastFileNum=int(matchLastFileNum.group(1))

    return ferro, sweepNumBeforeEq,lag,lastFileNum

# paraFileSelectedNum=15
def searchAfterEqFile(oneTFile):
    """

    :param oneTFile: one T directory
    :return: a list containing the summaryAfterEq.txt file, the list will be length 0 if the
    file does not exist
    """
    file=glob.glob(oneTFile+"/summaryAfterEq.txt")

    return file

def parseAfterEqFile(oneTFile):
    """

    :param oneTFile: one T directory
    :return: a list containing the summaryAfterEq.txt file, and the total sweep number after eq
    """
    fileList=searchAfterEqFile(oneTFile)
    sweepNumAfterEq=-1
    if len(fileList)!=0:
        fileName=fileList[0]
        fptr=open(fileName,"r")
        contents=fptr.readlines()
        for line in contents:
            match=re.search(r"total sweep number\s*:\s*(\d+)",line)
            if match:
                sweepNumAfterEq=int(match.group(1))
    return fileList,sweepNumAfterEq


def EAndSFilesSelected(oneTFile):
    """

    :param oneTFile: one T directory
    :return: E files and s files to be parsed
    """
    smrFile=oneTFile+"/summary.txt"
    ferro,sweepNumBeforeEq,lag,lastFileNum=parseSummaryFerro(smrFile)
    fileAfterEqList,sweepNumAfterEq=parseAfterEqFile(oneTFile)
    fileNumSelected=0#files' numbers to be parsed
    if ferro==1:
        fileNumSelected=1

    else:
        if len(fileAfterEqList)!=0:
            sweepNumToInclude=sweepNumInOneFlush*lastFileNum+sweepNumAfterEq
            fileNumSelected=int(np.ceil(sweepNumToInclude/sweepNumInOneFlush))
        else:
            fileNumSelected=lastFileNum
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
    return ferro, retEAllFileNames,retSAllFileNames,lag,fileNumSelected






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
    :return: combined values of E and s from each file, names of the parsed files
    """
    ferro, ELastFileNames,sLastFileNames,lag,fileNumSelected=EAndSFilesSelected(oneTFile)
    EVecValsCombined=parseEFile(ELastFileNames[0])
    for file in ELastFileNames[1:]:
        EVecValsCombined+=parseEFile(file)
    sMeanAbsVecCombined=sAbsAvg(parseSFile(sLastFileNames[0]))
    for file in  sLastFileNames[1:]:
        sMeanAbsNext=sAbsAvg(parseSFile(file))
        sMeanAbsVecCombined=np.r_[sMeanAbsVecCombined,sMeanAbsNext]

    return ferro,EVecValsCombined,sMeanAbsVecCombined,lag,fileNumSelected



# def lagVal(oneTFile):
#     ferro,EVecValsCombined,sMeanAbsVecCombined,lag=combineValues(oneTFile)
#     print(oneTFile)
#     if ferro==1:
#         # lag=-1
#         return ferro,EVecValsCombined,sMeanAbsVecCombined,lag,0.05
#     else:
#         eps=0.05
#         acfOfEVec=np.abs(sm.tsa.acf(EVecValsCombined))
#         acfOfSVec=np.abs(sm.tsa.acf(sMeanAbsVecCombined))
#         # print("ferro="+str(ferro))
#         # print(np.min(acfOfEVec))
#         # print(np.min(acfOfSVec))
#         # print("================")
#         minAcfEAll=np.min(acfOfEVec)
#         minAcfSAll=np.min(acfOfSVec)
#         # print(minAcfEAll)
#         # print(minAcfSAll)
#         if minAcfEAll>eps or minAcfSAll>eps:
#             eps=np.max([minAcfEAll,minAcfSAll])
#         # print(len(EVecValsCombined))
#         # print(eps)
#         lagEVec=np.where(acfOfEVec<=eps)[0][0]
#         lagSVec=np.where(acfOfSVec<=eps)[0][0]
#         lag=np.max([lagEVec,lagEVec])
#         # return ferro,EVecValsCombined,sMeanAbsVecCombined,lag
#
#         return ferro,EVecValsCombined,sMeanAbsVecCombined,lag,eps

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



outRoot=pathGroupRow
def diagnosticsAndObservables(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: diagnostic plots and observable values
    """
    tOneFileStart=datetime.now()
    #diagnostics
    ferro,EVecValsCombined,sMeanAbsVecCombined,lag,fileNumSelected=combineValues(oneTFile)#lagVal(oneTFile)

    TTmpMatch=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    matchPartNum=re.search(r"part(\d+)",oneTFile)
    if matchPartNum:
        partNum=int(matchPartNum.group(1))
    # EHistAllDir=outRoot+"part"+str(partNum)+"EHistAll"
    # sHistAllDir=outRoot+"part"+str(partNum)+"sHistAll"
    # EBlkMeanDir=outRoot+"part"+str(partNum)+"EBlkMean"
    # Path(EHistAllDir).mkdir(parents=True, exist_ok=True)
    # Path(sHistAllDir).mkdir(parents=True, exist_ok=True)
    # Path(EBlkMeanDir).mkdir(parents=True,exist_ok=True)
    if TTmpMatch:
        TTmp=float(TTmpMatch.group(1))
    ##############ferro case################################################################
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
        # plt.savefig(EHistAllDir+"/"+EHistOut)
        plt.close()

        plt.figure()
        plt.hist(sMeanAbsVecCombined,bins=nbins)
        meanS=np.mean(sMeanAbsVecCombined)
        sdS=np.sqrt(np.var(sMeanAbsVecCombined)/len(sMeanAbsVecCombined))
        plt.title("ferromagnetic, T="+str(TTmp)+", mean="+str(meanS)+", sd="+str(sdS),fontsize=10)
        plt.xlabel("$|s|$")
        sHistOut="T"+str(TTmp)+"sHist.png"
        plt.savefig(oneTFile+"/"+sHistOut)
        # plt.savefig(sHistAllDir+"/"+sHistOut)
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
        # plt.savefig(EBlkMeanDir+"/T"+str(TTmp)+"EBlk.png")
        plt.close()
        #observables
        chi_ps,hfInterval=JackknifeForChi(sMeanAbsVecCombined,TTmp)
        chiOutFileName="T"+str(TTmp)+"chi.txt"
        contents=["chi="+str(chi_ps)+"\n","hfLength="+str(hfInterval)+"\n"+"lag="+str(lag)\
                  +"\nlastElemNum="+str(numElem)]
        fptr1=open(oneTFile+"/"+chiOutFileName,"w+")
        fptr1.writelines(contents)
        fptr1.close()

        chiAllDir=outRoot+"chiAll"
        Path(chiAllDir).mkdir(exist_ok=True,parents=True)
        fptr2=open(chiAllDir+"/"+chiOutFileName,"w+")
        fptr2.writelines(contents)
        fptr2.close()
        #<|s|>, Jackknife
        ps_s=np.mean(sMeanAbsVecCombined)
        var_s=np.var(sMeanAbsVecCombined,ddof=1)
        hfLengh_s=1.96*np.sqrt(var_s/len(sMeanAbsVecCombined))
        outSFileName="T"+str(TTmp)+"s.txt"
        contents=["s="+str(ps_s)+"\n","hfLength="+str(hfInterval)+"\n"+"lag="+str(lag) \
                  +"\nlastElemNum="+str(numElem)]
        fptr3=open(oneTFile+"/"+outSFileName,"w+")
        fptr3.writelines(contents)
        fptr3.close()

        sAllDir=outRoot+"sAll"
        Path(sAllDir).mkdir(exist_ok=True,parents=True)
        fptr4=open(sAllDir+"/"+outSFileName,"w+")
        fptr4.writelines(contents)
        fptr4.close()


    ##############paramagnetic case################################################################
    else:
        # print("entering else")
        EPerSupercell=np.array(EVecValsCombined)/M



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
        fig.suptitle("T="+str(TTmp))
        plt.savefig(oneTFile+"/T"+str(TTmp)+"EBlk.png")
        # plt.savefig(EBlkMeanDir+"/T"+str(TTmp)+"EBlk.png")
        plt.close()







        # print("finish block mean")
        halfLength=int(len(EPerSupercell)/2)
        #
        # #histogram of distribution of epsilon
        EVecPart0=EPerSupercell[:halfLength]
        EvecPart1=EPerSupercell[halfLength:]



        sVecPart0=sMeanAbsVecCombined[:halfLength]
        sVecPart1=sMeanAbsVecCombined[halfLength:]

        ESelectedFromPart0=EVecPart0[::lag]
        ESelectedFromPart1=EvecPart1[::lag]
        # print("total E val numbers="+str(len(EPerSupercell)))
        # print("lag="+str(lag))
        # print("E val numbers selected from part0="+str(len(ESelectedFromPart0)))
        # print("E val numbers selected from part1="+str(len(ESelectedFromPart1)))
        # EComTmp=np.r_[ESelectedFromPart0,ESelectedFromPart1]
        # i1=1000
        # i2=14062
        # print(str(i1)+" th vec="+str(EComTmp[i1]*M))
        # print(str(i2)+" th vec="+str(EComTmp[i2]*M))

        #diagnostics of E
        meanEPart0=np.mean(ESelectedFromPart0)
        varEPart0=np.var(ESelectedFromPart0,ddof=1)
        hfLengthEPart0=1.96*np.sqrt(varEPart0/len(ESelectedFromPart0))

        meanEPart1=np.mean(ESelectedFromPart1)
        varEPart1=np.var(ESelectedFromPart1,ddof=1)
        hfLengthEPart1=1.96*np.sqrt(varEPart1/len(ESelectedFromPart1))

        nbins=100

        #histogram of E's part0 and E's part1
        fig=plt.figure()
        axE0=fig.add_subplot(1,2,1)
        (n0,_,_)= axE0.hist(ESelectedFromPart0,bins=nbins)
        meanEPart0=np.round(meanEPart0,4)
        hfLengthEPart0=np.round(hfLengthEPart0,4)
        axE0.set_title("part0, T="+str(np.round(TTmp,3)))
        axE0.set_xlabel("$\epsilon$")
        axE0.set_ylabel("#")
        xPosE0Text=(np.max(ESelectedFromPart0)-np.min(ESelectedFromPart1))*1/2+np.min(ESelectedFromPart0)
        yPosE0Text=np.max(n0)*2/3
        axE0.text(xPosE0Text,yPosE0Text,"mean="+str(meanEPart0)+"\nhalfLength="+str(hfLengthEPart0)+"\nlag="+str(lag))

        axE1=fig.add_subplot(1,2,2)
        (n1,_,_)=axE1.hist(ESelectedFromPart1,bins=nbins)
        meanEPart1=np.round(meanEPart1,4)
        hfLengthEPart1=np.round(hfLengthEPart1,4)
        axE1.set_title("part1, T="+str(np.round(TTmp,3)))
        axE1.set_xlabel("$\epsilon$")
        axE1.set_ylabel("#")
        xPosE1Text=(np.max(ESelectedFromPart1)-np.min(ESelectedFromPart1))*1/2+np.min(ESelectedFromPart1)
        yPosE1Text=np.max(n1)*2/3
        axE1.text(xPosE1Text,yPosE1Text,"mean="+str(meanEPart1)+"\nhalfLength="+str(hfLengthEPart1)+"\nlag="+str(lag))
        EHistOut="T"+str(TTmp)+"EHist.png"

        plt.savefig(oneTFile+"/"+EHistOut)
        # plt.savefig(EHistAllDir+"/"+EHistOut)
        plt.close()
        # print("diagnostics of E finished")
        # #diagnostics of s
        sSelectedFromPart0=sVecPart0[::lag]
        sSelectedFromPart1=sVecPart1[::lag]

        meanSPart0=np.mean(sSelectedFromPart0)
        varSPart0=np.var(sSelectedFromPart0,ddof=1)
        hfLengthSPart0=1.96*np.sqrt(varSPart0/len(sSelectedFromPart0))

        meanSPart1=np.mean(sSelectedFromPart1)
        varSPart1=np.var(sSelectedFromPart1,ddof=1)
        hfLengthSPart1=1.96*np.sqrt(varSPart1/len(sSelectedFromPart1))

        #histogram of s's part0 and s's part1
        fig=plt.figure()
        axS0=fig.add_subplot(1,2,1)
        (n0,_,_)=axS0.hist(sSelectedFromPart0,bins=nbins)
        meanSPart0=np.round(meanSPart0,4)
        hfLengthSPart0=np.round(hfLengthSPart0,4)
        axS0.set_title("part0, T="+str(np.round(TTmp,3)))
        axS0.set_xlabel("$|s|$")
        axS0.set_ylabel("#")
        xPosS0Text=(np.max(sSelectedFromPart0)-np.min(sSelectedFromPart0))*1/2+np.min(sSelectedFromPart0)
        yPosS0Text=np.max(n0)*2/3
        axS0.text(xPosS0Text,yPosS0Text,"mean="+str(meanSPart0)+"\nhalfLength="+str(hfLengthSPart0)+"\nlag="+str(lag))

        axS1=fig.add_subplot(1,2,2)
        (n1,_,_)=axS1.hist(sSelectedFromPart1,bins=nbins)
        meanSPart1=np.round(meanSPart1,4)
        hfLengthSPart1=np.round(hfLengthSPart1,4)
        axS1.set_title("part1, T="+str(np.round(TTmp,3)))
        axS1.set_xlabel("$|s|$")
        axS1.set_ylabel("#")
        xPosS1Text=(np.max(sSelectedFromPart1)-np.min(sSelectedFromPart1))*1/2+np.min(sSelectedFromPart1)
        yPosS1Text=np.max(n1)*2/3
        axS1.text(xPosS1Text,yPosS1Text,"mean="+str(meanSPart1)+"\nhalfLength="+str(hfLengthSPart1)+"\nlag="+str(lag))
        sHistOut="T"+str(TTmp)+"sHist.png"
        plt.savefig(oneTFile+"/"+sHistOut)
        # plt.savefig(sHistAllDir+"/"+sHistOut)
        plt.close()
        # print("diagnostics of s finished")
        #observavles

        #chi
        SVec=np.r_[sSelectedFromPart0,sSelectedFromPart1]
    #     # print(SVec)
    #     # print(len(SVec))
    #     i1=1000
    #     i2=14062
    #     print(str(i1)+" th vec="+str(SVec[i1]))
    #     print(str(i2)+" th vec="+str(SVec[i2]))
    #
    #
        # print("len(SVec)="+str(len(SVec)))
        chi_ps,hfInterval=JackknifeForChi(SVec,TTmp)
    #     # print(chi_ps)
    #     # print(hfInterval)
    #
        chiOutFileName="T"+str(TTmp)+"chi.txt"
        contents=["chi="+str(chi_ps)+"\n","hfLength="+str(hfInterval)+"\n"+"lag="+str(lag)\
                  +"\nlastFilesNum="+str(fileNumSelected)]
        fptr1=open(oneTFile+"/"+chiOutFileName,"w+")
        fptr1.writelines(contents)
        fptr1.close()

        chiAllDir=outRoot+"chiAll"
        Path(chiAllDir).mkdir(exist_ok=True,parents=True)
        fptr2=open(chiAllDir+"/"+chiOutFileName,"w+")
        fptr2.writelines(contents)
        fptr2.close()
        # print("chi finished")
        #<|s|>, Jackknife
        ps_s=np.mean(sMeanAbsVecCombined)
        var_s=np.var(sMeanAbsVecCombined,ddof=1)
        hfLengh_s=1.96*np.sqrt(var_s/len(sMeanAbsVecCombined))
        outSFileName="T"+str(TTmp)+"s.txt"
        contents=["s="+str(ps_s)+"\n","hfLength="+str(hfInterval)+"\n"+"lag="+str(lag) \
                  +"\nlastFilesNum="+str(fileNumSelected)]
        fptr3=open(oneTFile+"/"+outSFileName,"w+")
        fptr3.writelines(contents)
        fptr3.close()

        sAllDir=outRoot+"sAll"
        Path(sAllDir).mkdir(exist_ok=True,parents=True)
        fptr4=open(sAllDir+"/"+outSFileName,"w+")
        fptr4.writelines(contents)
        fptr4.close()
    tOneFileEnd=datetime.now()
    print("one file time: ",tOneFileEnd-tOneFileStart)







# diagnosticsAndObservables(inTFileNamesSorted[0])

tStart=datetime.now()
# # procNum=48
# # #parallel
# # # pool0=Pool(procNum)
# # #
# # # ret=pool0.map(diagnosticsAndObservables,inTFileNamesSorted)
# # #serial
for file in inTFileNamesSorted:
    diagnosticsAndObservables(file)
tEnd=datetime.now()
print("total time: ",tEnd-tStart)

