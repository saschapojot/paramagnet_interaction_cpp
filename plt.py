import matplotlib.pyplot as plt
import numpy
import numpy as np
import pandas as pd



part=0

inDir="./part"+str(part)+"/"

inFileName=inDir+"reduced.csv"


dtFrm=pd.read_csv(inFileName)

def str2float(vec):
    ret=[float(elem) for elem in vec]
    return vec

TAll=np.array(dtFrm["T"])

sAll=np.array(dtFrm["s"])
chiAll=np.array(dtFrm["chi"])
CAll=np.array(dtFrm["C"])

#plot <s> vs T
fig,ax=plt.subplots()
ax.scatter(TAll,sAll,color="black")
plt.xlabel("$T$")
plt.ylabel("<s>")
plt.title("Temperature from "+str(TAll[0])+" to "+str(TAll[-1]))
# ax.spines['left'].set_position('zero')
# ax.spines['bottom'].set_position('zero')
# # Hide top and right spines
# ax.spines['right'].set_color('none')
# ax.spines['top'].set_color('none')
plt.yticks([0,0.2,0.4,0.6,0.8,1])
ax.tick_params(axis='both', which='major', labelsize=6)
plt.savefig(inDir+"T"+str(TAll[0])+"toT"+str(TAll[-1])+"sAvg.png")
plt.close()

# plot chi vs T

fig,ax=plt.subplots()
ax.scatter(TAll,chiAll,color="red")
plt.title("Temperature from "+str(TAll[0])+" to "+str(TAll[-1]))
plt.xlabel("$T$")
plt.ylabel("$\chi$")
# ax.spines['left'].set_position('zero')
# ax.spines['bottom'].set_position('zero')
# # Hide top and right spines
# ax.spines['right'].set_color('none')
# ax.spines['top'].set_color('none')
ax.tick_params(axis='both', which='major', labelsize=6)
plt.savefig(inDir+"T"+str(TAll[0])+"toT"+str(TAll[-1])+"Chi.png")
plt.close()


#plot C vs T
fig,ax=plt.subplots()
ax.scatter(TAll,CAll,color="blue")
plt.title("Temperature from "+str(TAll[0])+" to "+str(TAll[-1]))
plt.xlabel("$T$")
plt.ylabel("$C$")
# ax.spines['left'].set_position('zero')
# ax.spines['bottom'].set_position('zero')
#
# # Hide top and right spines
# ax.spines['right'].set_color('none')
# ax.spines['top'].set_color('none')
ax.tick_params(axis='both', which='major', labelsize=6)
plt.savefig(inDir+"T"+str(TAll[0])+"toT"+str(TAll[-1])+"specificHeat.png")

