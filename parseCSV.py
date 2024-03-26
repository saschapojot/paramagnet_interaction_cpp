import pandas as pd
import sys

#this script parses csv containing parameters L,M,t,J,g

#python parseCSV.py groupNum rowNum
if len(sys.argv)!=3:
    print("wrong number of arguments")

groupNum=int(sys.argv[1])
rowNum=int(sys.argv[2])

inParamFileName="./group"+str(groupNum)+".csv"

df=pd.read_csv(inParamFileName)
oneRow=df.iloc[rowNum,:]

L=int(oneRow.loc["L"])
M=int(oneRow.loc["M"])
t=float(oneRow.loc["t"])
J=float(oneRow.loc["J"])
g=float(oneRow.loc["g"])
f=float(oneRow.loc["f"])
print("L"+str(L)+"M"+str(M)+"t"+str(t)+"J"+str(J)+"g"+str(g)+"f"+str(f))