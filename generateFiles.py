import numpy as np

#this script generates MCMC computing scripts

TAll=np.linspace(0.1,3.1,70)

counter=0
for T in TAll:
    bashContents = []
    bashContents.append("#!/bin/bash\n")
    bashContents.append("#SBATCH -n 12\n")
    bashContents.append("#SBATCH -N 1\n")
    bashContents.append("#SBATCH -t 0-40:00\n")
    bashContents.append("#SBATCH -p CLUSTER\n")
    bashContents.append("#SBATCH --mem=140GB\n")
    bashContents.append("#SBATCH -o outmc" + str(counter) + ".out\n")
    bashContents.append("#SBATCH -e outmc" + str(counter) + ".err\n")
    bashContents.append("cd /home/cywanag/liuxi/Documents/cppCode/paramagnet_interaction_cpp\n")
    command="./mc "+str(T)
    bashContents.append(command)
    bsFileName="./mcBash/mcsh"+str(counter)+".sh"
    fptrTmp=open(bsFileName,"w+")
    for oneline in bashContents:
        fptrTmp.write(oneline)
    fptrTmp.close()
    counter+=1