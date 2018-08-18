#!/bin/bash 
#$ -S /bin/bash

#$ -l h_rss=1000M,h_fsize=1000M,h_cpu=3:00:00,hw=x86_64
#$ -cwd 

WORKDIR=$HOME/Documents/code/noisedynamics/rungslcode
TARGETDIR=/data2/finite/wildan/noisedyn/datachapter4

EXEC=./a.out
PARAMETERS=(0 0 0 0 0 0 0 0 0 0)

echo ${WORKDIR}
hostname

##################################################################################
### setup simulation directory
if [ ! -d "/scratch/${USER}" ] 
then
        mkdir /scratch/${USER}
fi

for i in $(seq 1 1 500)
do
    if [ ! -d "/scratch/$USER/sim${i}" ]
    then
        mkdir /scratch/${USER}/sim${i}
        while [ $? -ne 0 ]; do
	    let i=i+1
	    mkdir /scratch/${USER}/sim${i}
	done
	break
    fi
done
echo ${i}
##################################################################################
### start program
cd /scratch/${USER}/sim${i}                      # change to scratch directory
${WORKDIR}/${EXEC} ${PARAMETERS[*]}
##################################################################################
### cleanup
if [ $? -eq 0 ]
then
    cp *.*   ${TARGETDIR}
    rm -f *
    cd ..
    rmdir sim${i}
else
    echo error
fi
##################################################################################
