#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8 
#$ -q comp.q
#$ -P prabha_group
#$ -V
./clean.sh
ip=`/sbin/ifconfig eth0 | awk '/inet addr/ {gsub("addr:", "", $2); print $2}' | cut -d'.' -f4`
echo "node$((ip-2))" > ifconfig.txt
ipd="19"
ipr="2"
dirname=`basename $PWD`
actdir=$PWD
mkdir -p /scratch/dilip/$actdir
python genPlumedDeca_opt.py $3 $1 $2 $4 $5
if [ $ip != $ipr ]
then
   cd "/scratch/dilip/$actdir"
   mpirun -np $1 mdrun_mpi -s $actdir/filesImplicit2/mxtc_topol -o /scratch/$actdir/traj  -e ./ener  -plumed $actdir/plumed.dat -multi $1 -nsteps $2
else
   mpirun -np $1 mdrun_mpi -s filesImplicit2/mxtc_topol -o ./traj  -e ./ener  -plumed -multi $1 -nsteps $2
fi
if [[ $ip != $ipd && $ip != $ipr ]];
then
    scp -r /scratch/dilip/$actdir $USER@172.31.1.19:/scratch/dilip/$actdir/../
    rm -r /scratch/dilip/$actdir
fi

