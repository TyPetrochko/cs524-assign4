#!/bin/bash
#PBS -l procs=1,tpn=1,mem=34gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

# cd $PBS_O_WORKDIR

module load Langs/Intel/15

# ./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata1 > ./serial_data1.txt
# ./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata2 > ./serial_data2.txt
# ./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata3 > ./serial_data3.txt
./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata4 > ./serial_data4.txt


#./fserial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/testdata1 > ./testdata1_f.out
