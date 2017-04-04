module load Langs/Intel/15
module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15

mpiexec -n 8 parallel < data/actualdata1 > ./parallel_data1
mpiexec -n 8 parallel < data/actualdata2 > ./parallel_data2
mpiexec -n 8 parallel < data/actualdata3 > ./parallel_data3
mpiexec -n 8 parallel < data/actualdata4 > ./parallel_data4

./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata1 > ./serial_data1.txt
./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata2 > ./serial_data2.txt
./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata3 > ./serial_data3.txt
./serial < /home/fas/hpcprog/ahs3/cpsc424/assignment4/data/actualdata4 > ./serial_data4.txt

