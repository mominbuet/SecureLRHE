#!/bin/bash

#$1 =data file
#$2 =number of covariates +1
#$3 =number of rows
#$4 = iteration (1/2 iterations)
#$5 = test data file
#$6 = record count from $5


rm -rf beta.encrypted data.encrypted labels.encrypted
RED='\033[0;31m'
NC='\033[0m'
ulimit -s unlimited
export OMP_STACKSIZE=8G

#let batchsize=$[$3/$4]-1
batchsize=20
if [ $batchsize -gt $[$[$2-1]*4] ]; then
	batchsize=$[$[$2-1]*4]
fi
if [ $[$3/$4] -gt $batchsize ]; then
	batchsize=$[$3/$4]
fi 


echo -e "\033[0;31mbatchsize $batchsize iteration $4\033[0m"
# Encryption part
g++ Encrypt.cpp -o main.o a.out -std=c++11  -I./include/ -I./include/nfl/ -I./include/nfl/prng/ -I./lib/prng/ -I./lib/params/ -I./include/nfl/opt/arch/ -lgmpxx -lgmp  -lmpfr -m64 -DNTT_AVX -DNTT_SSE
SECONDS=0
./main.o $1 $2 $3
duration=$SECONDS
echo -e "\033[0;31mEncryption time $(($duration / 60)) minutes and $(($duration % 60)) seconds\033[0m"

g++ Cloud.cpp -o main.o a.out -std=c++11  -I./include/ -I./include/nfl/ -I./include/nfl/prng/ -I./lib/prng/ -I./lib/params/ -I./include/nfl/opt/arch/ -lgmpxx -lgmp  -lmpfr -m64 -DNTT_AVX -DNTT_SSE -fopenmp
SECONDS=0
/usr/bin/time -v  ./main.o $2 $3 $batchsize $4
#g++ fv_lr.cpp -o main.o asm.o  -I./include/ -I./include/nfl/ -I./include/nfl/prng/ -I./lib/prng/ -I./lib/params/ -I./include/nfl/opt/arch/ -lgmpxx -lgmp  -lmpfr -m64 -DNTT_AVX -DNTT_SSE && /usr/bin/time -v ./main.o $1 $2 $3 $4 $5
duration=$SECONDS
echo -e "\033[0;31m Execution time $(($duration / 60)) minutes and $(($duration % 60)) seconds \033[0m"

g++ Decrypt.cpp -o main.o a.out -std=c++11  -I./include/ -I./include/nfl/ -I./include/nfl/prng/ -I./lib/prng/ -I./lib/params/ -I./include/nfl/opt/arch/ -lgmpxx -lgmp  -lmpfr -m64 -DNTT_AVX -DNTT_SSE
SECONDS=0
 ./main.o $5 beta.encrypted $2 $6 $4
Rscript AUC.R
duration=$SECONDS
echo -e "\033[0;31mDecryption time$(($duration / 60)) minutes and $(($duration % 60)) seconds \033[0m"
rm -rf beta.encrypted data.encrypted labels.encrypted
