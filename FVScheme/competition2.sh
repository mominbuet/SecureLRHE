#!/bin/bash

#$1 =data file
#$2 =number of covariates +1
#$3 =number of rows
#$4 = iteration (1/2 iterations)
#$5 = test data file
#$6 = record count from $5
#$7 = batchsize
#./competition2.sh 9_training_data.csv 20 1421 2 9_testing_data.csv 100 200
rm -rf beta_packed.encrypted data_packed.encrypted labels_packed.encrypted main.o
RED='\033[0;31m'
NC='\033[0m'
#let batchsize=$7
#max(n/10, min(20, m*4))
ulimit -s unlimited
export OMP_STACKSIZE=8G
batchsize=$7
#if [ $batchsize -gt $[$[$2-1]*4] ]; then
	#batchsize=$[$[$2-1]*4]
#fi
#if [ $[$3/10] -gt $batchsize ]; then
#	batchsize=$[$3/10]
#fi

let b=$3/$batchsize

echo -e "\033[0;31mbatchsize $batchsize packing $b covariates $[$2-1] records $3\033[0m"
## Encryption part
g++ EncryptPack.cpp -o main.o a.out  -I./include/ -I./include/nfl/ -I./include/nfl/prng/ -I./lib/prng/ -I./lib/params/ -I./include/nfl/opt/arch/ -std=c++11 -lgmpxx -lgmp  -lmpfr -m64 -DNTT_AVX -DNTT_SSE 
SECONDS=0
./main.o $1 $2 $3 $batchsize
duration=$SECONDS
echo -e "\033[0;31mEncryption time $(($duration / 60)) minutes and $(($duration % 60)) seconds\033[0m"


FILENAME=data_packed.encrypted
FILESIZE=$(stat -c%s "$FILENAME")
FILENAME=labels_packed.encrypted
FILESIZE=$[$FILESIZE+$(stat -c%s "$FILENAME")]
echo -e "\033[0;31mSize of encrypted data = $(($FILESIZE / (1024*1024))) MB\033[0m"


g++ CloudFusion2.cpp -o cloud.o a.out  -I./include/ -I./include/nfl/ -I./include/nfl/prng/ -I./lib/prng/ -I./lib/params/ -I./include/nfl/opt/arch/ -std=c++11 -lgmpxx -lgmp  -lmpfr -m64 -DNTT_AVX -DNTT_SSE  -fopenmp
SECONDS=0
/usr/bin/time -v ./cloud.o $2 $3 $batchsize $4
duration=$SECONDS
echo -e "\033[0;31mTraining time $(($duration / 60)) minutes and $(($duration % 60)) seconds\033[0m"

g++ DecryptPack.cpp -o main.o a.out  -I./include/ -I./include/nfl/ -I./include/nfl/prng/ -I./lib/prng/ -I./lib/params/ -I./include/nfl/opt/arch/ -std=c++11  -lgmpxx -lgmp  -lmpfr -m64 -DNTT_AVX -DNTT_SSE

SECONDS=0
./main.o $5 beta_packed.encrypted $2 $6 $4 $b
Rscript AUC.R
duration=$SECONDS
echo -e "\033[0;31mDecryption time $(($duration / 60)) minutes and $(($duration % 60)) seconds\033[0m"
