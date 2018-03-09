#!/bin/bash
#$1 =data file
#$2 =number of covariates +1
#$3 =number of rows
#$4 = iteration (2/3 iterations)
#$5 = test data file
#$6 = record count from $5
#$7=batchsize
ulimit -s unlimited
export OMP_STACKSIZE=8G
let bit_security=80
batchsize=$7
iteration=$4
#if [ $batchsize -lt $[$3/$4] ]; then
#	batchsize=$[$3/$4]
#fi
#if [ $batchsize -gt 200 ]; then
#	batchsize=200
#fi
#iteration=$[$3/$batchsize]
iteration=$4
#if [ $iteration -gt 2 ]; then #change later for higher iterations
#	iteration=2
#fi
echo "Batchsize $batchsize iterations allowed $iteration"
rm -rf cloud_matrix.data cloud_labels.data cloud.key secret.key your-program

g++ Encrypt.cpp circuits.cpp  -o your-program -ltfhe-spqlios-avx -std=c++11
SECONDS=0
./your-program $1 $2 $3 $bit_security
duration=$SECONDS
echo -e "\033[0;31mEncryption time $(($duration / 60)) minutes and $(($duration % 60)) seconds\033[0m"

FILENAME=cloud_labels.data
FILESIZE=$(stat -c%s "$FILENAME")
FILENAME=cloud_matrix.data
FILESIZE=$[$FILESIZE+$(stat -c%s "$FILENAME")]
echo -e "\033[0;31mSize of encrypted data = $(($FILESIZE / 1024)) kilobytes\033[0m"


g++ CloudThread.cpp circuits.cpp  -o your-program -ltfhe-spqlios-avx  -fopenmp -std=c++11
SECONDS=0
/usr/bin/time -v ./your-program $2 $3 $batchsize $iteration $5 $6
duration=$SECONDS
echo -e "\033[0;31mTraining time $(($duration / 60)) minutes and $(($duration % 60)) seconds\033[0m"
#/usr/bin/time -v

g++  Decrypt.cpp circuits.cpp -o your-program -ltfhe-spqlios-avx -std=c++11
SECONDS=0
./your-program $5 beta.data $2 $6
Rscript AUC.R
duration=$SECONDS
echo -e "\033[0;31mDecryption time $(($duration / 60)) minutes and $(($duration % 60)) seconds\033[0m"
rm -rf cloud_matrix.data cloud_labels.data #cloud.key secret.key your-program
