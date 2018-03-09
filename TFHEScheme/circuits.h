#ifndef CIRCUITS_H
#define CIRCUITS_H
#include <omp.h>
#include <stdio.h>
#include "tfhe/tfhe.h"
#include "tfhe/tfhe_io.h"


void add(LweSample* result, const LweSample* a, const LweSample* b, const TFheGateBootstrappingCloudKeySet* bk,const LweSample* carry);

void addOneMore(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk,const LweSample* carry);//the numberof bits are one more than the original one, regular addtion
void addOneMoreThreaded(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk,const LweSample* carry);//the numberof bits are one more than the original one, regular addtion

void add(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk,const LweSample* carry);

void add(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk);//the numberof bits are equal to the input, will not result the same sometimes

void sub(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk);

void bootsMUXGC(LweSample* result, const LweSample* x, const LweSample* y, const LweSample* c, const TFheGateBootstrappingCloudKeySet* bk);

void zeros(LweSample* result, int nb_bits, const TFheGateBootstrappingCloudKeySet* bk);

void copyVar(LweSample* result, LweSample* input,int nb_bits, const TFheGateBootstrappingCloudKeySet* bk);

void mul(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk);

void mulBinary(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk);

void mul(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits_a,const  int nb_bits_b, const TFheGateBootstrappingCloudKeySet* bk) ;
//right shift by b bits
void rightShift(LweSample* result, const LweSample* a, const int b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) ;

void leftShift(LweSample* result, const LweSample* a, const int bitshifts, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) ;

void twosComplement(LweSample* result, const LweSample* a,  const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk);

void transformCipherText(LweSample *result, LweSample *a, int size_a, int bitsize,
                         const TFheGateBootstrappingCloudKeySet *bk);
int decryptCheck(LweSample *answer, int nb_bits, TFheGateBootstrappingSecretKeySet *key);
#endif