#include "circuits.h"


void
NoBootsXOR(LweSample *result, const LweSample *ca, const LweSample *cb, const TFheGateBootstrappingCloudKeySet *bk) {
    //static const Torus32 MU = modSwitchToTorus32(1,8);
    const LweParams *in_out_params = bk->params->in_out_params;

    //result = new_LweSample(in_out_params);

    //compute: (0,1/4) + 2*(ca + cb)
    //static const Torus32 XorConst=modSwitchToTorus32(1,4);
    //lweNoiselessTrivial(temp_result, XorConst, in_out_params);
    lweAddMulTo(result, 2, ca, in_out_params);
    lweAddMulTo(result, 2, cb, in_out_params);
    //bootsCOPY(result,temp_result,bk);
    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    //tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);

    //delete_LweSample(temp_result);
}

void
NoBootsAND(LweSample *result, const LweSample *ca, const LweSample *cb, const TFheGateBootstrappingCloudKeySet *bk) {
    //static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;

    //LweSample* temp_result = new_LweSample(in_out_params);

    //compute: (0,-1/8) + ca + cb
    //static const Torus32 AndConst=modSwitchToTorus32(-1,8);
    //lweNoiselessTrivial(temp_result, AndConst, in_out_params);
    lweAddTo(result, ca, in_out_params);
    lweAddTo(result, cb, in_out_params);

    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    //tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);

    //delete_LweSample(temp_result);
}
void add(LweSample* result, const LweSample* a, const LweSample* b, const TFheGateBootstrappingCloudKeySet* bk,const LweSample* carry) {
    LweSample* t1 = new_gate_bootstrapping_ciphertext_array(1, bk->params);
    bootsXOR(t1,a,carry,bk);

    LweSample* t2 = new_gate_bootstrapping_ciphertext_array(1, bk->params);
    bootsXOR(t2,b,carry,bk);

    bootsXOR(&result[0], a, t2, bk);//not bootstrapped

    bootsAND(t1,t1,t2,bk);

    //LweSample* res1 = new_gate_bootstrapping_ciphertext_array(1, bk->params);
    bootsXOR(&result[1],carry,t1,bk);
    //bootsCOPY(&result[1],res1,bk);

    delete_gate_bootstrapping_ciphertext_array(1,t1);
    delete_gate_bootstrapping_ciphertext_array(1,t2);
}
void addOneMore(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk,const LweSample* carry) {
    LweSample* new_res = new_gate_bootstrapping_ciphertext_array(2, bk->params);
    add(new_res,&a[0],&b[0],bk,carry);
    bootsCOPY(&result[0],&new_res[0],bk);

    //LweSample* t = new_gate_bootstrapping_ciphertext_array(2, bk->params);
    for (unsigned i=0; i<nb_bits-2; i++) {//2 cause of BITS+1
        add(new_res,&a[i+1],&b[i+1],bk,&new_res[1]);
        bootsCOPY(&result[i+1],&new_res[0],bk);
    }

    bootsCOPY(&result[nb_bits-1],&new_res[1],bk);

    //printf("%d\n",bootsSymDecrypt(&result[BITS],bk))
    //delete_gate_bootstrapping_ciphertext_array(2,t);
    delete_gate_bootstrapping_ciphertext_array(2,new_res);
}

int decryptCheck(LweSample *answer, int nb_bits, TFheGateBootstrappingSecretKeySet *key) {
    int int_answer = 0;
    for (int i = 0; i < nb_bits; i++) {
        int ai = bootsSymDecrypt(&answer[i], key);
        //printf("%d,",ai);
        int_answer |= (ai << i);
    }
    if (int_answer > pow(2, nb_bits) / 2) {
        int_answer = -1 * (pow(2, nb_bits) - int_answer);

    }
    return int_answer;
    //delete_gate_bootstrapping_ciphertext_array(BITS, answer);
    //delete_gate_bootstrapping_secret_keyset(key);
}

void addOneMoreThreaded(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk,const LweSample* carry) {
    LweSample* new_res = new_gate_bootstrapping_ciphertext_array(2, bk->params);
    add(new_res,&a[0],&b[0],bk,carry);
    bootsCOPY(&result[0],&new_res[0],bk);

    //LweSample* t = new_gate_bootstrapping_ciphertext_array(2, bk->params);
    #pragma omp parallel num_threads(4)
        #pragma omp parallel for
        for (unsigned i=0; i<nb_bits-2; i++) {//2 cause of BITS+1
            add(new_res,&a[i+1],&b[i+1],bk,&new_res[1]);
            bootsCOPY(&result[i+1],&new_res[0],bk);
        }

    bootsCOPY(&result[nb_bits-1],&new_res[1],bk);

    //printf("%d\n",bootsSymDecrypt(&result[BITS],bk))
    //delete_gate_bootstrapping_ciphertext_array(2,t);
    delete_gate_bootstrapping_ciphertext_array(2,new_res);
}
void add(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk,const LweSample* carry) {
    LweSample* new_res = new_gate_bootstrapping_ciphertext_array(2, bk->params);
    add(new_res,&a[0],&b[0],bk,carry);
    bootsCOPY(&result[0],&new_res[0],bk);

    //LweSample* t = new_gate_bootstrapping_ciphertext_array(2, bk->params);
    for (unsigned i=0; i<nb_bits-1; i++) {
        add(new_res,&a[i+1],&b[i+1],bk,&new_res[1]);
        bootsCOPY(&result[i+1],&new_res[0],bk);
    }

    bootsCOPY(&result[nb_bits-1],&new_res[1],bk);
    //printf("%d\n",bootsSymDecrypt(&result[nb_bits],bk))
    //delete_gate_bootstrapping_ciphertext_array(2,t);
    delete_gate_bootstrapping_ciphertext_array(2,new_res);

}
void add(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) {
    //no carry for this operation, so carry = 0
    LweSample* carry = new_gate_bootstrapping_ciphertext_array(1, bk->params);
    bootsCONSTANT(&carry[0], 0, bk);
    //printf("%d\n",nb_bits);
    add(result, a, b, nb_bits, bk, carry);
}

void sub(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) {

    LweSample* notb = new_gate_bootstrapping_ciphertext_array(nb_bits-1, bk->params);
    for (int i=0; i<nb_bits-1; i++) bootsNOT(&notb[i], &b[i], bk);

    LweSample* carry = new_gate_bootstrapping_ciphertext_array(1, bk->params);
    bootsCONSTANT(carry, 1, bk);
    //printf("%d\n",nb_bits);
    addOneMore(result, a, notb, nb_bits, bk,carry);

    //bootsCONSTANT(&result[nb_bits-1],0, bk);
    //bootsCOPY(result, carry,bk);
    delete_gate_bootstrapping_ciphertext_array(nb_bits-1,notb);
    delete_gate_bootstrapping_ciphertext_array(1,carry);
}

void bootsMUXGC(LweSample* result, const LweSample* x, const LweSample* y, const LweSample* c, const TFheGateBootstrappingCloudKeySet* bk){
    bootsMUX(result,c,x,y,bk);
}

void twosComplement(LweSample* result, const LweSample* a,  const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) {
    LweSample* reachOne = new_gate_bootstrapping_ciphertext_array(1, bk->params);
    bootsCONSTANT(&reachOne[0], 0, bk);
    for (int i=0; i<nb_bits; i++){
        bootsXOR(&result[i],&a[i],&reachOne[0],bk);
        bootsOR(&reachOne[0],&reachOne[0],&a[i],bk);
    }
}

void zeros(LweSample* result, int nb_bits, const TFheGateBootstrappingCloudKeySet* bk){
for (int i=0; i<nb_bits; i++){
        bootsCONSTANT(&result[i], 0, bk);
    }
}

void transformCipherText(LweSample *result, LweSample *a, int size_a, int bitsize,
                         const TFheGateBootstrappingCloudKeySet *bk) {
    //LweSample* zero= new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);
    //zero = zeros(&zeros,nb_bits,&bk);
//    /std::cout<<"total bits"<<nb_bits<<std::endl;
    LweSample *tmp = new_gate_bootstrapping_ciphertext_array(bitsize, bk->params);
    zeros(tmp, bitsize, bk);
    for (int j = size_a - 1; j < bitsize; j++)
        bootsCOPY(&tmp[j], &a[size_a - 1], bk);

    for (int j = 0; j < size_a; j++)
        bootsCOPY(&tmp[j], &a[j], bk);

    copyVar(result, tmp, bitsize, bk);

    delete_gate_bootstrapping_ciphertext_array(bitsize, tmp);
    //std::cout<<"LS"<<j<<std::endl;

}
void copyVar(LweSample* result, LweSample* input,int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) {
    for (int i=0; i<nb_bits; i++){
        bootsCOPY(&result[i], &input[i], bk);
    }
}
void mulBinary(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk){

    for (int i=0; i<nb_bits; i++){
        bootsAND(&result[i],&a[0],&b[i],bk);
    }
}
void mul(LweSample* result, const LweSample* a, const LweSample* b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) {

    LweSample* toAdd =new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);
    LweSample* zeros = new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);
    for (int i=0; i<nb_bits; i++){
        bootsCONSTANT(&zeros[i], 0, bk);
        bootsCONSTANT(&toAdd[i], 0, bk);
    }
    //LweSample* toAdd = zeros;
    //LweSample* toAdd =zeros;
    for (int i=0; i<nb_bits; i++){
        bootsMUX(&toAdd[i],&b[0],&a[i],&zeros[i],bk);
    
        bootsCOPY(&result[i],&toAdd[i],bk);
    }
    for (int i=1; i<nb_bits; ++i) {
        
        for (int j=i; j<i+nb_bits-1; j++) 
            bootsCOPY(&toAdd[j-i],&result[j],bk);
        
        LweSample* tmp = new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);
        for (int j=0; j<nb_bits; j++)
            bootsMUX(&tmp[j],&b[i], &a[j],&zeros[j],bk);
        //LweSample* res = new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);
        add(toAdd,tmp,toAdd,nb_bits, bk);

        for (int j=0; j<nb_bits; j++) 
            bootsCOPY(&result[i+j],&toAdd[j],bk);
        delete_gate_bootstrapping_ciphertext_array(nb_bits,tmp);
        //delete_gate_bootstrapping_ciphertext_array(nb_bits,res);
    }
    /*delete_gate_bootstrapping_ciphertext_array(nb_bits,zeros);
    delete_gate_bootstrapping_ciphertext_array(nb_bits,toAdd);*/
}

void mul(LweSample* result,const LweSample* a,const LweSample* b,const int nb_bits_a,const int nb_bits_b, const TFheGateBootstrappingCloudKeySet* bk) {

    int nb_bits = nb_bits_a+nb_bits_b;

    LweSample* toAdd =new_gate_bootstrapping_ciphertext_array(nb_bits_a, bk->params);
    LweSample* zeros = new_gate_bootstrapping_ciphertext_array(nb_bits_a, bk->params);

    for (int i=0; i<nb_bits_a; i++){
        bootsCONSTANT(&zeros[i], 0, bk);
        bootsCONSTANT(&toAdd[i], 0, bk);
    }

    //LweSample* toAdd = zeros;
    //LweSample* toAdd =zeros;
    for (int i=0; i<nb_bits_a; i++){
        bootsMUX(&toAdd[i],&b[0],&a[i],&zeros[i],bk);

        bootsCOPY(&result[i],&toAdd[i],bk);
    }

    for (int i=1; i<nb_bits_b; ++i) {

        for (int j=i; j<i+nb_bits_a; j++)
            bootsCOPY(&toAdd[j-i],&result[j],bk);

        LweSample* tmp = new_gate_bootstrapping_ciphertext_array(nb_bits_a, bk->params);

        for (int j=0; j<nb_bits_a; j++)
            bootsMUX(&tmp[j],&b[i], &a[j],&zeros[j],bk);
        //LweSample* res = new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);

        add(toAdd,tmp,toAdd,nb_bits_a, bk);

        for (int j=0; j<nb_bits_a; j++)
            bootsCOPY(&result[i+j],&toAdd[j],bk);

        delete_gate_bootstrapping_ciphertext_array(nb_bits_a,tmp);
        //delete_gate_bootstrapping_ciphertext_array(nb_bits,res);
    }


    delete_gate_bootstrapping_ciphertext_array(nb_bits_a,zeros);
    //delete_gate_bootstrapping_ciphertext_array(nb_bits_a,toAdd);
}
//right shift by b bits
void rightShift(LweSample* result, const LweSample* a, const int b, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) {
    //LweSample* zero= new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);
    //zero = zeros(&zeros,nb_bits,&bk);
    LweSample *tmp = new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);
    zeros(tmp, nb_bits, bk);

    for (int j=0; j<nb_bits-b; j++)
        bootsCOPY(&tmp[j], &a[j + b], bk);
    for (int j = nb_bits - b; j < nb_bits; j++)//for signed
        bootsCOPY(&tmp[j], &a[nb_bits - 1], bk);
    //for (int j=nb_bits-b; j<nb_bits; j++)
    //bootsCOPY(&tmp[nb_bits-1],&a[nb_bits-1],bk);
    copyVar(result, tmp, nb_bits, bk);

    delete_gate_bootstrapping_ciphertext_array(nb_bits, tmp);
}
//left shift by b bits
//expands the ciphertext
void leftShift(LweSample* result, const LweSample* a, const int bitshifts, const int nb_bits, const TFheGateBootstrappingCloudKeySet* bk) {

    LweSample *tmp = new_gate_bootstrapping_ciphertext_array(nb_bits, bk->params);
    zeros(tmp, nb_bits, bk);

    for (int j = 0; j < nb_bits - bitshifts; j++)
        bootsCOPY(&tmp[j + bitshifts], &a[j], bk);

    copyVar(result, tmp, nb_bits, bk);

    delete_gate_bootstrapping_ciphertext_array(nb_bits, tmp);
}
