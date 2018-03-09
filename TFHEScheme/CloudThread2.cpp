
#include <chrono>
#include <iostream>
#include <bitset>
#include <iterator>
#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <omp.h>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include "circuits.h"

using namespace std;
int INPUTBITS = 1;
int FIXEDBITSIZE = 16;
int CPU =8;
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

void rightshiftTest(const TFheGateBootstrappingCloudKeySet *bk, TFheGateBootstrappingSecretKeySet *key) {
    LweSample *a1 = new_gate_bootstrapping_ciphertext_array(16, bk->params);
    zeros(a1, 16, bk);
    bootsCONSTANT(&a1[12], 1, bk);
    twosComplement(a1, a1, 16, bk);
    std::cout << decryptCheck(a1, 16, key) << " " << 16 << std::endl;
    LweSample *a2 = new_gate_bootstrapping_ciphertext_array(16, bk->params);
    rightShift(a2, a1, 10, 16, bk);

    std::cout << decryptCheck(a2, 16, key) << " " << 16 << std::endl;
}

/**
 * make a to bitsize by padding 0s on the left
 * @param result
 * @param a
 * @param size_a
 * @param bitsize
 * @param bk
 */

LweSample* getInnerProduct( LweSample *data,LweSample *beta, const TFheGateBootstrappingCloudKeySet* bk){
	LweSample*  result = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE,bk->params);	
	mulBinary(result, data, beta, FIXEDBITSIZE, bk);
	return result;
}

int main(int argc, char *argv[]) {
	system("date");
    //char *datafile="cloud_matrix.data";
    int m, n, MAX_ITER = 1;
    int batchsize = 25;
    if (argc >= 4) {
        //datafile = argv[1];
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        batchsize = (n > atoi(argv[3])) ? atoi(argv[3]) : n;
        MAX_ITER = atoi(argv[4]);
        std::cout << "Starting with covariates " << (m-1) << " data " << n << " batchsize " << batchsize
                  << std::endl;
    } else {
        printf("usage: EncryptedFile covariatecount datacount\n");
        return 1;
    }
/*    if(batchsize%2!=0)
    {
    	printf("batchsize needs to be even.\n");
        return 1;
    }
	if(batchsize%CPU!=0)
    {
    	std::cout<<"batchsize needs to be a factor of "<<CPU<<".\n";
        return 1;
    }*/
	
    //int DOUBLEFIXEDBITSIZE = 2 * FIXEDBITSIZE;
    FILE *cloud_key = fopen("cloud.key", "rb");
    TFheGateBootstrappingCloudKeySet *bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
    fclose(cloud_key);

    /**
     * remove this later on
     */
    FILE *secret_key = fopen("secret.key", "rb");
    TFheGateBootstrappingSecretKeySet *secretkey = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
    fclose(secret_key);

    /**
     * to here
     */

    /**
     * TEST block
     *
     */


    //rightshiftTest(bk,secretkey);
    //return 0;
/**
 * init block
 */


    LweSample *beta[m - 1];
    for (int j = 0; j < m - 1; j++) {
        beta[j] = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, bk->params);//sizing with the currentbitsize
        zeros(beta[j], FIXEDBITSIZE, bk);
        bootsCONSTANT(&beta[j][0], 0, bk);
    }
    LweSample *five106 = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE,
                                                                 bk->params);//fixed 8,change later!! (2days later) after spending 1.5 days on this later :(
    zeros(five106, FIXEDBITSIZE, bk);
    for (int l = 0; l < FIXEDBITSIZE; l++) {//00110010
        bootsCONSTANT(&five106[l],(l==1||l==4||l==5) ?1:0,bk);
        //bootsCONSTANT(&five106[l], (l == 5 || l == 8 || l == 13 || (l > 14 && l < 19)) ? 1 : 0, bk);//500
    }
    //std::cout << "five 10^6 " << decryptCheck(five106, FIXEDBITSIZE, secretkey) << std::endl;
    twosComplement(five106, five106, FIXEDBITSIZE, bk);

    LweSample *carry = new_gate_bootstrapping_ciphertext_array(INPUTBITS, bk->params);
    bootsCONSTANT(&carry[0], 0, bk);

    /**
     * computation block
     */
	//std::chrono::milliseconds span (100);
    //const TFheGateBootstrappingParameterSet* params = bk->params;
    //int batchsize = (n>20)?20:n;
    int start = 0;int labelSize = 8;
    FILE *cloud_data = fopen("cloud_matrix.data", "rb");
    FILE *cloud_lables_data = fopen("cloud_labels.data", "rb");
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        std::cout << "started iteration " << iter << " start " << start << std::endl;
        //LweSample *inbetween;//= new_gate_bootstrapping_ciphertext_array(1+FIXEDBITSIZE, bk->params);

        LweSample *data[batchsize][m-1];
        LweSample *labels[batchsize];
        LweSample *ciphertext;
        for (unsigned y = start; y < start + batchsize; y++) {
            for (unsigned x = 0; x < m; x++) {
                ciphertext = new_gate_bootstrapping_ciphertext_array(INPUTBITS, bk->params);
                for (int i = 0; i < INPUTBITS; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext[i], bk->params);

                if (x == m - 1) {
                    labels[y - start] = &ciphertext[0];
                } else {
                    data[y - start][x] = &ciphertext[0];
                }
                //
            }
        }

        LweSample *labelsMinus[batchsize];
        
        for (unsigned y = start; y < start + batchsize; y++) {
            LweSample *tmp = new_gate_bootstrapping_ciphertext_array(labelSize, bk->params);
            //zeros(labelsMinusFile[y],1+FIXEDBITSIZE,bk);
            for (int i = 0; i < labelSize; i++) //from Encrypt.cpp
                import_gate_bootstrapping_ciphertext_fromFile(cloud_lables_data, &tmp[i], bk->params);

            labelsMinus[y - start] = new_gate_bootstrapping_ciphertext_array(1 + FIXEDBITSIZE, bk->params);

            transformCipherText(labelsMinus[y - start], tmp, labelSize, 1 + FIXEDBITSIZE, bk);
            //std::cout <<"labels "<<y<<" "<<decryptCheck(labelsMinus[y-start],1+FIXEDBITSIZE,secretkey)<<std::endl;
            delete_gate_bootstrapping_ciphertext_array(labelSize, tmp);
        }


        //std::cout << "Grabbed Inputs" << std::endl;

        /*for (unsigned y = start; y < start + batchsize; y++) {
            for (unsigned x = 0; x < m - 1; x++) {
                std::cout << decryptCheck(data[y - start][x], INPUTBITS, secretkey);
            }
            std::cout << std::endl;
        }*/
        /**
         * X times beta
         */
        //std::cout << "Calcing Xbeta" << std::endl;
        //LweSample *labelsMinus[batchsize];
		//int y=start,x=0;
		
			//#pragma omp parallel num_threads(CPU)
			{
				if (iter > 0) {
				//#pragma omp parallel for
				for (int  y = start; y < start + batchsize; y++) {

                        LweSample *inbetween = new_gate_bootstrapping_ciphertext_array(1+FIXEDBITSIZE, bk->params);
						zeros(inbetween, 1+FIXEDBITSIZE, bk);
						LweSample *tmp;
						//std::future<LweSample*> fut[m-1] ;
						for (unsigned x = 0; x < (m - 1); x+=1) {
							tmp = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, bk->params);
							zeros(tmp, FIXEDBITSIZE, bk);
							mulBinary(tmp, data[y - start][x], beta[x], FIXEDBITSIZE, bk);

							/*std::cout << "mul " << decryptCheck(data[y - start][x], INPUTBITS, secretkey) << "*"
									  << decryptCheck(beta[x], FIXEDBITSIZE, secretkey) << "="
									  << decryptCheck(tmp, FIXEDBITSIZE, secretkey);*/
							if (x == 0) copyVar(inbetween, tmp, FIXEDBITSIZE, bk);
							else addOneMore(inbetween, inbetween, tmp, 1+FIXEDBITSIZE, bk, carry);

							delete_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, tmp);
							/*for (int i = 0; i < ((m-1)-x<CPU?(m-1)-x:CPU); ++i)
							{
								fut[x+i] = std::async(getInnerProduct,data[y - start][x+i],beta[x+i],bk);
							}
							for (int i = 0; i < ((m-1)-x<CPU?(m-1)-x:CPU); ++i)
							{
								if (x == 0) copyVar(inbetween, fut[x+i].get(), FIXEDBITSIZE, bk);
								else addOneMore(inbetween, inbetween, fut[x+i].get(), 1+FIXEDBITSIZE, bk, carry);
							}*/
							
							
						}
						/*for (int x = 0; x < m-1; x++)
						{
							
							//tmp[x] = fut[x].get();
							if (x == 0) copyVar(inbetween, fut[x].get(), FIXEDBITSIZE, bk);
							else addOneMore(inbetween, inbetween, fut[x].get(), 1+FIXEDBITSIZE, bk, carry);
							//std::cout << " add  " << decryptCheck(inbetween, FIXEDBITSIZE, secretkey) << std::endl;
							//delete_gate_bootstrapping_ciphertext(tmp[x]);
						}*/

						

						/**
						* Sigmoid Approx
						*/

						/*LweSample *square = new_gate_bootstrapping_ciphertext_array(DOUBLEFIXEDBITSIZE, bk->params);//y^2
						zeros(square, DOUBLEFIXEDBITSIZE, bk);
						mul(square, inbetween, inbetween, FIXEDBITSIZE, bk);*/
						//std::cout << "square " << decryptCheck(square, DOUBLEFIXEDBITSIZE, secretkey);


						LweSample *twenty4XB = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, bk->params);//25*1000^2
						//LweSample *twenty4XB2 = new_gate_bootstrapping_ciphertext_array(DOUBLEFIXEDBITSIZE, bk->params);//25*1000^2
						zeros(twenty4XB, FIXEDBITSIZE, bk);
						//zeros(twenty4XB2, DOUBLEFIXEDBITSIZE, bk);
						//copyVar(twenty4XB, inbetween, FIXEDBITSIZE, bk);
						rightShift(twenty4XB, inbetween, 2, FIXEDBITSIZE, bk);//.24~1/4
						//std::cout << ".24XB " << decryptCheck(twenty4XB, FIXEDBITSIZE, secretkey);
						twosComplement(twenty4XB, twenty4XB, FIXEDBITSIZE, bk);


						LweSample *sigmoidResult = new_gate_bootstrapping_ciphertext_array(1 + FIXEDBITSIZE, bk->params);
						bootsCONSTANT(&sigmoidResult[FIXEDBITSIZE],0,bk);
						copyVar(sigmoidResult,twenty4XB,FIXEDBITSIZE,bk);
						//addOneMore(sigmoidResult,square,labelsMinus[y-start],1+DOUBLEFIXEDBITSIZE,bk,carry);
						//addOneMore(sigmoidResult, square, twenty4XB, 1 + DOUBLEFIXEDBITSIZE, bk, carry);

						//addOneMore(sigmoidResult, sigmoidResult, five106, 1 + FIXEDBITSIZE, bk, carry);
						//rightShift(sigmoidResult, sigmoidResult, 10, 1 + FIXEDBITSIZE, bk);

						//std::cout << " sigmoidResult  " << decryptCheck(sigmoidResult, FIXEDBITSIZE, secretkey);


						delete_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, twenty4XB);
						//delete_gate_bootstrapping_ciphertext_array(DOUBLEFIXEDBITSIZE, square);
						/**
						 * Labels + x2-24x-0.5
						 */
						//labelsMinus[y - start] = new_gate_bootstrapping_ciphertext_array(1 + FIXEDBITSIZE, bk->params);//y^2
						//zeros(labelsMinus[y - start], 1 + FIXEDBITSIZE, bk);
						//copyVar(labelsMinus[y - start], labels[y - start], INPUTBITS, bk);//can be cheaper
						//leftShift(labelsMinus[y - start], labelsMinus[y - start], 10, 1 + FIXEDBITSIZE, bk);
						//std::cout << "labels  " << decryptCheck(labelsMinus[y - start], 1 + FIXEDBITSIZE, secretkey);

						addOneMore(labelsMinus[y - start], labelsMinus[y - start], sigmoidResult, 1 + FIXEDBITSIZE, bk, carry);
						//std::cout << " labels minus " << decryptCheck(labelsMinus[y - start], FIXEDBITSIZE, secretkey) << std::endl;
					/*} else {
						//copyVar(labelsMinus[y-start],labelsMinusFile[y-start],1+FIXEDBITSIZE,bk);
					}*/
                    delete_gate_bootstrapping_ciphertext(inbetween);
				}
			}
		
			//#pragma omp parallel for
			for (int y = 0; y < m - 1; y++) {
                LweSample *inbetween = new_gate_bootstrapping_ciphertext_array(1 + FIXEDBITSIZE, bk->params);
				zeros(inbetween, 1 + FIXEDBITSIZE, bk);
				LweSample *tmp;
				for (int  x = start; x < start + batchsize; x++) {
					tmp = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, bk->params);
					zeros(tmp, FIXEDBITSIZE, bk);
					mulBinary(tmp, data[x - start][y], labelsMinus[x - start], FIXEDBITSIZE, bk);
					/*std::cout << "mul " << decryptCheck(data[x - start][y], INPUTBITS, secretkey) << "*"
							  << decryptCheck(labelsMinus[x - start], FIXEDBITSIZE, secretkey) << "="
							  << decryptCheck(tmp, FIXEDBITSIZE, secretkey);*/
					if (x == start) copyVar(inbetween, tmp, FIXEDBITSIZE, bk);
					else addOneMore(inbetween, inbetween, tmp, 1 + FIXEDBITSIZE, bk, carry);
					
					delete_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, tmp);
					//std::cout << " add  " << decryptCheck(inbetween, FIXEDBITSIZE, secretkey) << std::endl;
				}
				//#pragma omp critical
				if (iter == 0)
					copyVar(beta[y], inbetween, FIXEDBITSIZE, bk);
				else {
					leftShift(beta[y], beta[y], 6, FIXEDBITSIZE, bk);//as we did not divide in the last phase
					//std::cout << "before add " << decryptCheck(beta[y], FIXEDBITSIZE - 1, secretkey);
					bootsCONSTANT(&inbetween[FIXEDBITSIZE], 0, bk);

					addOneMore(beta[y], beta[y], inbetween, FIXEDBITSIZE, bk, carry);
					//std::cout << "after add " << decryptCheck(beta[y], FIXEDBITSIZE - 1, secretkey) << std::endl;
				}
				

				rightShift(beta[y], beta[y], 6, FIXEDBITSIZE - 1, bk);//multiply with .01
				//std::cout<< "beta " << y <<":" << decryptCheck(beta[y], FIXEDBITSIZE - 1, secretkey)<< std::endl;
				//copyVar(beta[y],tmp1,FIXEDBITSIZE,bk);
                delete_gate_bootstrapping_ciphertext(inbetween);
			}
		}
        /**
         * batching enabled
         */
        if (start + 2 * batchsize < n)
            start += batchsize;
        else if (start + batchsize == n){
            rewind(cloud_data);
            rewind(cloud_lables_data);
           
            start = 0;
        }
        else
            start = n - (start + batchsize);
        /*for (int y = 0; y < batchsize; y++) {
            delete_gate_bootstrapping_ciphertext_array(1+FIXEDBITSIZE, labelsMinus[y]);
        }*/


        //
        delete_gate_bootstrapping_ciphertext_array(INPUTBITS, ciphertext);
        FILE *beta_data = fopen("beta.data", "wb");
        cout << "Saving beta " << iter << " to file!" << endl;
		//#pragma omp critical
        for (int k = 0; k < m - 1; ++k) {
            LweSample *ciphertext1 = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE - 1, bk->params);
            copyVar(ciphertext1, beta[k], FIXEDBITSIZE - 1, bk);
            for (int i = 0; i < FIXEDBITSIZE - 1; i++) {
                //bootsSymEncrypt(&ciphertext1[i], (rawData[y][x] >> i) & 1, key);
                export_gate_bootstrapping_ciphertext_toFile(beta_data, &ciphertext1[i], bk->params);
            }
            delete_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE - 1, ciphertext1);

        }
        fclose(beta_data);
    }
    /*std::cout << "beta " << std::endl;
    for (int k = 0; k < m - 1; ++k) {
        std::cout << decryptCheck(beta[k], FIXEDBITSIZE - 1, secretkey) << std::endl;
        delete_gate_bootstrapping_ciphertext(beta[k]);
    }*/
    fclose(cloud_data);
    fclose(cloud_lables_data);

    return 0;
}
