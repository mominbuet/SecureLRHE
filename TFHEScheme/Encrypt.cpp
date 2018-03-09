
#include <chrono>
#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include "tfhe/tfhe.h"
#include "tfhe/tfhe_io.h"
#include "circuits.h"


using namespace std;
int BITS = 1;



void LoadFile(char *filename, int m, int n, int **rawData, int *labels) {
    int y = 0;
    //bool ret[][] = new bool[m][n];
    // int ret[m][n];

    std::string output;
    std::ifstream in(filename);
    if (in.is_open()) {
        getline(in, output);//skipping first line
        while (!in.eof()) {

            //vector<int> data(m-1);
            //getline(in, output);
            /*char *dup = strdup(output.c_str());
            char * token = strtok(dup, ",");*/
            getline(in, output);
            std::stringstream ss(output);
            for (int x = 0; x < m; x++) {
                getline(ss, output, ',');
                int tmp = atoi(output.c_str());
                //std::cout<<output<<",";
                if (x == 0)
                    labels[y] = tmp;
                else {
                    rawData[y][x - 1] = tmp;
                }
            }
            y++;
            if (y >= n)
                return;

        }
    } else {
        printf("cannot open file\n");
        return;
    }
    printf("read file, rows %d columns %d\n",n,m );
    in.close();
    return;
}


int main(int argc, char *argv[]) {

    char *datafile;
    int m, n, sec = 80;


    if (argc >= 4) {
        datafile = argv[1];
        m = atoi(argv[2]);
        n = atoi(argv[3]);
        sec = atoi(argv[4]);
        //cout << "Starting " << datafile << " covariates " << m << " data " << n << endl;
    } else {
        printf("usage: Filename covariatecount datacount");
        return 1;
    }
    int **rawData = new int*[n];
    for (int i = 0; i < n; i++)
        rawData[i] = new int[m - 1];

    int *labels = new int[n];
    LoadFile(datafile,m,n,rawData,labels);


    const int minimum_lambda = sec;//80 bit security
    TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(minimum_lambda);

    //generate a random key
    uint32_t seed[] = { 314, 1592, 657 };
    tfhe_random_generator_setSeed(seed,3);
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);

    //export the secret key to file for later use
    FILE* secret_key = fopen("secret.key","wb");
    export_tfheGateBootstrappingSecretKeySet_toFile(secret_key, key);
    fclose(secret_key);

    //export the cloud key to a file (for the cloud)
    FILE* cloud_key = fopen("cloud.key","wb");
    export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
    fclose(cloud_key);
    FILE* cloud_data = fopen("cloud_matrix.data","wb");
    printf("Saving data to file!\n");

    for (int y = 0; y < n; y++) {
        //std::cout <<endl<< labels[y] << ",";
        LweSample* ciphertext1 = new_gate_bootstrapping_ciphertext_array(BITS, params);
        for (int i=0; i<BITS; i++) {
            bootsSymEncrypt(&ciphertext1[i], (labels[y]>>i)&1, key);
            export_gate_bootstrapping_ciphertext_toFile(cloud_data, &ciphertext1[i], params);
        }
        for (int x = 0; x < m - 1; x++) {
            //std::cout << rawData[y][x] << ",";
            LweSample*  ciphertext1 = new_gate_bootstrapping_ciphertext_array(BITS, params);
            for (int i=0; i<BITS; i++) {
                bootsSymEncrypt(&ciphertext1[i], (rawData[y][x] >> i) & 1, key);
                export_gate_bootstrapping_ciphertext_toFile(cloud_data, &ciphertext1[i], params);
            }
            delete_gate_bootstrapping_ciphertext_array(BITS, ciphertext1);
        }
        //std::cout << std::endl;

    }
    fclose(cloud_data);


    FILE *cloud_data_lables = fopen("cloud_labels.data", "wb");
    int FIXEDBITSIZE = 8 * BITS;
    LweSample *carry = new_gate_bootstrapping_ciphertext_array(BITS, params);
    bootsCONSTANT(&carry[0], 0, &key->cloud);

    LweSample *five100 = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE,
                                                                 params);//fixed 8,change later!! (2days later) after spending 1.5 days on this later :(
    zeros(five100, FIXEDBITSIZE, &key->cloud);
    for (int l = 0; l < FIXEDBITSIZE; l++) {//00110010
        bootsCONSTANT(&five100[l],(l==1||l==4||l==5) ?1:0,&key->cloud);
        //bootsCONSTANT(&five100[l], (l == 2 || (l > 3 && l < 9)) ? 1 : 0, &key->cloud);//500
    }

    LweSample *five100C = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, params);
    //leftShift(thousand,five100,1,FIXEDBITSIZE,&key->cloud);


    twosComplement(five100C, five100, FIXEDBITSIZE, &key->cloud);

    for (int y = 0; y < n; y++) {
        LweSample *ciphertext1 = new_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, params);
        zeros(ciphertext1, FIXEDBITSIZE, &key->cloud);
        //LweSample* tmp1 = new_gate_bootstrapping_ciphertext_array(1+FIXEDBITSIZE, params);
        //zeros(tmp1,1+FIXEDBITSIZE,&key->cloud);

        if (labels[y] == 1)
            copyVar(ciphertext1, five100, FIXEDBITSIZE, &key->cloud);
        else
            copyVar(ciphertext1, five100C, FIXEDBITSIZE, &key->cloud);
        //add(tmp1,ciphertext1,five100, FIXEDBITSIZE,&key->cloud);
        for (int i = 0; i < FIXEDBITSIZE; i++) {
            export_gate_bootstrapping_ciphertext_toFile(cloud_data_lables, &ciphertext1[i], params);
        }

        //std::cout<<decryptCheck(ciphertext1, FIXEDBITSIZE, key)<<std::endl;
        delete_gate_bootstrapping_ciphertext_array(FIXEDBITSIZE, ciphertext1);
        //delete_gate_bootstrapping_ciphertext_array(1+FIXEDBITSIZE, tmp1);
    }
    fclose(cloud_data_lables);

    for (int i = 0; i < n; ++i) {
        delete[] rawData[i];
    }

    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);
    printf ("Generated Keys!\n");
}