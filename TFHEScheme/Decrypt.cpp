
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
int BITS = 15;

void LoadFile(char *filename, int m, int n, int **rawData, int *labels) {
    int y = 0;
    std::string output;
    std::ifstream in(filename);
    if (in.is_open()) {
        getline(in, output);//skipping first line
        while (!in.eof()) {
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
                //ret[x][y]=atoi(output.c_str());
                //cout<<ret[x][y]<<endl;
            }

            //free(dup);
            //getline(in, output);
            //std::cout<<std::endl;
            y++;
            if (y >= n)
                return;

            //cout << atoi(output)<<endl;
        }
    } else {
        printf("cannot open file\n");
        return;
    }
    //cout << "read file, rows" << rawData.NumRows() << " columns  " << rawData.NumCols() << endl;
    printf("read file, rows %d columns %d\n", n, m);
    in.close();
    return;
    //fclose(infile);
    //return ret;
}


int main(int argc, char *argv[]) {
    char *betafile, *datafile;
    int m, n;
    if (argc >= 2) {
        datafile = argv[1];
        betafile = argv[2];
        m = atoi(argv[3]);
        n = atoi(argv[4]);
        //cout << "Starting " << datafile << " covariates " << m << " data " << n << endl;
    } else {
        printf("usage: TestDataFilename BetaFilename covariatecount records_count");
        return 1;
    }
    FILE *cloud_key = fopen("cloud.key", "rb");
    TFheGateBootstrappingCloudKeySet *bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
    fclose(cloud_key);

    FILE *secret_key = fopen("secret.key", "rb");
    TFheGateBootstrappingSecretKeySet *secretkey = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
    fclose(secret_key);

    FILE *beta_data = fopen(betafile, "rb");

    std::cout << "beta " << std::endl;
    double beta[m-1];
    for (int k = 0; k < m-1; ++k) {
        LweSample *tmp = new_gate_bootstrapping_ciphertext_array(BITS, bk->params);
        zeros(tmp, BITS, bk);

        for (int i = 0; i < BITS; i++) //from Encrypt.cpp
            import_gate_bootstrapping_ciphertext_fromFile(beta_data, &tmp[i], bk->params);


        //std::cout <<"labels "<<y<<" "<<decryptCheck(labelsMinusFile[y],1+FIXEDBITSIZE,secretkey)<<std::endl;

        beta[k] = decryptCheck(tmp, BITS, secretkey)/100.0;
        std::cout << beta[k] << ",";
        delete_gate_bootstrapping_ciphertext_array(BITS, tmp);
    }
    cout << endl;
    fclose(beta_data);

    /**
    * Test script
     */
    int **rawData = new int *[n];
    for (int i = 0; i < n; i++)
        rawData[i] = new int[m - 1];

    int *labels = new int[n];
    LoadFile(datafile, m, n, rawData, labels);
    double res[n];
    //double mean = 0.0;
    double min=0.0,max=0.0;
   
    for (unsigned y = 0; y < n; y++) {
        res[y] = 0.0;
        for (unsigned x = 0; x < m-1; x++) {
            res[y] += rawData[y][x] * beta[x];
        }
         if(res[y]>max)
            max=res[y];
        if(res[y]<min)
            min = res[y];
        //mean += res[y];
    }
    //mean = mean / n;
    //cout<<mean<<endl;
    ofstream myoutfile;
    myoutfile.open("output.txt");
    for (unsigned y = 0; y < n; y++) {
        //res[y] = res[y] >= 0 ? 1 : 0;
        //cout<<res[y]<<",";
        res[y]=(res[y]-min)/(max-min);
        myoutfile << res[y] << ((y == n - 1) ? "" : ",");
    }
    myoutfile << "\n";
    for (unsigned y = 0; y < n; y++) {
        myoutfile << labels[y] << ((y == n - 1) ? "" : ",");
    }
    myoutfile << "\n";
    //cout << endl;
    myoutfile.close();

    return 0;
}