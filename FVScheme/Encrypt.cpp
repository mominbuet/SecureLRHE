#include <chrono>
#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstddef>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <stdlib.h>
#include <nfl.hpp>
#include "lib/params/params.cpp"
#include <sstream>
#include "lib/prng/fastrandombytes.cpp"
#include "lib/prng/randombytes.cpp"
//#include "nfl/prng/crypto_stream_salsa20.h"

const int MODULOUS = 62 * 20;
const int DEGREE = 1 << 8;
const double SIGMA = 3.0;
const std::string T = "461168601831223323332";

namespace FV {
    namespace params {
        using poly_t = nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS>;
        template <typename T>
        struct plaintextModulus;
        template <>
        struct plaintextModulus<mpz_class> {
            static mpz_class value() {//19
                /*poly_t  p = (poly_t)1;
                return  mpz_class(p.moduli_product());*/
                return mpz_class(T);
            }
            //static mpz_class value(int val) { return mpz_class(val); }
        };

        using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
        using gauss_t = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;
        gauss_t fg_prng_sk(SIGMA, 128, 1 << 14);
        gauss_t fg_prng_evk(SIGMA, 128, 1 << 14);
        gauss_t fg_prng_pk(SIGMA, 128, 1 << 14);
        gauss_t fg_prng_enc(SIGMA, 128, 1 << 14);
    }
}  // namespace FV::params
#include "FV.hpp"
using namespace FV;
using namespace std;


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
                    labels[y]=tmp;
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
    printf("read file, rows %d columns %d\n",n,m );
    in.close();
    return;
    //fclose(infile);
    //return ret;
}


void saveSecretKey(const char* filename, sk_t &sk){
    ofstream myfile;
    myfile.open (filename);
    sk.value.serialize_manually(myfile);
    /*std::array <mpz_t, DEGREE> tmp =  sk.value.poly2mpz();
    for (int i = 0; i < DEGREE; ++i) {
        myfile<<tmp[i]<<((i<DEGREE-1)?",":"");
    }
    myfile<<"\n";*/
    myfile.close();

}

void savePublicKey(const char* filename, pk_t &pk){
    ofstream myfile;
    myfile.open (filename);
    pk.a.serialize_manually(myfile);
    pk.b.serialize_manually(myfile);
    /*std::array <mpz_t, DEGREE> tmp =  pk.a.poly2mpz();
    for (int i = 0; i < DEGREE; ++i) {
        myfile<<tmp[i]<<((i<DEGREE-1)?",":"");
    }*/
    //myfile<<"\n";

    /*tmp =  pk.b.poly2mpz();
    for (int i = 0; i < DEGREE; ++i) {
        myfile<<tmp[i]<<((i<DEGREE-1)?",":"");
    }
    myfile<<"";*/

    myfile.close();

}

void saveEvalKey(const char* filename, evk_t &evk){
    ofstream myfile;
    myfile.open (filename);

    for (unsigned j = 0; j < evk.ell; ++j) {
        evk.values[j][0].serialize_manually(myfile);
        evk.values[j][1].serialize_manually(myfile);
        /*std::array<mpz_t, DEGREE> tmp = evk.values[j][0].poly2mpz();
        for (int i = 0; i < DEGREE; ++i) {
            myfile << tmp[i] << ((i < DEGREE - 1) ? "," : "");
        }
        myfile << "\n";
        tmp = evk.values[j][1].poly2mpz();
        for (int i = 0; i < DEGREE; ++i) {
            myfile << tmp[i] << ((i < DEGREE - 1) ? "," : "");
        }
        myfile << "\n";*/

    }

    //myfile << pk->b_shoup<<"\n";
    /*tmp =  pk.delta.poly2mpz();
    for (int i = 0; i < DEGREE; ++i) {
        myfile<<tmp[i]<<((i<DEGREE-1)?",":"");
    }
    myfile<<"\n";*/
    //myfile << pk->delta_shoup<<"\n";
    myfile.close();

}

nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS>  get_array(string line){
    nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS>  res;
    //unsigned long long int* res = new unsigned long long int[8192];
    line.replace(0,1,"");
    line.replace(line.end()-1,line.end(),"");
    //char* tmp;
    //strcpy(tmp,line);
    char * c = (char*)line.c_str();
    char* pch = strtok (c,",");

    int i=0;
    char * pEnd;
    std::array <mpz_t, DEGREE> arr=res.poly2mpz();
    while (pch != NULL)
    {
        mpz_t &s = arr[i++];
        mpz_set_ui(s,strtoull(pch,&pEnd,10));

        //res[i++] = strtoull(pch,&pEnd,10);
        //cout<<s<<endl;
        pch = strtok (NULL, ",");
        if(i>DEGREE-1) break;
    }
    cout<<i<<endl;

    res.mpz2poly(arr);
    return res;

}

void printEncValPrint(ciphertext_t &ct, sk_t &sk,pk_t &pk){
    mess_t m_dec;

    decrypt(m_dec, sk, pk,ct);
    //std::cout <<m_dec;
    mpz_t tmp;
    FV::params::poly_t pt =(FV::params::poly_t) m_dec.getValue();
    std::array <mpz_t, DEGREE> arr = pt.poly2mpz();
    mpz_t &s=arr[0];
    mpz_init(tmp);
    mpz_cdiv_q_ui(tmp,FV::params::plaintextModulus<mpz_class>::value().get_mpz_t() ,2);
    //std::cout<<s<<std::endl;
    if(mpz_cmp(s , tmp) >0) {
        mpz_sub(tmp, FV::params::plaintextModulus<mpz_class>::value().get_mpz_t() , s);
        //mpz_tdiv_q(tmp,tmp,quantization);
        //std::cout << "beta " << i << ": -" <<tmp << std::endl;
        std::cout <<"-" <<tmp ;
    }else {
        //mpz_tdiv_q(tmp,tmp,quantization);
        //std::cout << "beta " << i << ": " << m_dec << std::endl;
        std::cout << m_dec ;
    }
    mpz_clear(tmp);

    //return;
}
void readPublicKey( pk_t &pk){
    string line;
    ifstream myfile ("public.key");
    if (myfile.is_open())
    {
        int lineno=0;
        while ( getline (myfile,line) )
        {
            //unsigned long long int *res=get_array(line);
            if(!line.empty()) {
                switch (lineno) {
                    case 0:
                        pk.a = get_array(line);
                        pk.a_shoup = nfl::compute_shoup(pk.a);
                        break;
                    case 1:
                        pk.b = get_array(line);
                        pk.b_shoup = nfl::compute_shoup(pk.b);
                        break;
                    /*case 2:
                        pk.delta = get_array(line);
                        pk.delta_shoup = nfl::compute_shoup(pk.delta);
                        break;*/
                }
                lineno++;
            }
        }
        myfile.close();
    }
    else cout << "Unable to open file";
}
int main(int argc, char *argv[]) {
    //std::cout <<"started"<<std::endl;

    sk_t sk;
    string secret_key_file = "secret.key";
    saveSecretKey(secret_key_file.c_str(),sk);
    string eval_key_file = "eval.key";

    int word_size = 1 << 6;
    evk_t evk(sk, word_size);
    saveEvalKey(eval_key_file.c_str(),evk);

    pk_t pk(sk, evk);
    //pk_t pk;
    string public_key_file = "public.key";
    savePublicKey(public_key_file.c_str(),pk);
    //readPublicKey(pk);

    std::cout <<"Keygeneration done"<<std::endl;

    char *datafile;
    int m, n;
    if (argc >= 4) {
        datafile = argv[1];
        m = atoi(argv[2]);
        n = atoi(argv[3]);
        //cout << "Starting " << datafile << " covariates " << m << " data " << n << endl;
    } else {
        printf("usage: Filename covariatecount(+1) datacount \n");
        return 1;
    }
    int **rawData = new int*[n];
    for(int i = 0; i < n; ++i)
        rawData[i] = new int[m-1];

    int *labels = new int[n];



    LoadFile(datafile,m,n,rawData,labels);

    //FV::params::poly_p quantization = getPolyPredefined(1000);


    ofstream labelsEncrypted,dataEncrypted;
    labelsEncrypted.open ("labels.encrypted");
    dataEncrypted.open("data.encrypted");


    ciphertext_t rawDataPoly[n][m - 1];
    ciphertext_t labelsPoly[n];
    //poly_t rawDataPolyT[n][m - 1];
    for (unsigned y = 0; y < n; y++) {
        for (unsigned x = 0; x < m - 1; x++) {

            //std::cout << rawData[y][x];

            FV::encrypt_integer(rawDataPoly[y][x],pk, rawData[y][x]);

            rawDataPoly[y][x].c0.serialize_manually(dataEncrypted);
            rawDataPoly[y][x].c1.serialize_manually(dataEncrypted);

        }

        FV::encrypt_integer(labelsPoly[y], pk, labels[y]);
        //std::cout << labels[y]<<std::endl;
        labelsPoly[y].c0.serialize_manually(labelsEncrypted);
        labelsPoly[y].c1.serialize_manually(labelsEncrypted);

    }
    std::cout<<"log2(q)="<<labelsPoly[0].c0.bits_in_moduli_product()<<std::endl;
    std::cout<<"degree="<<labelsPoly[0].c0.degree<<std::endl;
    std::cout<<"q="<<labelsPoly[0].c0.moduli_product()<<std::endl;
    labelsEncrypted.close();dataEncrypted.close();
/*
    cout << "Decrypting records" << endl;

    for (unsigned y = 0; y < n; y++) {
        for (unsigned x = 0; x < m - 1; x++) {
            printEncValPrint(rawDataPoly[y ][x], sk, pk);
        }
        printEncValPrint(labelsPoly[y ], sk, pk);
        cout << endl;
    }*/
    cout << "Encryption done" << endl;
    return  1;
}