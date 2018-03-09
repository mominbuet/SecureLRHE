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
        template<typename T>
        struct plaintextModulus;

        template<>
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

void readPublicKey(pk_t &pk) {

    ifstream myfile("public.key");
    if (myfile.is_open()) {
        pk.a.deserialize_manually(myfile);
        pk.a_shoup = nfl::compute_shoup(pk.a);
        pk.b.deserialize_manually(myfile);
        pk.b_shoup = nfl::compute_shoup(pk.b);
        /*
        int lineno=0;
        while ( getline (myfile,line) )
        {
            //unsigned long long int *res=get_array(line);
            if (!line.empty()) {
                switch(lineno){
                    case 0:pk.a = get_array(line);break;
                    case 1: pk.b = get_array(line);pk.b_shoup = nfl::compute_shoup(pk.b);break;
                    //case 2: pk.delta =get_array(line);pk.delta_shoup = nfl::compute_shoup(pk.delta);break;
                }
                lineno++;
            }
        }*/
        myfile.close();
    } else cout << "Unable to open file";
}

void readSecretKey(sk_t &sk) {
    ifstream myfile("secret.key");
    if (myfile.is_open()) {
        sk.value.deserialize_manually(myfile);
        sk.value_shoup = nfl::compute_shoup(sk.value);
        /*int lineno=0;
        while ( getline (myfile,line) ) {
            //unsigned long long int *res=get_array(line);
            if (!line.empty()) {
                sk.value = get_array(line);
                sk.value_shoup = nfl::compute_shoup(sk.value);

                lineno++;
                break;
            }
        }*/

        myfile.close();
    } else cout << "Unable to open file";
}

void readEvalKey(evk_t &evk) {
    ifstream myevalfile("eval.key");

    if (myevalfile.is_open()) {
        for (unsigned j = 0; j < evk.ell; ++j) {
            evk.values[j] = new FV::params::poly_p[2];
            evk.values_shoup[j] = new FV::params::poly_p[2];
            evk.values[j][0].deserialize_manually(myevalfile);
            evk.values_shoup[j][0] = nfl::compute_shoup(evk.values[j][0]);
            evk.values[j][1].deserialize_manually(myevalfile);
            evk.values_shoup[j][1] = nfl::compute_shoup(evk.values[j][1]);
        }
        /*int lineno=0;
        while ( getline (myevalfile,line) ) {
            if (!line.empty()) {

                if(lineno%2==0) {
                    evk.values[lineno / 2] = new FV::params::poly_p[2];evk.values_shoup[lineno/2]= new FV::params::poly_p[2];
                }
                evk.values[lineno/2][lineno%2] = get_array(line);
                evk.values_shoup[lineno/2][lineno%2] = nfl::compute_shoup(evk.values[lineno/2][lineno%2]);
                //cout<<lineno<<endl;

                lineno++;

            }else break;
        }*/

        myevalfile.close();
    } else cout << "Unable to eval open file";
}

int main(int argc, char *argv[]) {


    char *betafile, *datafile;
    int m, n, iter;
    if (argc >= 2) {
        datafile = argv[1];
        betafile = argv[2];
        m = atoi(argv[3]);
        n = atoi(argv[4]);
        iter = atoi(argv[5]);
        //cout << "Starting " << datafile << " covariates " << m << " data " << n << endl;
    } else {
        printf("usage: TestDataFilename BetaFilename covariatecount records_count");
        return 1;
    }
    sk_t sk;
    readSecretKey(sk);
    evk_t evk(1 << 6);


    readEvalKey(evk);
    pk_t pk1;
    readPublicKey(pk1);
    //pk.evk=&evk;
    pk_t pk = pk_t(pk1.a, pk1.a_shoup, pk1.b, pk1.b_shoup, evk);
    std::cout << "Read all the keys" << std::endl;

    double quant = 10000.0;
    if (iter > 1) {
        for (int i = 1; i < iter; ++i) {
            quant = quant * quant * 100*100;
        }
    }
    ifstream myfile(betafile);
    double beta[m-1];
    mess_t m_dec;
    ciphertext_t ct;
    cout << "beta :";
    for (int i = 0; i < m-1; ++i) {
        ct.c0.deserialize_manually(myfile);
        ct.c1.deserialize_manually(myfile);
        ct.isnull = false;
        ct.pk = &pk;

        decrypt(m_dec, sk, pk, ct);
        mpz_t tmp;
        FV::params::poly_t pt = (FV::params::poly_t) m_dec.getValue();
        std::array<mpz_t, DEGREE> arr = pt.poly2mpz();
        mpz_t &s = arr[0];
        mpz_init(tmp);
        mpz_cdiv_q_ui(tmp, FV::params::plaintextModulus<mpz_class>::value().get_mpz_t(), 2);
        //std::cout<<s<<std::endl;
        if (mpz_cmp(s, tmp) > 0) {
            mpz_sub(tmp, FV::params::plaintextModulus<mpz_class>::value().get_mpz_t(), s);
            //mpz_tdiv_q(tmp,tmp,quantization);
            //std::cout << "beta " << i << ": -" <<tmp << std::endl;
            beta[i] = mpz_get_ui(tmp) / (quant * (-1));

        } else {
            //mpz_tdiv_q(tmp,tmp,quantization);
            //std::cout << "beta " << i << ": " << m_dec << std::endl;
            beta[i] = mpz_get_ui(s) / (quant);

        }
        cout << beta[i];
        mpz_clear(tmp);
        if (i < m - 2)
            cout << ",";
    }
    cout << endl;
    

    myfile.close();

    /**
     * Test script
     */
    int **rawData = new int *[n];
    for (int i = 0; i < n; i++)
        rawData[i] = new int[m - 1];

    int *labels = new int[n];
    LoadFile(datafile, m, n, rawData, labels);
    double res[n];
    double mean = 0.0,min=0.0,max=0.0;
    for (unsigned y = 0; y < n; y++) {
        res[y] = 0.0;
        for (unsigned x = 0; x < m-1; x++) {
            res[y] += rawData[y][x] * beta[x];
        }
        mean += res[y];
        if(res[y]>max)
            max=res[y];
        if(res[y]<min)
            min = res[y];
    }
    mean = mean / n;
    //cout<<avg<<endl;
    ofstream myoutfile;
    myoutfile.open("output.txt");
    for (unsigned y = 0; y < n; y++) {
        //res[y] = res[y] >=avg ? 1 : 0;
        res[y]=(res[y]-min)/(max-min);
        //cout<<res[y]<<",";
        myoutfile << res[y] << ((y == n - 1) ? "" : ",");
    }
    myoutfile << "\n";
    for (unsigned y = 0; y < n; y++) {
        myoutfile << labels[y] << ((y == n - 1) ? "" : ",");
    }
    myoutfile << "\n";
    cout << endl;
    myoutfile.close();

    return 1;
}