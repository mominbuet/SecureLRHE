
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
#include <omp.h>
#include <nfl.hpp>
#include "lib/params/params.cpp"
#include "lib/prng/fastrandombytes.cpp"
#include "lib/prng/randombytes.cpp"
//#include "nfl/prng/crypto_stream_salsa20.h"
using namespace std;

int CPU =4;

#include "fvnamespace.h"
#include "FV_CW.hpp"

using namespace FV;

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
    // Seed (for deterministic values)
    srand(0);

    // Keygen
    FV::sk_t sk;
    int word_size = 1 << 6;
    FV::evk_t evk(sk, word_size);
    FV::pk_t pk(sk, evk);
    string secret_key_file = "secret.key";
    saveSecretKey(secret_key_file.c_str(),sk);
    string eval_key_file = "eval.key";

    //evk_t evk(sk, word_size);
    saveEvalKey(eval_key_file.c_str(),evk);

    //pk_t pk(sk, evk);
    //pk_t pk;
    string public_key_file = "public.key";
    savePublicKey(public_key_file.c_str(),pk);
    //readPublicKey(pk);

    std::cout <<"Keygeneration done"<<std::endl;

    char *datafile;
    int m, n,batchsize;
    if (argc >= 4) {
        datafile = argv[1];
        m = atoi(argv[2]);
        n = atoi(argv[3]);
        batchsize = atoi(argv[4]);
        //cout << "Starting " << datafile << " covariates " << m << " data " << n << " batchsize "<<batchsize<<endl;
    } else {
        printf("usage: Filename covariatecount(+1) datacount \n");
        return 1;
    }
    int **rawData = new int*[n];
    for(int i = 0; i < n; ++i)
        rawData[i] = new int[m-1];

    int *labels = new int[n];



    LoadFile(datafile,m,n,rawData,labels);

    // Polynomials
    //FV::params::poly_pl m[n/batchsize];

    ciphertext_t rawDataPoly[batchsize][m - 1];
    ciphertext_t labelsPoly[batchsize];
    ofstream labelsEncrypted,dataEncrypted;
    labelsEncrypted.open ("labels_packed.encrypted");
    dataEncrypted.open("data_packed.encrypted");

    for (size_t k = 0; k <batchsize ; ++k) {//10

        for (size_t l = 0; l <m-1 ; l++) {
            std::array<mpz_t, DEGREE> ntt_polym0;
            for (size_t j = 0; j < n / batchsize; j++) {
                //cout<< k+j*batchsize<<",";
                mpz_init_set_ui(ntt_polym0[j], rawData[k+j*batchsize][l]);
            }
            //cout<<"-"<<l<<endl;
            for (size_t j = n/batchsize ;j<DEGREE; j++) {
                mpz_init_set_ui(ntt_polym0[j], 0);
            }
            FV::params::poly_pl m;
            m.mpz2poly(ntt_polym0);
            FV::encrypt_poly(rawDataPoly[k][l], pk, m);
            rawDataPoly[k][l].c0.serialize_manually(dataEncrypted);
            rawDataPoly[k][l].c1.serialize_manually(dataEncrypted);
            //clear
            for (size_t i = 0; i < DEGREE; i++) {
                mpz_clears(ntt_polym0[i], nullptr);
            }

            //for decrypt/check
            /*std::array<mpz_t, DEGREE> polym0;
            for (size_t i = 0; i < DEGREE; i++) {
                mpz_inits(polym0[i],  nullptr);
            }
            FV::decrypt_poly(polym0, sk, pk, rawDataPoly[k][l], true);
            for (size_t i = 0; i < n / batchsize; i++) {
                std::cout << mpz_class(polym0[i]).get_str()<<" ";
            }
            cout<<endl;
            for (size_t i = 0; i < DEGREE; i++) {
                mpz_clears(polym0[i], nullptr);
            }
            */
        }

        std::array<mpz_t, DEGREE> ntt_polym0;
        for (size_t j = 0; j < n / batchsize; j++) {
            //cout<< k+j*batchsize<<",";
            mpz_init_set_ui(ntt_polym0[j], labels[k+j*batchsize]);
        }
        //cout<<endl;
        for (size_t j = n/ batchsize ;j<DEGREE; j++) {
            mpz_init_set_ui(ntt_polym0[j], 0);
        }
        FV::params::poly_pl m;
        m.mpz2poly(ntt_polym0);
        FV::encrypt_poly(labelsPoly[k], pk, m);
        labelsPoly[k].c0.serialize_manually(labelsEncrypted);
        labelsPoly[k].c1.serialize_manually(labelsEncrypted);
        //clear
        for (size_t i = 0; i < DEGREE; i++) {
            mpz_clears(ntt_polym0[i], nullptr);
        }

    }
    labelsEncrypted.close();dataEncrypted.close();
    cout << "Packed Encryption done" << endl;
    /*for (size_t i = 0; i < DEGREE; i++) {
        mpz_init_set_si(ntt_polym0[i],i+1);
    }

    m[0]={1,2,3,4,5,6,7,8};//.mpz2poly(ntt_polym0);// = {1, 2, 3, 4, 5, 6, 7, 8};

    for (size_t i=0;i<8;i++)
        mpz_set_si(ntt_polym0[i],8-i);
    m[1].mpz2poly(ntt_polym0);// = {8, 7, 6, 5, 4, 3, 2, 1};

    for (size_t i = 0; i < 8; i++) {
        mpz_clear(ntt_polym0[i]);
    }

    // Encrypt
    std::array<FV::ciphertext_t, 2> c;
    FV::encrypt_poly(c[0], public_key, m[0]);
    FV::encrypt_poly(c[1], public_key, m[1]);

    // Initialize polym
    std::array<mpz_t, 8> polym0, polym1;
    for (size_t i = 0; i < 8; i++) {
        mpz_inits(polym0[i], polym1[i], nullptr);
    }

    // decrypt to the polym
    FV::decrypt_poly(polym0, secret_key, public_key, c[0], true);
    FV::decrypt_poly(polym1, secret_key, public_key, c[1], true);

    // Script sage
    std::cout << "# Sage script for the verification" << std::endl;
    std::cout << "p=" << FV::params::plaintextModulus<mpz_class>::value().get_str() << std::endl;
    std::cout << "K.<X> = QuotientRing(Integers(p)[x], Integers(p)[x].ideal(x^8 + 1));"
              << std::endl;

    std::cout << "m0 = ";
    for (size_t i = 0; i < 8; i++) {
        std::cout << mpz_class(polym0[i]).get_str() << "*X^" << i
                  << (i == 7 ? ";\n" : "+");
    }
    std::cout << "m1 = ";
    for (size_t i = 0; i < 8; i++) {
        std::cout << mpz_class(polym1[i]).get_str() << "*X^" << i
                  << (i == 7 ? ";\n" : "+");
    }

    // Multiplication and decryption
    FV::ciphertext_t prod = c[0] * c[1];
    FV::decrypt_poly(polym0, secret_key, public_key, prod, true);

    std::cout << "m = ";
    for (size_t i = 0; i < 8; i++) {
        std::cout << mpz_class(polym0[i]).get_str() << "*X^" << i
                  << (i == 7 ? ";\n" : "+");
    }

    std::cout << "m == m0*m1" << std::endl;

    // Clean
    for (size_t i = 0; i < 8; i++) {
        mpz_clears(polym0[i], polym1[i], nullptr);
    }*/

    return 0;
}
