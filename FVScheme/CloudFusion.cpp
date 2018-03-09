
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

#include <nfl.hpp>
#include "lib/params/params.cpp"

#include "lib/prng/fastrandombytes.cpp"
#include "lib/prng/randombytes.cpp"
//#include "nfl/prng/crypto_stream_salsa20.h"
using namespace std;


const int MODULOUS = 62 * 5;
const int DEGREE = 1 << 12;
const double SIGMA = 3.0;
const std::string T = "461168601832";
/// include the FV homomorphic encryption library
namespace FV {
    namespace params {
        using poly_t = nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS>;
        template <typename T>
        struct plaintextModulus;
        template <>
        struct plaintextModulus<mpz_class> {
            static mpz_class value() {//19
                /*poly_t  p = (poly_t)1;
                return  mpz_class(p.moduli_product());*///4611686018326724609
                return mpz_class(T);//21267647931552827693776735476788494337
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

void LoadFile(char *filename, int m, int n, int **rawData, int *labels) {
    int x, y=0;
    //bool ret[][] = new bool[m][n];
    // int ret[m][n];

    std::string output;
    std::ifstream in(filename);
    if (in.is_open()) {
        getline(in, output);
        while (!in.eof()) {

            //vector<int> data(m-1);
            for (x = 0; x < m; x++) {
                getline(in, output, ',');
                int tmp = atoi(output.c_str());
                //std::cout<<output<<" ";
                if (x == 0)
                    labels[y]=tmp;
                else {
                    rawData[y][x - 1] = tmp;
                }
                //ret[x][y]=atoi(output.c_str());
                //cout<<ret[x][y]<<endl;
            }

            getline(in, output);
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

mess_t printEncVal(ciphertext_t &ct, sk_t &sk,pk_t &pk){
    mess_t m_dec;
    decrypt(m_dec, sk, pk,ct);
    return m_dec;
}



nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS>  get_array(string line){
    nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS>  res;
    /*line.replace(0,1,"");
    line.replace(line.end()-1,line.end(),"");*/
    char * c = (char*)line.c_str();
    char* pch = strtok (c,",");

    int i=0;
    char * pEnd;
    std::array <mpz_t, DEGREE> arr;
    for (size_t j = 0; j < DEGREE; j++) {
        mpz_init(arr[j]);
    }
    while (pch != NULL)
    {
        //mpz_t &s = arr[i++];
        if(pch!="") {
            mpz_set_str(arr[i++], pch, 10);

            //res[i++] = strtoull(pch,&pEnd,10);
            //cout<<pch<<endl;
            pch = strtok(NULL, ",");
            if(i>DEGREE) break;
        }
    }
    //cout<<i<<endl;

    res.mpz2poly(arr);
    for (size_t j = 0; j < DEGREE; j++) {
        mpz_clear(arr[j]);
    }
    return res;

}
FV::params::poly_p getPolyPredefined(int val){

    FV::params::poly_p pp = (FV::params::poly_p) val;
    std::array <mpz_t, DEGREE> arr  = pp.poly2mpz();
    for (unsigned j = 0; j < arr.size(); j++) {
        mpz_t &s = arr[j];
        mpz_set_ui(s,val);
    }
    pp.mpz2poly(arr);
    return  pp;
}

void readPublicKey( pk_t  &pk){
    string line;
    ifstream myfile ("public.key");
    if (myfile.is_open())
    {
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
    }

    else cout << "Unable to open file";
}
void readSecretKey( sk_t  &sk){
    string line;
    ifstream myfile ("secret.key");
    if (myfile.is_open())
    {
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
    }

    else cout << "Unable to open file";
}

void readEvalKey( evk_t  &evk){
    string line;
    ifstream myevalfile ("eval.key");

    /*int ell = floor(mpz_sizeinbase(P::moduli_product(), 2) / word_size) + 1;
    const int N = ell;
    std::array <P,  N> arr;*/
    cout<< evk.ell<<endl;
    if (myevalfile.is_open())
    {
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
    }

    else cout << "Unable to eval open file";
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

void calculateIteration( ciphertext_t *betaPoly,ciphertext_t** rawDataPoly,ciphertext_t* labelsPoly,int batchsize, int m, unsigned long long int quantization, int iter,pk_t  &pk){
    ciphertext_t xbeta[batchsize];


    for (unsigned y = 0; y < batchsize; y++) {
        //std::cout <<"here "<<std::endl;
        if (iter > 0) {
            ciphertext_t inbetween = rawDataPoly[y ][0] * betaPoly[0];
            for (unsigned x = 1; x < m - 1; x++) {

                //printEncVal(rawDataPoly[y][x],sk,pk);

                /*
                printEncVal(betaPoly[x],sk,pk);*/
                inbetween += (rawDataPoly[y][x] * betaPoly[x]);
                //if(iter>1) {
                /*printEncValPrint(rawDataPoly[y - start][x], sk, pk);
                cout << "*";
                printEncValPrint(betaPoly[x], sk, pk);
                cout << "=";
                ciphertext_t tmp = rawDataPoly[y - start][x] * betaPoly[x];
                printEncValPrint(tmp, sk, pk);
                cout << endl;*/
                //}
            }

            xbeta[y] = inbetween * (24 * quantization );
            //xbeta[y-start] = inbetween*(twenty4*quantization/hundred);
            xbeta[y] = xbeta[y] + 50 * quantization * quantization ;
            xbeta[y] = xbeta[y ] - (inbetween * inbetween);//take care of this later on


        } else {
            FV::encrypt_integer(xbeta[y], pk, 5);//.5 to 5
        }
    }


    for (unsigned y = 0; y <  batchsize; y++) {
        //average+=xbeta[y - start];
        xbeta[y] = labelsPoly[y ] * (quantization * quantization *((iter == 0) ? 10 : 100)) -//((iter == 0) ? hundred : 1))
                           xbeta[y ];

    }
    for (unsigned y = 0; y < m - 1; y++) {
        ciphertext_t inbetween = rawDataPoly[0][y] * xbeta[0];
        for (unsigned x = 1; x <  batchsize; x++) {
            inbetween = inbetween + rawDataPoly[x ][y] * xbeta[x];
        }
        /*printEncValPrint(inbetween, sk, pk);
        cout<<endl;*/
        if (iter > 0) {
            /*std::cout << "adding ";
            printEncValPrint(inbetween, sk, pk);*/

            betaPoly[y] = betaPoly[y] * (quantization * quantization * 10*100) + inbetween;//hundred for alpha
            /*std::cout << " beta ";
            printEncValPrint(betaPoly[y] , sk, pk);
            std::cout << std::endl;*/
        } else
            betaPoly[y] = inbetween;
        //cout<<endl;
    }
}
int main(int argc, char *argv[]) {
    std::cout <<"started"<<std::endl;

    sk_t sk;
    readSecretKey(sk);
    evk_t evk(1 << 6);


    readEvalKey(evk);
    pk_t pk1;
    readPublicKey(pk1);
    //pk.evk=&evk;
    pk_t pk = pk_t(pk1.a,pk1.a_shoup,pk1.b,pk1.b_shoup, evk);
    std::cout <<"Read all the keys"<<std::endl;

    char *datafile;
    int m, n,max_iter=2;
    int batchsize = 0;
    if (argc >= 4) {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        batchsize = (n>atoi(argv[3]))?atoi(argv[3]):n;
        max_iter = atoi(argv[4]);
        //cout << "Starting " << datafile << " covariates " << m << " data " << n << endl;
    } else {
        printf("usage: Filename covariatecount datacount iteration \n");
        return 1;
    }

    ciphertext_t betaPoly[m - 1];
    for (unsigned y = 0; y < m - 1; y++) {
        FV::encrypt_integer(betaPoly[y], pk, 0);
        //printEncValPrint(betaPoly[y],sk,pk);
    }
    //max_iter = batchsize;
    //FV::params::poly_p quantization = getPolyPredefined(1000);
    std::cout<<"log2(q)="<<betaPoly[0].c1.bits_in_moduli_product()<<std::endl;
    std::cout<<"q="<<betaPoly[0].c1.moduli_product()<<std::endl;
    unsigned long long int quantization = 1;//lower starting point
    int fifty = 50;
    int hundred = 100;//alpha
    //int ten = 10;//alpha
    //int thousand=1000;
    int twenty4 = 24;
    int start=0;
    ifstream myfile("labels.encrypted");
    ifstream myfile2("data.encrypted");
    if (myfile.is_open() && myfile2.is_open()) {
        for (int iter = 0; iter < max_iter; iter++) {
            std::cout << "Iteration " << iter << " data " << start << " batchsize " << batchsize <<  std::endl;
            ciphertext_t betaPolyLocal[m - 1];
            for (unsigned y = 0; y < m - 1; y++) {
                FV::encrypt_integer(betaPolyLocal[y], pk, 0);
                //printEncValPrint(betaPoly[y],sk,pk);
            }
            ciphertext_t **rawDataPoly;//[batchsize][m - 1];
            rawDataPoly = new ciphertext_t*[batchsize];
            for (int i = 0; i <batchsize; ++i) {
                rawDataPoly[i] =  new ciphertext_t[m-1];
            }
            ciphertext_t labelsPoly[batchsize];

            for (unsigned y = start; y < start + batchsize; y++) {
                for (unsigned x = 0; x < m - 1; x++) {
                    rawDataPoly[y - start][x].c0.deserialize_manually(myfile2);
                    rawDataPoly[y - start][x].c1.deserialize_manually(myfile2);
                    rawDataPoly[y - start][x].isnull = false;
                    rawDataPoly[y - start][x].pk = &pk;
                }
                labelsPoly[y - start].c0.deserialize_manually(myfile);
                labelsPoly[y - start].c1.deserialize_manually(myfile);
                labelsPoly[y - start].isnull = false;
                labelsPoly[y - start].pk = &pk;
            }
            cout << "Decrypting records" << endl;

            for (unsigned y = start; y < start + batchsize; y++) {
                for (unsigned x = 0; x < m - 1; x++) {
                    printEncValPrint(rawDataPoly[y - start][x], sk, pk);
                }
                printEncValPrint(labelsPoly[y - start], sk, pk);
                cout << endl;
            }

            std::cout << "Finished reading encrypted data" << std::endl;

            //cout<<pk.delta<<endl;

            calculateIteration(betaPolyLocal,rawDataPoly,labelsPoly,batchsize,m, quantization,0,pk);
            calculateIteration(betaPolyLocal,rawDataPoly,labelsPoly,batchsize,m, 10*100,1,pk);

            //free(labelsPoly);
            //for (unsigned y = 0; y < batchsize; y++)  free(rawDataPoly[y]);

            for (unsigned x = 0; x < m - 1; x++) {
                cout << "beta local " << x << " ";
                printEncValPrint(betaPolyLocal[x], sk, pk);
                betaPoly[x]+=betaPolyLocal[x];
                cout << endl;
            }
            /*if (iter == 0)
                quantization = 10 * hundred;//step size =10^-2
            else
                quantization = quantization * quantization * hundred*hundred;*/

            /*mpz_mul_ui(quantization , quantization ,1000);//for step
            mpz_clear(this_iter_quant);*/
            if (start +  batchsize <= n)
                start += batchsize;
            else {
                cout<<start<<" "<<batchsize<<endl;
                start = 0;
                myfile.clear();myfile.seekg(0);
                myfile2.clear();myfile2.seekg(0);
            }/*else
                start = n - (start + batchsize);*/
            //FV::params::plaintextModulus<mpz_class>::value().set_str("21267647931552827693776735476788494337",10);
        }
        myfile.close();
        myfile2.close();
    } else
        cout << "Unable to open labels/data.encrypted file";
    /***
     * decryption block
     */
    ofstream betaEncrypted;
    betaEncrypted.open("beta.encrypted");

    mpz_t tmp;
    for (unsigned i = 0; i < m - 1; i++) {
        betaPoly[i].c0.serialize_manually(betaEncrypted);
        betaPoly[i].c1.serialize_manually(betaEncrypted);
        printEncValPrint(betaPoly[i], sk, pk);
        if (i < m - 2)
            cout << ",";

    }

    betaEncrypted.close();

    std::cout << endl << "quantization " << quantization << endl;

    std::cout <<"ended"<<std::endl;
    return 0;
}
