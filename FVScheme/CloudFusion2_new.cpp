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
#include <thread>
#include <future>
#include <nfl.hpp>
#include "lib/params/params.cpp"

#include "lib/prng/fastrandombytes.cpp"
#include "lib/prng/randombytes.cpp"
//#include "nfl/prng/crypto_stream_salsa20.h"
using namespace std;

int CPU =4;

const int DEGREE = 1 << 10;
const double SIGMA = 3.0;

const int MODULOUS = 62 * 8;
const std::string T = "461168601832";
/// include the FV homomorphic encryption library

typedef int64_t uint48_t;// __attribute__((mode(TI)));
typedef int64_t int48_t;// __attribute__((mode(TI)));
/// include the FV homomorphic encryption library

namespace nfl {
    template<>
    struct params<uint48_t> {
        typedef unsigned int uint128_t __attribute__((mode(TI)));
        typedef uint48_t value_type;
        typedef int48_t signed_value_type;
        typedef uint128_t greater_value_type;

        typedef value_type* poly_t;

        // The moduli used in each 64 bit block (40 bits long each)
        // They are of the form p = 2**40 - i*2*kMaxPolyDegree + 1 for increasing i
        static constexpr unsigned int kMaxNbModuli = 3;
        static constexpr value_type P[kMaxNbModuli] = { 970662608897U, 901943132161U, 850403524609U};

        // The associated lower word of their Newton quotients
        static constexpr value_type Pn[kMaxNbModuli] = {2335245935078U, 3853526466848U, 755220107940U};//{276493118631687U, 21445712413896U, 278631795018165U};

        static constexpr unsigned int kModulusBitsize = 40;
        static constexpr unsigned int kModulusRepresentationBitsize = 42;

        // A primitive 2**15 root of unity for each one of the moduli
        static constexpr value_type primitive_roots[kMaxNbModuli] = { 178079742449U, 899552175130U, 260620111698U};

        // Inverses of kMaxPolyDegree (for the other degrees it can be derived easily)
        // for the different moduli
        static constexpr value_type invkMaxPolyDegree[kMaxNbModuli] = {970603364353U, 901888081921U, 850351620097U};

        // Polynomial related data
        static constexpr unsigned int kMaxPolyDegree = 16384;
    };

    constexpr typename params<uint48_t>::value_type params<uint48_t>::P[];
    constexpr typename params<uint48_t>::value_type params<uint48_t>::Pn[];
    constexpr typename params<uint48_t>::value_type params<uint48_t>::primitive_roots[];
    constexpr typename params<uint48_t>::value_type params<uint48_t>::invkMaxPolyDegree[];

}
namespace FV {
    namespace params {
        using poly_t = nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS>;
        using poly_c = nfl::poly_from_modulus<uint32_t, DEGREE, 2*30>;
        template<typename T>
        struct plaintextModulus;

        template<>
        struct plaintextModulus<mpz_class> {
            static mpz_class value() {
                return mpz_class(poly_c::moduli_product());
            }
        };

        using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
        using gauss_t = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;
        gauss_t fg_prng_sk(SIGMA, 128, 1 << 10);
        gauss_t fg_prng_evk(SIGMA, 128, 1 << 10);
        gauss_t fg_prng_pk(SIGMA, 128, 1 << 10);
        gauss_t fg_prng_enc(SIGMA, 128, 1 << 10);
    }
}

#include "FV_JUAN.hpp"

using namespace FV;

void LoadFile(char *filename, int m, int n, int **rawData, int *labels) {
    int x, y = 0;
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
                    labels[y] = tmp;
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
    printf("read file, rows %d columns %d\n", n, m);
    in.close();
    return;
    //fclose(infile);
    //return ret;
}

mess_t printEncVal(ciphertext_t &ct, sk_t &sk, pk_t &pk) {
    mess_t m_dec;
    decrypt(m_dec, sk, pk, ct);
    return m_dec;
}


nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS> get_array(string line) {
    nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS> res;
    /*line.replace(0,1,"");
    line.replace(line.end()-1,line.end(),"");*/
    char *c = (char *) line.c_str();
    char *pch = strtok(c, ",");

    int i = 0;
    char *pEnd;
    std::array<mpz_t, DEGREE> arr;
    for (size_t j = 0; j < DEGREE; j++) {
        mpz_init(arr[j]);
    }
    while (pch != NULL) {
        //mpz_t &s = arr[i++];
        if (pch != "") {
            mpz_set_str(arr[i++], pch, 10);

            //res[i++] = strtoull(pch,&pEnd,10);
            //cout<<pch<<endl;
            pch = strtok(NULL, ",");
            if (i > DEGREE) break;
        }
    }
    //cout<<i<<endl;

    res.mpz2poly(arr);
    for (size_t j = 0; j < DEGREE; j++) {
        mpz_clear(arr[j]);
    }
    return res;

}

FV::params::poly_p getPolyPredefined(int val) {

    FV::params::poly_p pp = (FV::params::poly_p) val;
    std::array<mpz_t, DEGREE> arr = pp.poly2mpz();
    for (unsigned j = 0; j < arr.size(); j++) {
        mpz_t &s = arr[j];
        mpz_set_ui(s, val);
    }
    pp.mpz2poly(arr);
    return pp;
}

void readPublicKey(pk_t &pk) {
    string line;
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
    string line;
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
    string line;
    ifstream myevalfile("eval.key");

    /*int ell = floor(mpz_sizeinbase(P::moduli_product(), 2) / word_size) + 1;
    const int N = ell;
    std::array <P,  N> arr;*/
    //cout << evk.ell << endl;
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

void printEncValPrint(ciphertext_t &ct, sk_t &sk, pk_t &pk, int times) {
    std::array<mpz_t, DEGREE> arr;
    for (size_t i = 0; i < DEGREE; i++) {
        mpz_inits(arr[i],  nullptr);
    }
    FV::decrypt_poly(arr, sk, pk, ct, true);

    for (int i = 0; i < times; ++i) {

        mpz_t tmp;
        mpz_t &s = arr[i];
        mpz_init(tmp);
        mpz_cdiv_q_ui(tmp, FV::params::plaintextModulus<mpz_class>::value().get_mpz_t(), 2);
        //std::cout<<s<<std::endl;
        if (mpz_cmp(s, tmp) > 0) {
            mpz_sub(tmp, FV::params::plaintextModulus<mpz_class>::value().get_mpz_t(), s);
            //mpz_tdiv_q(tmp,tmp,quantization);
            //std::cout << "beta " << i << ": -" <<tmp << std::endl;
            std::cout << "-" << tmp;
        } else {
            //mpz_tdiv_q(tmp,tmp,quantization);
            //std::cout << "beta " << i << ": " << m_dec << std::endl;
            std::cout << s;
        }
        cout<<" ";
        mpz_clear(tmp);
    }
    for (size_t i = 0; i < DEGREE; i++) {
        mpz_clears(arr[i],  nullptr);
    }
    //return;
}
ciphertext_t calculateInnerProduct( ciphertext_t data1, ciphertext_t data2){
	return data1*data2;
}

ciphertext_t* calculateIteration(int m,ciphertext_t *betaPoly, ciphertext_t **rawDataPoly, ciphertext_t *labelsPoly, int batchsize, unsigned long long int quantization, int iter, pk_t &pk,sk_t &sk) {
    
	ciphertext_t *xbeta;
    xbeta = new ciphertext_t [batchsize];
	
	unsigned y=0,x=0;
	//#pragma omp parallel num_threads(CPU)
	{
		#pragma omp parallel for
		for ( y = 0; y < batchsize; y++) {
			std::cout <<"here "<<std::endl;
			if (iter > 0) {
				ciphertext_t inbetween = rawDataPoly[y][0] * betaPoly[0];
				for ( x = 1; x < m - 1; x++) {
					//printEncVal(rawDataPoly[y][x],sk,pk);
					/*
					printEncVal(betaPoly[x],sk,pk);*/
					inbetween += rawDataPoly[y][x]* betaPoly[x];
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

				xbeta[y] = ciphertext_t(inbetween * (24 * quantization));
				//xbeta[y-start] = inbetween*(twenty4*quantization/hundred);
				xbeta[y] = xbeta[y] + 50 * quantization * quantization;
				xbeta[y] = xbeta[y] - (inbetween * inbetween);//take care of this later on


			}
			/*else {
				
				//FV::encrypt_integer(xbeta[y], pk, 5);//.5 to 5
				xbeta[y] = fives[y];
				
			}*/
			//cout<<"here"<<batchsize<<endl;
			//printEncValPrint(labelsPoly[y], sk, pk,2);
			//#pragma omp critical
			if(iter==0)
				xbeta[y] =  (labelsPoly[y] * 10) -5;
			else
				xbeta[y] =  labelsPoly[y] * (quantization * quantization * ((iter == 0) ? 10 : 100)) //((iter == 0) ? hundred : 1))
						-xbeta[y];
				
		}
		//cout<<endl;
	}
    	cout<<"done 1st iter"<<endl;

	//#pragma omp parallel num_threads(CPU)
	{
		#pragma omp parallel for
		for ( y = 0; y < m - 1; y++) {
			ciphertext_t inbetween = rawDataPoly[0][y] * xbeta[0];
			for ( x = 1; x < batchsize; x++) {
				inbetween = inbetween + rawDataPoly[x][y] * xbeta[x];
			}
			//#pragma omp critical
			if (iter > 0) {
				betaPoly[y] = betaPoly[y] * (quantization * quantization * 10 * 100) + inbetween;//hundred for alpha
			} else
				betaPoly[y] = inbetween;
		}
	}
    return betaPoly;
}


int main(int argc, char *argv[]) {
    std::cout << "started" << std::endl;

    sk_t sk;
    readSecretKey(sk);
    evk_t evk(1 << 6);


    readEvalKey(evk);
    pk_t pk1;
    readPublicKey(pk1);
    //pk.evk=&evk;
    pk_t pk = pk_t(pk1.a, pk1.a_shoup, pk1.b, pk1.b_shoup, evk);
    std::cout << "Read all the keys" << std::endl;

    char *datafile;
    int m, n, max_iter = 2;
    int batchsize = 0;
    if (argc >= 4) {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        batchsize = (n > atoi(argv[3])) ? atoi(argv[3]) : n;
        max_iter = (atoi(argv[4])>2?2:atoi(argv[4]));
        //cout << "Starting " << datafile << " covariates " << m << " data " << n << endl;
    } else {
        printf("usage: Filename covariatecount datacount iteration \n");
        return 1;
    }

    ciphertext_t *betaPoly;
    betaPoly = new ciphertext_t [m-1];
    for (unsigned y = 0; y < m - 1; y++) {
        std::array<mpz_t, DEGREE> ntt_polym0;

        for (size_t j = 0; j < DEGREE; j++) {
            mpz_init_set_ui(ntt_polym0[j], 0);
        }
        FV::params::poly_pl m;
        m.mpz2poly(ntt_polym0);
        FV::encrypt_poly(betaPoly[y], pk, m);

        for (size_t i = 0; i < DEGREE; i++) {
            mpz_clears(ntt_polym0[i], nullptr);
        }
        //printEncValPrint(betaPoly[y],sk,pk);
    }
	/*ciphertext_t *fives;
    fives = new ciphertext_t [batchsize];
    for (unsigned y = 0; y < batchsize; y++) {
        std::array<mpz_t, DEGREE> ntt_polym0;

        for (size_t j = 0; j < DEGREE; j++) {
            mpz_init_set_ui(ntt_polym0[j], 5);
        }
        FV::params::poly_pl m;
        m.mpz2poly(ntt_polym0);
        FV::encrypt_poly(fives[y], pk, m);

        for (size_t i = 0; i < DEGREE; i++) {
            mpz_clears(ntt_polym0[i], nullptr);
        }
        //printEncValPrint(fives[y],sk,pk,2);
    }*/
	
    //max_iter = batchsize;
    //FV::params::poly_p quantization = getPolyPredefined(1000);
    std::cout << "log2(q)=" << betaPoly[0].c1.bits_in_moduli_product() << std::endl;
    std::cout << "q=" << betaPoly[0].c1.moduli_product() << std::endl;
    unsigned long long int quantization = 1;//lower starting point
    int fifty = 50;
    int hundred = 100;//alpha
    //int ten = 10;//alpha
    //int thousand=1000;
    int twenty4 = 24;
    //int start = 0;

    ciphertext_t **rawDataPoly;//[batchsize][m - 1];
    rawDataPoly = new ciphertext_t *[batchsize];
    for (int i = 0; i < batchsize; ++i) {
        rawDataPoly[i] = new ciphertext_t[m - 1];
    }
    ciphertext_t labelsPoly[batchsize];

    ifstream myfile("labels_packed.encrypted");
    ifstream myfile2("data_packed.encrypted");
    if (myfile.is_open() && myfile2.is_open()) {
        for (unsigned y = 0; y < batchsize; y++) {
            for (unsigned x = 0; x < m - 1; x++) {
                rawDataPoly[y][x].c0.deserialize_manually(myfile2);
                rawDataPoly[y][x].c1.deserialize_manually(myfile2);
                rawDataPoly[y][x].isnull = false;
                rawDataPoly[y][x].pk = &pk;
            }
            labelsPoly[y].c0.deserialize_manually(myfile);
            labelsPoly[y].c1.deserialize_manually(myfile);
            labelsPoly[y].isnull = false;
            labelsPoly[y].pk = &pk;
        }
    }else{
        cout<<"Cannot open encrypted file labels_packed.encrypted/data_packed.encrypted \n";
        return 1;
    }
    myfile.close();
    myfile2.close();
    std::cout << "Finished reading encrypted data" << std::endl;
    for (int iter = 0; iter < max_iter; iter++) {
        std::cout << "Iteration " << iter << " batchsize " << batchsize << std::endl;

        //cout<<pk.delta<<endl;

        betaPoly = calculateIteration( m,betaPoly, rawDataPoly, labelsPoly, batchsize, quantization, iter, pk,sk);
        if (iter == 0)
            quantization = 10 * hundred;//step size =10^-2
        else
            quantization = quantization * quantization * hundred*hundred;
		//cout<<"quant "<<quantization<<endl;
        //calculateIteration(betaPoly, rawDataPoly, labelsPoly, batchsize, m, 10 * 100, 1, pk);
        /*for (int i = 0; i <m-1 ; ++i) {
            printEncValPrint(betaPoly[i],sk,pk, n/batchsize);
            cout<<",";
        }
        cout<<endl;*/
    }


    /***
     * decryption block
     */
    ofstream betaEncrypted;
    betaEncrypted.open("beta_packed.encrypted");

    for (unsigned i = 0; i < m - 1; i++) {
        betaPoly[i].c0.serialize_manually(betaEncrypted);
        betaPoly[i].c1.serialize_manually(betaEncrypted);
    }

    betaEncrypted.close();

    //std::cout << endl << "quantization " << quantization << endl;

    std::cout << "ended" << std::endl;
    return 0;
}
