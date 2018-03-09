#include <cstddef>

#include <gmpxx.h>
#include <chrono>
#include <iostream>

#include <thread>
#include <vector>
#include <nfl.hpp>
#include "lib/params/params.cpp"

#include "lib/prng/fastrandombytes.cpp"
#include "lib/prng/randombytes.cpp"
#include "nfl/prng/crypto_stream_salsa20.h"
/// include the FV homomorphic encryption library
const int MODULOUS = 248;
const int DEGREE = 1<<4;
namespace FV {
    namespace params {
        using poly_t = nfl::poly_from_modulus<uint64_t, DEGREE, MODULOUS>;
        template <typename T>
        struct plaintextModulus;
        template <>
        struct plaintextModulus<mpz_class> {
            static mpz_class value() {
                return mpz_class("123456789");
            }
        };
        using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
        using gauss_t = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;
        gauss_t fg_prng_sk(8.0, 128, 1 << 14);
        gauss_t fg_prng_evk(8.0, 128, 1 << 14);
        gauss_t fg_prng_pk(8.0, 128, 1 << 14);
        gauss_t fg_prng_enc(8.0, 128, 1 << 14);
    }
}  // namespace FV::params
#include "FV.hpp"

int main() {
    // Seed (for deterministic values)
    srand(0);

    // Keygen
    FV::sk_t secret_key;
    FV::evk_t evaluation_key(secret_key, 32);
    FV::pk_t public_key(secret_key, evaluation_key);

    // Polynomials
    FV::params::poly_p m[3];

    m[0] = {1,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0};//, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    m[1] = {1,0,1,0, 1,0,0,0, 0,0,0,0, 0,0,0,0};
    m[2] = {0,0,1,0, 1,0,0,0, 0,0,0,0, 0,0,0,0};


    FV::params::poly_p summation;
    summation= {0,1,0,0,  0,0,0,0, 0,0,0,0,  0,0,0,0};
    FV::ciphertext_t summationC ;
    FV::encrypt_poly(summationC, public_key,summation);
    /*m[2] = {1,0,0,0, 0,0,0,0};
    m[3] = {1,1,1,1,1,1,1,123456788};//c2
    m[1] = {123456788,123456788,123456788,123456788, 123456788,123456788,123456788,123456788};//c1+
    */

    FV::params::poly_p beta[3];
    for (int j = 0; j < 3; ++j) {
        beta[j] = {1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    }

    FV::params::poly_p IM[3];

    IM[0] = {1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    IM[1] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,123456788,0,0};
    IM[2] = {0,0,0,0, 0,0,0,0, 0,123456788,0,0, 0,0,0,0};
    //m[1].ntt_pow_phi();

    // Encrypt

    std::array<FV::ciphertext_t, 3> c;
    std::array<FV::ciphertext_t, 3> IMc;
    std::array<FV::ciphertext_t, 3> cbeta;
    for (int j = 0; j < 3; ++j) {
        FV::encrypt_poly(c[j], public_key, m[j]);
        FV::encrypt_poly(cbeta[j], public_key, beta[j]);
        FV::encrypt_poly(IMc[j], public_key, IM[j]);
    }
    //mat mult
    FV::ciphertext_t tmp;
    for (int j = 0; j < 3; ++j) {
        tmp+=c[j]*cbeta[j];
    }
    //std::cout<<public_key.<<std::endl;
    for (int j = 0; j < 3; ++j) {
        cbeta[j] = c[j] * tmp;
        if (j > 0)
        cbeta[j] *= summationC*IMc[j];

    }
        //cbeta[j]*=summationC;
    }
    // Initialize polym
    std::array<mpz_t, DEGREE> polym0, polym1,polym2;
    for (size_t i = 0; i < DEGREE; i++) {
        mpz_inits(polym0[i], polym1[i],polym2[i], nullptr);
    }

    // decrypt to the polym
    FV::decrypt_poly(polym0, secret_key, public_key, c[0]);
    FV::decrypt_poly(polym1, secret_key, public_key, c[1]);
    FV::decrypt_poly(polym2, secret_key, public_key, c[2]);

    // Script sage
    /*std::cout << "# Sage script for the verification" << std::endl;
    std::cout << "p=" << FV::params::plaintextModulus<mpz_class>::value().get_str() << std::endl;
    std::cout << "K.<X> = QuotientRing(GF(p)[x], GF(p)[x].ideal(x^DEGREE+1 + 1));"
              << std::endl;
*/
    std::cout << "m0 = ";
    for (size_t i = 0; i < DEGREE; i++) {
        std::cout << mpz_class(polym0[i]).get_str() << " x^" << i
                  << (i == DEGREE-1 ? ";\n" : "+ ");
    }
    std::cout << "m1 = ";
    for (size_t i = 0; i < DEGREE; i++) {
        std::cout << mpz_class(polym1[i]).get_str() << "x^" << i
                  << (i == DEGREE-1 ? ";\n" : "+");
    }

    // Multiplication and decryption
    //FV::ciphertext_t prod = c[0] *c[2] ;
    //prod += c[2] *c[3];
    //prod = prod-c[0]*c[3];
    //prod = prod+prod2;
    FV::decrypt_poly(polym0, secret_key, public_key, cbeta[1]);
    /*FV::params::poly_p tmp;
    tmp.mpz2poly(polym0);
    tmp.ntt_pow_phi();
    polym0 = tmp.poly2mpz();*/

    std::cout << "m = ";
    for (size_t i = 0; i < DEGREE; i++) {

        std::cout << mpz_class(polym0[i]).get_str() << " *X^" << i
                  << (i == DEGREE-1 ? ";\n" : "+ ");
    }

    std::cout << "m = m0*m1 + m2*m1 " << std::endl;

    // Clean
    for (size_t i = 0; i < DEGREE; i++) {
        mpz_clears(polym0[i], polym1[i], nullptr);
    }

    return 0;
}