//
// Created by tlangminung on 21/10/17.
//
#include <gmpxx.h>
#include <nfl.hpp>


const int DEGREE = 1 << 12;//13-4,12-2
const double SIGMA = 3.0;

const int MODULOUS = 62 * 2;//4 to 2
const int PlainTEXTMODULOUS = 30 * 1;

const std::string T = "461168602";
/// include the FV homomorphic encryption library
typedef int64_t uint48_t;// __attribute__((mode(TI)));
typedef int64_t int48_t;
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
        static constexpr unsigned int kMaxNbModuli = 2;
        static constexpr value_type P[kMaxNbModuli] = {206158430209U, 77309411329U};//{ 970662608897U, 901943132161U, 850403524609U};

        // The associated lower word of their Newton quotients
        static constexpr value_type Pn[kMaxNbModuli] = {1466015503587U, 977343668325U};//{2335245935078U, 3853526466848U, 755220107940U};

        static constexpr unsigned int kModulusBitsize = 39;
        static constexpr unsigned int kModulusRepresentationBitsize = 41;

        // A primitive 2**15 root of unity for each one of the moduli
        static constexpr value_type primitive_roots[kMaxNbModuli] = {146448211898U, 70400686702U};//{ 178079742449U, 899552175130U, 260620111698U};

        // Inverses of kMaxPolyDegree (for the other degrees it can be derived easily)
        // for the different moduli
        static constexpr value_type invkMaxPolyDegree[kMaxNbModuli] = {206145847297U, 77304692737U};//{970603364353U, 901888081921U, 850351620097U};

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
        using poly_c = nfl::poly_from_modulus<uint32_t, DEGREE, PlainTEXTMODULOUS>;
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
        gauss_t fg_prng_sk(SIGMA, 128, 1 << 12);
        gauss_t fg_prng_evk(SIGMA, 128, 1 << 12);
        gauss_t fg_prng_pk(SIGMA, 128, 1 << 12);
        gauss_t fg_prng_enc(SIGMA, 128, 1 << 12);
    }
}
