/**
// RIV[CLHASH-256, XOR-CTR[AES-128]] reference code.
// Note: This uses CLHASH with four iterations of the 
//       Toeplitz matrix extension.
// Note: This version might be susceptible to side-channel attacks.
// 
// Author: Eik List. October 2015.
*/
#ifdef DEBUG
    #include<stdio.h>
#endif
#include <emmintrin.h>
#include <smmintrin.h>
#include <wmmintrin.h>
#include <mm_malloc.h>
#include <stdint.h>
#include <string.h>
#include "riv.h"

// ---------------------------------------------------------------------

#define DOMAIN_0        1
#define DOMAIN_1        2

// ---------------------------------------------------------------------
// Load, Store, Helpers
// ---------------------------------------------------------------------

#define loadu(p)         _mm_loadu_si128((__m128i*)p)
#define load(p)          _mm_load_si128((__m128i*)p)
#define storeu(p,x)      _mm_storeu_si128((__m128i*)p, x)
#define store(p,x)       _mm_store_si128((__m128i*)p, x)

#define zero             _mm_setzero_si128()
#define one              _mm_setr_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1) 
#define eight            _mm_setr_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8) 
#define seight           _mm_setr_epi8(0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0) 
#define PERM_MASK        _mm_setr_epi8(0,1,2,3,4,5,6,7,15,14,13,12,11,10,9,8)

#define vaesenc(x,y)     _mm_aesenc_si128(x,y) 
#define vaesdec(x,y)     _mm_aesdec_si128(x,y) 
#define vaesdeclast(x,y) _mm_aesdeclast_si128(x,y) 
#define vinversemc(x)    _mm_aesimc_si128(x)

#define vadd(x,y)        _mm_add_epi64(x,y)
#define vand(x,y)        _mm_and_si128(x,y)
#define vandnot(x,y)     _mm_andnot_si128(x,y)
#define vor(x,y)         _mm_or_si128(x,y)
#define vxor(x,y)        _mm_xor_si128(x,y)
#define vxor3(x,y,z)     _mm_xor_si128(x,_mm_xor_si128(y,z))

// ---------------------------------------------------------------------

static __m128i load_partial(const void *p, size_t n)
{
    if (n == 0) {
        return zero;
    } else if (n % 16 == 0) {
        return _mm_loadu_si128((__m128i*)p);
    } else {
        __m128i tmp;
        size_t i;

        for (i = 0; i < n; ++i) {
            ((char*)&tmp)[i] = ((char*)p)[i];
        }

        return tmp;
    }
}

// ---------------------------------------------------------------------

static void store_partial(const void *p, __m128i x, size_t n)
{
    if (n == 0) {
        return;
    } else if (n >= BLOCKLEN) {
        storeu(p, x);
    } else {
        size_t i;
        uint8_t* p_ = (uint8_t*)p;
        uint8_t* x_ = (uint8_t*)&x;

        for (i = 0; i < n; ++i) {
            p_[i] = x_[i];
        }
    }
}

// ---------------------------------------------------------------------
// Print
// ---------------------------------------------------------------------

#ifdef DEBUG
static void print_128(const char* label, __m128i var)
{
    uint8_t val[BLOCKLEN];
    store((void*)val, var);
    printf("%s\n", label);
    size_t i;

    for (i = 0; i < BLOCKLEN; ++i) {
        printf("%02x ", val[i]);
    }

    puts("\n");
}

// ---------------------------------------------------------------------

static inline void print_hex_var(const char *label, 
                                 const uint8_t *x, 
                                 const int len)
{
    printf("%s: \n", label);
    int i;

    for (i = 0; i < len; i++) {
        if ((i != 0) && (i % 16 == 0)) {
            puts("");
        }

        printf("%02x ", x[i]);
    }

    puts("");
}

// ---------------------------------------------------------------------

static inline void print_hex(const char* label, const uint8_t *c)
{
    print_hex_var(label, c, BLOCKLEN);
}
#endif

// ---------------------------------------------------------------------
// AES
// ---------------------------------------------------------------------

static __m128i aes_keygen_assist(__m128i temp1, __m128i temp2)
{
    __m128i temp3;
    temp2 = _mm_shuffle_epi32(temp2, 0xff);
    temp3 = _mm_slli_si128(temp1, 0x4);
    temp1 = vxor(temp1, temp3);
    temp3 = _mm_slli_si128(temp3, 0x4);
    temp1 = vxor(temp1, temp3);
    temp3 = _mm_slli_si128(temp3, 0x4);
    temp1 = vxor(temp1, temp3);
    temp1 = vxor(temp1, temp2);
    return temp1;
}

// ---------------------------------------------------------------------

#define aes_expand_round_key(round, rcon) \
    tmp = _mm_aeskeygenassist_si128(key[round-1], rcon); \
    key[round] = aes_keygen_assist(key[round-1], tmp)

// ---------------------------------------------------------------------

static void aes_expand_key(__m128i userkey, AES_KEY key)
{
    __m128i tmp;
    key[0] = userkey;
    aes_expand_round_key(1, 0x01);
    aes_expand_round_key(2, 0x02);
    aes_expand_round_key(3, 0x04);
    aes_expand_round_key(4, 0x08);
    aes_expand_round_key(5, 0x10);
    aes_expand_round_key(6, 0x20);
    aes_expand_round_key(7, 0x40);
    aes_expand_round_key(8, 0x80);
    aes_expand_round_key(9, 0x1b);
    aes_expand_round_key(10, 0x36);
}

// ---------------------------------------------------------------------

inline __attribute__((always_inline))
static __m128i aes_encrypt(__m128i in, __m128i* k)
{
    __m128i x = _mm_xor_si128(in, k[0]);
    x = _mm_aesenc_si128(x, k[1]);
    x = _mm_aesenc_si128(x, k[2]);
    x = _mm_aesenc_si128(x, k[3]);
    x = _mm_aesenc_si128(x, k[4]);
    x = _mm_aesenc_si128(x, k[5]);
    x = _mm_aesenc_si128(x, k[6]);
    x = _mm_aesenc_si128(x, k[7]);
    x = _mm_aesenc_si128(x, k[8]);
    x = _mm_aesenc_si128(x, k[9]);
    return _mm_aesenclast_si128(x, k[10]);
}

// ---------------------------------------------------------------------
// Deoxys constants
// ---------------------------------------------------------------------

#define H_PERMUTATION          _mm_setr_epi8(7,0,13,10, 11,4,1,14, 15,8,5,2, 3,12,9,6)
#define MSB_MASK               _mm_set1_epi8(0x80)
#define TRIVIAL_PERMUTATION    _mm_setr_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
#define SIMPLY_1B              _mm_set1_epi8(0x1b)
#define KILL_SHIFT             _mm_set1_epi8(0xfe)
#define DOMAIN_MASK            _mm_setr_epi8(0x1f,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff)
#define DOMAIN_ENC             _mm_setr_epi8(0x20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

static const unsigned char RCON[17] = {
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 
    0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91, 0x39, 
    0x72
};

// ---------------------------------------------------------------------
// Deoxys key/tweak setup
// ---------------------------------------------------------------------

#define gf2_8_double_bytes(out, in) \
    __m128i msbits = vandnot(in, MSB_MASK);  \
    __m128i multi_mask = vor(msbits, TRIVIAL_PERMUTATION);  \
    __m128i rot_cons = _mm_shuffle_epi8(SIMPLY_1B, multi_mask);  \
    __m128i tmp = _mm_slli_epi32(in, 1); \
    tmp = vand(tmp, KILL_SHIFT);  \
    out = vxor(tmp, rot_cons)

// ---------------------------------------------------------------------

#define tweakey_round(out, in) \
    gf2_8_double_bytes(out, in); \
    out = _mm_shuffle_epi8(out, H_PERMUTATION)

// ---------------------------------------------------------------------

void deoxys_keysetup(DEOXYS_KEY subkeys, const __m128i key)
{
    subkeys[0] = key;

    for (size_t i = 0; i < DEOXYS_ROUNDS; ++i) {
        tweakey_round(subkeys[i+1], subkeys[i]);
    }

    for (size_t i = 0; i <= DEOXYS_ROUNDS; ++i) {
        const __m128i rcon = _mm_setr_epi8(
            1,2,4,8,RCON[i],RCON[i],RCON[i],RCON[i],0,0,0,0, 0,0,0,0
        );
        subkeys[i] = vxor(subkeys[i], rcon);
    }
}

// ---------------------------------------------------------------------
// Deoxys en/decryption
// ---------------------------------------------------------------------

#define permute(x)                       _mm_shuffle_epi8(x, PERM_MASK)
#define permute_tweak(x)                 _mm_shuffle_epi8(x, H_PERMUTATION)
#define set_domain_in_tweak(tweak,mask)  vor(vand(DOMAIN_MASK, tweak), mask)

#define add_to_tweak(t,x) \
    t = permute(vadd(permute(t), x))

// y = permute(t, ( 0,1,2,3,4,5,6,7,  15,14,13,12,11,10,9,8 ));
// z = vadd   (y, ( 0,0,0,0, 0,0,0,0, 8,0,0,0, 0,0,0,0 ));
// t = permute(z, ( 0,1,2,3,4,5,6,7,  15,14,13,12,11,10,9,8 ));

// ---------------------------------------------------------------------

static inline void prepare_tweak_counters(__m128i* tweak_ctrs) {
    __m128i tmp = one;

    for (size_t round = 0; round < DEOXYS_ROUND_KEYS; ++round) {
        tweak_ctrs[round*8]   = zero;
        tweak_ctrs[round*8+1] = vadd(tweak_ctrs[round*8+0], tmp);
        tweak_ctrs[round*8+2] = vadd(tweak_ctrs[round*8+1], tmp);
        tweak_ctrs[round*8+3] = vadd(tweak_ctrs[round*8+2], tmp);
        tweak_ctrs[round*8+4] = vadd(tweak_ctrs[round*8+3], tmp);
        tweak_ctrs[round*8+5] = vadd(tweak_ctrs[round*8+4], tmp);
        tweak_ctrs[round*8+6] = vadd(tweak_ctrs[round*8+5], tmp);
        tweak_ctrs[round*8+7] = vadd(tweak_ctrs[round*8+6], tmp);
        tmp = permute_tweak(tmp);
    }
}

// ---------------------------------------------------------------------

static inline void xor_eight(__m128i* states, 
                             const __m128i tmp, 
                             const __m128i* tweak_ctrs) 
{
    states[0] = vxor(tmp, tweak_ctrs[0]);
    states[1] = vxor(tmp, tweak_ctrs[1]);
    states[2] = vxor(tmp, tweak_ctrs[2]);
    states[3] = vxor(tmp, tweak_ctrs[3]);
    states[4] = vxor(tmp, tweak_ctrs[4]);
    states[5] = vxor(tmp, tweak_ctrs[5]);
    states[6] = vxor(tmp, tweak_ctrs[6]);
    states[7] = vxor(tmp, tweak_ctrs[7]);
}

// ---------------------------------------------------------------------

// #define deoxys_enc_round_eight(states, tweak, tweak_ctrs, k) \
//     tmp = vxor(tweak, k); \
//     states[0] = vaesenc(states[0], vxor(*(tweak_ctrs+0), tmp)); \
//     states[1] = vaesenc(states[1], vxor(*(tweak_ctrs+1), tmp)); \
//     states[2] = vaesenc(states[2], vxor(*(tweak_ctrs+2), tmp)); \
//     states[3] = vaesenc(states[3], vxor(*(tweak_ctrs+3), tmp)); \
//     states[4] = vaesenc(states[4], vxor(*(tweak_ctrs+4), tmp)); \
//     states[5] = vaesenc(states[5], vxor(*(tweak_ctrs+5), tmp)); \
//     states[6] = vaesenc(states[6], vxor(*(tweak_ctrs+6), tmp)); \
//     states[7] = vaesenc(states[7], vxor(*(tweak_ctrs+7), tmp)); \
//     tweak = permute_tweak(tweak)

// ---------------------------------------------------------------------

// #define deoxys_enc_eight(states, tweak, tweak_ctrs, k, n) \
//     tmp = vxor(n, tweak); \
//     xor_eight(states, tmp, tweak_ctrs); \
//     tweak = permute_tweak(tweak); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(1*8), k[1]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(2*8), k[2]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(3*8), k[3]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(4*8), k[4]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(5*8), k[5]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(6*8), k[6]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(7*8), k[7]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(8*8), k[8]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(9*8), k[9]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(10*8), k[10]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(11*8), k[11]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(12*8), k[12]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(13*8), k[13]); \
//     deoxys_enc_round_eight(states, tweak, tweak_ctrs+(14*8), k[14])

// ---------------------------------------------------------------------

#define deoxys_enc_round_eight(states, tweak, tweak_ctrs, k) do {\
    tmp = vxor(tweak, k); \
    states[0] = vaesenc(states[0], vxor(*(tweak_ctrs+0), tmp)); \
    states[1] = vaesenc(states[1], vxor(*(tweak_ctrs+1), tmp)); \
    states[2] = vaesenc(states[2], vxor(*(tweak_ctrs+2), tmp)); \
    states[3] = vaesenc(states[3], vxor(*(tweak_ctrs+3), tmp)); \
    states[4] = vaesenc(states[4], vxor(*(tweak_ctrs+4), tmp)); \
    states[5] = vaesenc(states[5], vxor(*(tweak_ctrs+5), tmp)); \
    states[6] = vaesenc(states[6], vxor(*(tweak_ctrs+6), tmp)); \
    states[7] = vaesenc(states[7], vxor(*(tweak_ctrs+7), tmp)); \
} while (0)

// ---------------------------------------------------------------------

#define deoxys_enc_eight(states, tweaks, tweak_ctrs, k, n) do {\
    tmp = vxor(n, tweaks[0]); \
    xor_eight(states, tmp, tweak_ctrs); \
    deoxys_enc_round_eight(states, tweaks[1], tweak_ctrs+(1*8), k[1]); \
    deoxys_enc_round_eight(states, tweaks[2], tweak_ctrs+(2*8), k[2]); \
    deoxys_enc_round_eight(states, tweaks[3], tweak_ctrs+(3*8), k[3]); \
    deoxys_enc_round_eight(states, tweaks[4], tweak_ctrs+(4*8), k[4]); \
    deoxys_enc_round_eight(states, tweaks[5], tweak_ctrs+(5*8), k[5]); \
    deoxys_enc_round_eight(states, tweaks[6], tweak_ctrs+(6*8), k[6]); \
    deoxys_enc_round_eight(states, tweaks[7], tweak_ctrs+(7*8), k[7]); \
    deoxys_enc_round_eight(states, tweaks[8], tweak_ctrs+(8*8), k[8]); \
    deoxys_enc_round_eight(states, tweaks[9], tweak_ctrs+(9*8), k[9]); \
    deoxys_enc_round_eight(states, tweaks[10], tweak_ctrs+(10*8), k[10]); \
    deoxys_enc_round_eight(states, tweaks[11], tweak_ctrs+(11*8), k[11]); \
    deoxys_enc_round_eight(states, tweaks[12], tweak_ctrs+(12*8), k[12]); \
    deoxys_enc_round_eight(states, tweaks[13], tweak_ctrs+(13*8), k[13]); \
    deoxys_enc_round_eight(states, tweaks[14], tweak_ctrs+(14*8), k[14]); \
} while (0)

// ---------------------------------------------------------------------

static inline void load_xor_and_store_eight(__m128i* out, 
                                            const __m128i* in, 
                                            const __m128i* states) 
{
    out[0] = vxor(states[0], loadu(in  ));
    out[1] = vxor(states[1], loadu(in+1));
    out[2] = vxor(states[2], loadu(in+2));
    out[3] = vxor(states[3], loadu(in+3));
    out[4] = vxor(states[4], loadu(in+4));
    out[5] = vxor(states[5], loadu(in+5));
    out[6] = vxor(states[6], loadu(in+6));
    out[7] = vxor(states[7], loadu(in+7));
}

// ---------------------------------------------------------------------

#define print_foo(state, tweak, k, round) \
    printf("Round %d \n", round); \
    print_128("state", state); \
    print_128("tweak", tweak); \
    print_128("key", k)

// ---------------------------------------------------------------------

#define deoxys_enc_round(state, tweak, tweak_ctr, k, round) \
    state = vaesenc(state, vxor3(tweak, tweak_ctr, k)); \
    tweak = permute_tweak(tweak)

// ---------------------------------------------------------------------

#define deoxys_enc(state, tweak, tweak_ctrs, k) \
    state = vxor3(state, tweak, *(tweak_ctrs)); \
    tweak = permute_tweak(tweak); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+1*8), k[1], 1); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+2*8), k[2], 2); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+3*8), k[3], 3); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+4*8), k[4], 4); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+5*8), k[5], 5); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+6*8), k[6], 6); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+7*8), k[7], 7); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+8*8), k[8], 8); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+9*8), k[9], 9); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+10*8), k[10], 10); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+11*8), k[11], 11); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+12*8), k[12], 12); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+13*8), k[13], 13); \
    deoxys_enc_round(state, tweak, *(tweak_ctrs+14*8), k[14], 14)

// ---------------------------------------------------------------------

#define deoxys_encrypt_round_single(state, tweak, k) \
    state = vaesenc(state, vxor(tweak, k)); \
    tweak = permute_tweak(tweak)

// ---------------------------------------------------------------------

#define deoxys_encrypt(state, tweak, k) \
    state = vxor3(state, k[0], tweak); \
    tweak = permute_tweak(tweak); \
    deoxys_encrypt_round_single(state, tweak, k[1]); \
    deoxys_encrypt_round_single(state, tweak, k[2]); \
    deoxys_encrypt_round_single(state, tweak, k[3]); \
    deoxys_encrypt_round_single(state, tweak, k[4]); \
    deoxys_encrypt_round_single(state, tweak, k[5]); \
    deoxys_encrypt_round_single(state, tweak, k[6]); \
    deoxys_encrypt_round_single(state, tweak, k[7]); \
    deoxys_encrypt_round_single(state, tweak, k[8]); \
    deoxys_encrypt_round_single(state, tweak, k[9]); \
    deoxys_encrypt_round_single(state, tweak, k[10]); \
    deoxys_encrypt_round_single(state, tweak, k[11]); \
    deoxys_encrypt_round_single(state, tweak, k[12]); \
    deoxys_encrypt_round_single(state, tweak, k[13]); \
    state = vaesenc(state, vxor(tweak, k[14]))

// ---------------------------------------------------------------------

static inline void load_xor_store_n(__m128i* out, 
                                    const __m128i* in, 
                                    const __m128i* states, 
                                    const size_t num_blocks) 
{
    for (size_t i = 0; i < num_blocks; ++i) {
        out[i] = vxor(states[i], loadu(in+i));
    }
}

// ---------------------------------------------------------------------

static inline void deoxys_enc_n(__m128i* states, 
                                const __m128i tweak, 
                                const __m128i* tweak_ctrs, 
                                const __m128i* k, 
                                const size_t num_blocks, 
                                const __m128i n)
{
    size_t i, j;
    __m128i tmp = vxor(n, tweak);
    __m128i tmp_tweak = tweak;

    for(i = 0; i < num_blocks; i++) {
        states[i] = vxor(tmp, tweak_ctrs[i]);
    }

    tmp_tweak = permute_tweak(tmp_tweak);

    for(j = 1; j < DEOXYS_ROUND_KEYS; j++) {
        tmp = vxor(tmp_tweak, k[j]);

        for(i = 0; i< num_blocks; i++) {
            states[i] = vaesenc(states[i], vxor(tweak_ctrs[j*8+i], tmp));
        }

        tmp_tweak = permute_tweak(tmp_tweak);
    }
}

// ---------------------------------------------------------------------
// RIV
// ---------------------------------------------------------------------

void keysetup(riv_context_t* ctx, const unsigned char key[KEYLEN])
{
    AES_KEY expanded_sk;
    aes_expand_key(loadu(key), expanded_sk);
    
    const __m128i k = aes_encrypt(zero, expanded_sk);
    store(ctx->enc_key, k);

    deoxys_keysetup(ctx->expanded_key, k);

    uint8_t prf_key[CLHASH_KEYLEN];
    const size_t num_blocks = CLHASH_KEYLEN / BLOCKLEN;
    size_t j = 0;

    for (size_t i = 1; i <= num_blocks; ++i) {
        storeu((prf_key + j), aes_encrypt(
            _mm_setr_epi8(i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), expanded_sk));
        j += BLOCKLEN;
    }

    clhash_keysetup(&(ctx->prf_context), prf_key);
}

// ---------------------------------------------------------------------

/**
 * The \Psi_3 3-round tweakable Feistel construction from CDMS'09
 * Coron, Dodis, Mandal, and Seurin: "A Domain Extender for the Ideal Cipher".
 */
static void cdms(uint8_t out[TAGLEN], const uint8_t in[TAGLEN], __m128i* key)
{
    __m128i s = load(in);
    __m128i r = load(in + 1);
    deoxys_encrypt(s, r, key); // S <- E_K^{R}(L)

    __m128i tmp_s = s;
    __m128i y = load(in + 1);
    deoxys_encrypt(y, tmp_s, key); // Y <- E_K^{S}(R)
    store(out + 1, y);

    __m128i x = s;
    deoxys_encrypt(x, y, key); // X <- E_K^{Y}(S)
    store(out, x);
}

// ---------------------------------------------------------------------
// Encryption
// ---------------------------------------------------------------------

static inline uint64_t ceil(const uint64_t x, const uint64_t y)
{
    if (x == 0 || y == 0) {
        return 0;
    }

    return ((x-1) / y) + 1;
}

// ---------------------------------------------------------------------

static inline void sct_mode(riv_context_t* ctx, 
                            const unsigned char iv[DEOXYS_IVLEN],
                            const __m128i* in,
                            const uint64_t length, 
                            __m128i* out)
{
    const __m128i* k = ctx->expanded_key;
    uint64_t len = length;
    
    // ---------------------------------------------------------------------
    // The nonce serves as input to each call of the block cipher.
    // ---------------------------------------------------------------------
    
    const __m128i n = vxor(loadu(iv), k[0]);
    
    // ---------------------------------------------------------------------
    // We use r+1 tweaks to store the tweaks t_0, t_1, ..., t_r for one block
    // for r rounds:
    // tweak_ctr[r][i] = pi^{r}(i)
    // tweak_ctr[r][0] = pi^{r}(0) = 0
    // In each round, we then simply have to have the subtweakey:
    // K[r] xor pi^r(T) xor pi^{r}(i)
    // ---------------------------------------------------------------------
    __m128i tweak_ctrs[DEOXYS_ROUND_KEYS*8];
    prepare_tweak_counters(tweak_ctrs);
    
    // ---------------------------------------------------------------------
    // T, the initial tweak
    // We encode domain the into four least significant bits: 
    // tweak = (0001 || tag).
    // ---------------------------------------------------------------------
    
    const __m128i initial_tweak = set_domain_in_tweak(
        loadu((iv+BLOCKLEN)), DOMAIN_ENC
    );
    __m128i tweak_ctr_base = zero;
    __m128i tweaks[15]; // The permuted tweak for the individual rounds.
    __m128i states[8];
    __m128i tmp;
    uint64_t j = 0;
    
    while (len >= 8*BLOCKLEN) {
        // tweak = vxor(initial_tweak, tweak_ctr_base);
        // deoxys_enc_eight(states, tweak, tweak_ctrs, k, n);
        tweaks[0] = vxor(initial_tweak, tweak_ctr_base);

        for (size_t i = 1; i < 8; ++i) {
            tweaks[i] = permute_tweak(tweaks[i-1]);
        }

        for (size_t i = 8; i <= 14; ++i) {
            tweaks[i] = tweaks[i-8];
        }

        deoxys_enc_eight(states, tweaks, tweak_ctrs, k, n);
        load_xor_and_store_eight(out, in, states);

        len -= 8*BLOCKLEN;
        in += 8;
        out += 8;
        j += 8;

        // Every 256-th block, we have an overflow in the first byte and 
        // have to update the next highest bytes in the counter. 
        if ((j & 0xFF) == 0) { 
            add_to_tweak(tweak_ctr_base, seight);
        } else { // No overflow, increment only the lowest byte in the counter.
            tweak_ctr_base = vadd(tweak_ctr_base, eight);
        }
    }

    tweaks[0] = vxor(initial_tweak, tweak_ctr_base);
    
    const size_t ceil_num_blocks = ceil(len, BLOCKLEN);
    const size_t num_blocks = len / BLOCKLEN;
    const size_t last_block = len % BLOCKLEN;

    deoxys_enc_n(states, tweaks[0], tweak_ctrs, k, ceil_num_blocks, n);
    load_xor_store_n(out, in, states, num_blocks);

    if (last_block != 0) {    
        in += num_blocks;
        out += num_blocks;
        store_partial(out, 
            vxor(states[num_blocks], load_partial(in, last_block)), last_block);
    }
}

// ---------------------------------------------------------------------

static inline void xor_bytes(uint8_t* out, const uint8_t* a, const uint8_t* b, 
                             const size_t len)
{
    for(size_t i = 0; i < len; ++i) {
        out[i] = a[i] ^ b[i];
    }
}

// ---------------------------------------------------------------------

void encrypt_final(riv_context_t* ctx, 
                   const unsigned char* plaintext,
                   const unsigned long long plaintext_length, 
                   const unsigned char* header,
                   const unsigned long long header_length, 
                   unsigned char* ciphertext, 
                   unsigned char tag[TAGLEN])
{
    ALIGN(16) uint8_t iv[TAGLEN];
    ALIGN(16) uint8_t s[TAGLEN];

    clhash(&(ctx->prf_context), header, header_length, DOMAIN_0, 
        plaintext, plaintext_length, iv);

    cdms(iv, iv, ctx->expanded_key);

    ALIGN(16) uint8_t iv_prime[TAGLEN];
    cdms(iv_prime, iv, ctx->expanded_key);

    sct_mode(ctx, iv_prime, (const __m128i*)plaintext, 
        plaintext_length, (__m128i*)ciphertext);
    
    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_1, ciphertext, plaintext_length, s);
    
    cdms(s, s, ctx->expanded_key);
    xor_bytes(tag, s, iv, TAGLEN);
}

// ---------------------------------------------------------------------

int decrypt_final(riv_context_t* ctx, 
                  const unsigned char* ciphertext,
                  const unsigned long long ciphertext_length, 
                  const unsigned char* header,
                  const unsigned long long header_length, 
                  const unsigned char tag[TAGLEN], 
                  unsigned char* plaintext)
{
    ALIGN(16) uint8_t iv[TAGLEN];
    ALIGN(16) uint8_t iv_prime[TAGLEN];

    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_1, ciphertext, ciphertext_length, iv);

    cdms(iv, iv, ctx->expanded_key);
    xor_bytes(iv, iv, tag, TAGLEN);

    cdms(iv_prime, iv, ctx->expanded_key);
    sct_mode(ctx, iv_prime, (const __m128i*)ciphertext, 
        ciphertext_length, (__m128i*)plaintext);
    
    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_0, plaintext, ciphertext_length, iv_prime);

    cdms(iv_prime, iv_prime, ctx->expanded_key);
    return (_mm_testc_si128(load(iv), load(iv_prime)) - 1)
        | (_mm_testc_si128(load((iv+BLOCKLEN)), load((iv_prime+BLOCKLEN))) - 1);
}
