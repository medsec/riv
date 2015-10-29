/**
// RIV[HS1, Counter[AES-128]] x86-64-only code.
// Note: This uses HS1 Draft v2 with four iterations of the 
//       Toeplitz matrix extension, like in HS1-SIV (not -LO, not -HI).
// Note: This version might be susceptible to side-channel attacks.
// 
// Author: Eik List. June 2015.
// This uses and modifies the AVX2-optimized HS1v2 code by Romain Dolbeau who 
// optimized the HS1 code by Ted Krovetz.
// 
// This is free and unencumbered software released into the public domain.
//
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
//
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// For more information, please refer to <http://unlicense.org/>
// The authors know of no intellectual property claims relevant to this work.
*/
#ifdef DEBUG
    #include<stdio.h>
#endif
#include <emmintrin.h>
#include <smmintrin.h>
#include <wmmintrin.h>
#include <stdint.h>
#include <string.h>
#include "riv.h"

// ---------------------------------------------------------------------

#define DOMAIN_0        1
#define DOMAIN_1        2

// ---------------------------------------------------------------------
// Load, Store, Helpers
// ---------------------------------------------------------------------

#define loadu(p)       _mm_loadu_si128((__m128i*)p)
#define load(p)        _mm_load_si128((__m128i*)p)
#define storeu(p,x)    _mm_storeu_si128((__m128i*)p, x)
#define store(p,x)     _mm_store_si128((__m128i*)p, x)
#define vxor(x,y)      _mm_xor_si128(x,y)
#define vxor3(x,y,z)   _mm_xor_si128(_mm_xor_si128(x,y),z)
#define vadd(a,b)      _mm_add_epi64(a,b)

#define zero           _mm_setzero_si128()
#define one            _mm_set_epi32(0, 0, 0, 1)
#define two            _mm_set_epi32(0, 0, 0, 2)
#define eight          _mm_set_epi32(0, 0, 0, 8)
#define all_eight      _mm_set_epi32(8, 8, 8, 8)

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

static void print_hex(const char *label, const uint8_t *c, const size_t len)
{
    printf("%s: \n", label);
    size_t i;

    for (i = 0; i < len; i++) {
        printf("%02x ", c[i]);
    }

    puts("\n");
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
// RIV
// ---------------------------------------------------------------------

void keysetup(riv_context_t* ctx, const unsigned char key[KEYLEN])
{
    AES_KEY expanded_sk;
    __m128i sk = loadu(key);
    __m128i k;

    aes_expand_key(sk, expanded_sk);
    k = aes_encrypt(zero, expanded_sk);
    store(ctx->enc_key, k);

    aes_expand_key(k, ctx->expanced_enc_key);

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
// Encryption
// ---------------------------------------------------------------------

#define xor_eight(ciphertext, states, plaintext, k) \
    ciphertext[k+0] = vxor(states[0], plaintext[k+0]); \
    ciphertext[k+1] = vxor(states[1], plaintext[k+1]); \
    ciphertext[k+2] = vxor(states[2], plaintext[k+2]); \
    ciphertext[k+3] = vxor(states[3], plaintext[k+3]); \
    ciphertext[k+4] = vxor(states[4], plaintext[k+4]); \
    ciphertext[k+5] = vxor(states[5], plaintext[k+5]); \
    ciphertext[k+6] = vxor(states[6], plaintext[k+6]); \
    ciphertext[k+7] = vxor(states[7], plaintext[k+7])

// ---------------------------------------------------------------------

#define aes_encrypt_round_eight(states, k) \
    states[0] = _mm_aesenc_si128(states[0], k); \
    states[1] = _mm_aesenc_si128(states[1], k); \
    states[2] = _mm_aesenc_si128(states[2], k); \
    states[3] = _mm_aesenc_si128(states[3], k); \
    states[4] = _mm_aesenc_si128(states[4], k); \
    states[5] = _mm_aesenc_si128(states[5], k); \
    states[6] = _mm_aesenc_si128(states[6], k); \
    states[7] = _mm_aesenc_si128(states[7], k)

// ---------------------------------------------------------------------

#define aes_encrypt_last_round_eight(states, k) \
    states[0] = _mm_aesenclast_si128(states[0], k); \
    states[1] = _mm_aesenclast_si128(states[1], k); \
    states[2] = _mm_aesenclast_si128(states[2], k); \
    states[3] = _mm_aesenclast_si128(states[3], k); \
    states[4] = _mm_aesenclast_si128(states[4], k); \
    states[5] = _mm_aesenclast_si128(states[5], k); \
    states[6] = _mm_aesenclast_si128(states[6], k); \
    states[7] = _mm_aesenclast_si128(states[7], k)

// ---------------------------------------------------------------------

#define aes_eight(states, keys) \
    aes_encrypt_round_eight(states, keys[1]); \
    aes_encrypt_round_eight(states, keys[2]); \
    aes_encrypt_round_eight(states, keys[3]); \
    aes_encrypt_round_eight(states, keys[4]); \
    aes_encrypt_round_eight(states, keys[5]); \
    aes_encrypt_round_eight(states, keys[6]); \
    aes_encrypt_round_eight(states, keys[7]); \
    aes_encrypt_round_eight(states, keys[8]); \
    aes_encrypt_round_eight(states, keys[9]); \
    aes_encrypt_last_round_eight(states, keys[10])
    
// ---------------------------------------------------------------------

static inline void aes_encrypt_n(__m128i *text, int num_blocks,
                                 __m128i *keys)
{
    int i, j;

    for(j = 1; j < 10 ; j++) {
        for(i = 0; i< num_blocks; i++) {
            text[i] = _mm_aesenc_si128(text[i], keys[j]);
        }
    }

    for(i = 0; i < num_blocks; i++) {
        text[i] = _mm_aesenclast_si128(text[i], keys[j]);
    }
}

// ---------------------------------------------------------------------

static inline void counter_mode(riv_context_t* ctx, 
                                __m128i iv, 
                                __m128i* plaintext, 
                                uint64_t len, 
                                __m128i* ciphertext)
{
    __m128i ctr = zero;
    __m128i states[8];
    unsigned int i, k, num_blocks, num_chunks, lastblock, remaining_blocks;

    num_blocks = len / BLOCKLEN;   // len / 16
    lastblock = len % BLOCKLEN; // len mod 16

    if (lastblock != 0) {
        num_blocks++;
    }

    num_chunks = num_blocks >> 3;
    remaining_blocks = num_blocks % 8;

    iv  = vxor(iv, ctx->expanced_enc_key[0]);
    k = 0;

    for(i = 0; i != num_chunks; i++) {
        states[0] = vxor(ctr,iv); ctr = vadd(ctr,one);
        states[1] = vxor(ctr,iv); ctr = vadd(ctr,one);
        states[2] = vxor(ctr,iv); ctr = vadd(ctr,one);
        states[3] = vxor(ctr,iv); ctr = vadd(ctr,one);
        states[4] = vxor(ctr,iv); ctr = vadd(ctr,one);
        states[5] = vxor(ctr,iv); ctr = vadd(ctr,one);
        states[6] = vxor(ctr,iv); ctr = vadd(ctr,one);
        states[7] = vxor(ctr,iv); ctr = vadd(ctr,one);

        aes_eight(states, ctx->expanced_enc_key);

        xor_eight(ciphertext, states, plaintext, k);
        k += 8;
    }

    if (remaining_blocks != 0) {
        k = num_chunks * 8; // position
        ciphertext += k;
        plaintext += k;
        
        for(i = 0; i < remaining_blocks; i++) {
            states[i] = vxor(ctr, iv); ctr = vadd(ctr, one);
        }
        
        aes_encrypt_n(states, remaining_blocks, ctx->expanced_enc_key);
        
        for (i = 0; i < remaining_blocks-1; i++) {
            ciphertext[i] = vxor(states[i], plaintext[i]);
        }
        
        if (lastblock == 0) { // Last block is full
            ciphertext[i] = vxor(states[i], plaintext[i]);
        } else {
            store_partial(ciphertext+i, 
                vxor(
                    load_partial((const void*)(plaintext+i), lastblock), 
                    states[i]
                ), lastblock
            );
        }
    }
}

// ---------------------------------------------------------------------

static inline void encrypt(riv_context_t* ctx, 
                    const __m128i iv,
                    const uint8_t* plaintext,
                    const uint64_t plaintext_length, 
                    uint8_t* ciphertext)
{
    __m128i* m = (__m128i*)plaintext;
    __m128i* c = (__m128i*)ciphertext;
    // __m128i iv_encrypted = aes_encrypt(iv, ctx->expanced_enc_key);
    counter_mode(ctx, iv, m, plaintext_length, c);
}

// ---------------------------------------------------------------------

static inline void decrypt(riv_context_t* ctx, 
                    const __m128i iv, 
                    uint8_t* plaintext,
                    const uint64_t ciphertext_length, 
                    const uint8_t* ciphertext)
{
    __m128i* m = (__m128i*)plaintext;
    __m128i* c = (__m128i*)ciphertext;
    // __m128i iv_encrypted = aes_encrypt(iv, ctx->expanced_enc_key);
    counter_mode(ctx, iv, c, ciphertext_length, m);
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
    ALIGN(16) 
    uint8_t iv[BLOCKLEN];

    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_0, plaintext, plaintext_length, iv);
    const __m128i iv_ = aes_encrypt(load(iv), ctx->expanced_enc_key);
    
    encrypt(ctx, iv_, plaintext, plaintext_length, ciphertext);
    storeu(tag, iv_);
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
    const __m128i iv = loadu(tag);
    decrypt(ctx, iv, plaintext, ciphertext_length, ciphertext);
    
    ALIGN(16) 
    uint8_t iv_prime[BLOCKLEN];
    
    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_0, plaintext, ciphertext_length, iv_prime);
    const __m128i iv_prime_ = aes_encrypt(load(iv_prime), ctx->expanced_enc_key);
    return _mm_testc_si128(iv, iv_prime_) - 1;
}
