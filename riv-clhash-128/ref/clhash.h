#pragma once

// ---------------------------------------------------------------------

#include <stdint.h>

// ---------------------------------------------------------------------

#if __GNUC__
#define ALIGN(n)      __attribute__ ((aligned(n)))
#elif _MSC_VER
#define ALIGN(n)      __declspec(align(n))
#else
#define ALIGN(n)
#endif

// ---------------------------------------------------------------------

#define CL_NH_BLOCKLEN         1024 // = #bytes per NH block
#define CL_NH_WORDLEN             8 // = #bytes/64 bit
#define CL_NH_QWORDLEN           16 // = #bytes/128 bit
#define NUM_TOEPLITZ_ITERATIONS   2 // 
#define CL_NH_NUM_WORDS          64 // = CL_NH_BLOCKLEN / 16
#define CL_NH_KEYLEN           1040 // = CL_NH_BLOCKLEN + (2 * CL_NH_WORDLEN) * (NUM_TOEPLITZ_ITERATIONS - 1)
#define CL_K_LEN                 32 // = 2 * CL_NH_WORDLEN * NUM_TOEPLITZ_ITERATIONS
#define CL_K_PRIME_LEN           32 // = 2 * CL_NH_WORDLEN * NUM_TOEPLITZ_ITERATIONS
#define CL_K_PPRIME_LEN          16 // =     CL_NH_WORDLEN * NUM_TOEPLITZ_ITERATIONS
#define CLHASH_KEYLEN          1120 // = CL_NH_KEYLEN + CL_K_LEN + CL_K_PRIME_LEN + CL_K_PPRIME_LEN
#define CLHASH_HASHLEN           16 // 8 * NUM_TOEPLITZ_ITERATIONS

ALIGN(16)
typedef struct {
    uint8_t nh_key[CL_NH_KEYLEN];
    uint64_t k[2 * NUM_TOEPLITZ_ITERATIONS]; 
    uint64_t k_prime[2 * NUM_TOEPLITZ_ITERATIONS];
    uint64_t k_pprime[NUM_TOEPLITZ_ITERATIONS];
} clhash_ctx_t;

// ---------------------------------------------------------------------

void clhash_keysetup(clhash_ctx_t* ctx, 
                   const unsigned char key[CLHASH_KEYLEN]);

// ---------------------------------------------------------------------

void clhash(const clhash_ctx_t* ctx, 
            const unsigned char* header, 
            const unsigned long long num_header_bytes, 
            const unsigned char domain, 
            const unsigned char* message, 
            const unsigned long long num_message_bytes, 
            unsigned char target[CLHASH_HASHLEN]);

// ---------------------------------------------------------------------
