#pragma once

// ---------------------------------------------------------------------

#include <emmintrin.h>
#include <wmmintrin.h>
#include <smmintrin.h>
#include <stdint.h>

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
#define CL_NUM_PRECOMPUTED_POWERS_OF_K 4
#define CL_NUM_BYTES_BETWEEN_K   64 // = 2 * CL_NH_WORDLEN * CL_NUM_PRECOMPUTED_POWERS_OF_K

#if __GNUC__
#define ALIGN(n)      __attribute__ ((aligned(n)))
#elif _MSC_VER
#define ALIGN(n)      __declspec(align(n))
#else
#define ALIGN(n)
#endif

// ---------------------------------------------------------------------

ALIGN(16)
typedef struct {
	uint8_t nh_key[CL_NH_KEYLEN];
	// We precompute K, K^2, K^3, K^4
	uint8_t k[CL_NUM_PRECOMPUTED_POWERS_OF_K*CL_NH_QWORDLEN*NUM_TOEPLITZ_ITERATIONS]; 
	uint64_t k_prime[2*NUM_TOEPLITZ_ITERATIONS];
	uint64_t k_pprime[NUM_TOEPLITZ_ITERATIONS];
} clhash_ctx_t;

// ---------------------------------------------------------------------

void clhash(const clhash_ctx_t* ctx, 
            const uint8_t* header, 
            const uint64_t num_header_bytes, 
            const uint8_t domain, 
            const uint8_t* message, 
            const uint64_t num_message_bytes, 
            uint8_t target[CLHASH_HASHLEN]);

void clhash_keysetup(clhash_ctx_t* ctx, const uint8_t key[CLHASH_KEYLEN]);

// ---------------------------------------------------------------------
