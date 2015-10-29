#ifndef _RIV_H_
#define _RIV_H_

// ---------------------------------------------------------------------

#include <emmintrin.h>
#include "clhash.h"
#include "riv.h"
#include "api.h"

// ---------------------------------------------------------------------

#define BLOCKLEN      CRYPTO_NPUBBYTES
#define KEYLEN        CRYPTO_KEYBYTES
#define TAGLEN        CRYPTO_ABYTES
#define ROUNDS        10
#define ROUND_KEYS    11

typedef __m128i AES_KEY[ROUND_KEYS];
typedef unsigned char block[BLOCKLEN];
typedef int boolean;

// ---------------------------------------------------------------------

typedef struct {
    unsigned char enc_key[BLOCKLEN];
    AES_KEY expanced_enc_key;
    clhash_ctx_t prf_context;
} riv_context_t;

// ---------------------------------------------------------------------

void keysetup(riv_context_t* ctx, const unsigned char key[KEYLEN]);

void encrypt_final(riv_context_t* ctx, 
                   const unsigned char* plaintext,
                   const unsigned long long plaintext_length, 
                   const unsigned char* header,
                   const unsigned long long header_length, 
                   unsigned char* ciphertext, 
                   unsigned char tag[TAGLEN]);

int decrypt_final(riv_context_t* ctx, 
                  const unsigned char* ciphertext,
                  const unsigned long long ciphertext_length, 
                  const unsigned char* header,
                  const unsigned long long header_length, 
                  const unsigned char tag[TAGLEN], 
                  unsigned char* plaintext);

// ---------------------------------------------------------------------

#endif
