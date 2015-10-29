#pragma once

// ---------------------------------------------------------------------

#include "aes.h"
#include "api.h"
#include "clhash.h"

// ---------------------------------------------------------------------

#define BLOCKLEN             16
#define KEYLEN               CRYPTO_KEYBYTES
#define KEYLEN_BITS          8*KEYLEN
#define TAGLEN               CRYPTO_ABYTES
#define DEOXYS_IVLEN         32
#define DEOXYS_ROUNDS        14
#define DEOXYS_ROUND_KEYS    15
#define DOMAIN_0              1
#define DOMAIN_1              2

// ---------------------------------------------------------------------

typedef unsigned char block[BLOCKLEN];
typedef unsigned char DEOXYS_KEY[DEOXYS_ROUND_KEYS*BLOCKLEN];
typedef int boolean;
typedef struct {
    unsigned char secret_key[KEYLEN];
    unsigned char enc_key[BLOCKLEN];
    DEOXYS_KEY expanded_key;
    clhash_ctx_t prf_context;
} riv_context_t;

// ---------------------------------------------------------------------

void keysetup(riv_context_t* ctx, const unsigned char key[KEYLEN]);

// ---------------------------------------------------------------------

void encrypt_final(riv_context_t* ctx, 
                   const unsigned char* plaintext,
                   const unsigned long long plaintext_length, 
                   const unsigned char* header,
                   const unsigned long long header_length, 
                   unsigned char* ciphertext, 
                   unsigned char tag[TAGLEN]);

// ---------------------------------------------------------------------

int decrypt_final(riv_context_t* ctx, 
                  const unsigned char* ciphertext,
                  const unsigned long long ciphertext_length, 
                  const unsigned char* header,
                  const unsigned long long header_length, 
                  const unsigned char tag[TAGLEN], 
                  unsigned char* plaintext);
