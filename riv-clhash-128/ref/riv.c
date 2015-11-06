/**
// RIV[CLHASH-128, XOR-CTR[AES-128]] reference code.
// Note: This uses CLHASH with two iterations of the 
//       Toeplitz matrix extension.
// Note: This version might be susceptible to side-channel attacks.
// 
// Author: Eik List. October 2015.
*/
#ifdef DEBUG
  #include<stdio.h>
#endif
#include <stdint.h>
#include <string.h>
#include "clhash.h"
#include "riv.h"

// ---------------------------------------------------------------------

#ifdef DEBUG
static void print_hex(const char *label, const uint8_t *c, const size_t len)
{
    printf("%s: \n", label);
    
    for (size_t i = 0; i < len; i++) {
        printf("%02x ", c[i]);
    }

    puts("\n");
}
#endif

// ---------------------------------------------------------------------

static inline void vxor(uint8_t* result, 
                        const uint8_t* left, 
                        const uint8_t* right, 
                        const size_t len)
{
    for (size_t i = 0; i < len; ++i) {
        result[i] = left[i] ^ right[i];
    }
}

// ---------------------------------------------------------------------

static inline void xor_block(block result, const block left, const block right)
{
    vxor(result, left, right, BLOCKLEN);
}

// ---------------------------------------------------------------------

static inline void to_array(uint8_t* result, const uint64_t src)
{
    for (size_t i = 0; i < 8; ++i) {
        result[i] = (src >> (8*i)) & 0xFF;
    }
}

// ---------------------------------------------------------------------

void keysetup(riv_context_t* ctx, const unsigned char key[KEYLEN])
{
    AES_KEY aes_enc;
    memcpy(ctx->secret_key, key, KEYLEN);
    aes_expand_enc_key(ctx->secret_key, KEYLEN_BITS, &aes_enc);

    block ctr;
    memset(ctr, 0, BLOCKLEN);
    
    aes_encrypt(ctr, ctx->enc_key, &aes_enc);
    aes_expand_enc_key(ctx->enc_key, KEYLEN_BITS, &(ctx->expanced_enc_key));
    
    uint8_t prf_key[CLHASH_KEYLEN];
    const size_t num_prf_key_blocks = CLHASH_KEYLEN / BLOCKLEN;
    
    for (size_t i = 1, j = 0; i <= num_prf_key_blocks; ++i, j += BLOCKLEN) {
        ctr[0] = i;
        aes_encrypt(ctr, &(prf_key[j]), &aes_enc);
    }
    
    clhash_keysetup(&(ctx->prf_context), prf_key);
}

// ---------------------------------------------------------------------

static void xor_with_counter(block out, const block iv, const uint64_t counter)
{
    to_array(out, counter);
    size_t i;

    for (i = 0; i < 8; i++) {
        out[i] ^= iv[i];
    }

    for (i = 8; i < BLOCKLEN; i++) {
        out[i] = iv[i];
    }
}

// ---------------------------------------------------------------------

static void counter_mode(riv_context_t* ctx, 
                         const uint8_t iv[BLOCKLEN],
                         const uint8_t* input,
                         const uint64_t length, 
                         const uint8_t* output)
{
    uint64_t len = length;
    uint8_t* in = (uint8_t*)input;
    uint8_t* out = (uint8_t*)output;
    uint64_t counter = 0L;

    block iv_encrypted;
    aes_encrypt(iv, iv_encrypted, &ctx->expanced_enc_key);

    block cipher_input;
    block cipher_output;

    while(len > BLOCKLEN) {
        xor_with_counter(cipher_input, iv_encrypted, counter);
        aes_encrypt(cipher_input, cipher_output, &ctx->expanced_enc_key);
        xor_block(out, cipher_output, in);
        in += BLOCKLEN;
        out += BLOCKLEN;
        len -= BLOCKLEN;
        counter++;
    }
    
    if (len > 0) {
        xor_with_counter(cipher_input, iv_encrypted, counter);
        aes_encrypt(cipher_input, cipher_output, &ctx->expanced_enc_key);
        vxor(out, cipher_output, in, len);
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
    uint8_t iv[TAGLEN];
    uint8_t s[TAGLEN];
    memset(iv, 0x00, TAGLEN);
    memset(s, 0x00, TAGLEN);

    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_0, plaintext, plaintext_length, iv);
    
    aes_encrypt(iv, iv, &ctx->expanced_enc_key);

    counter_mode(ctx, iv, plaintext, plaintext_length, ciphertext);
    
    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_1, ciphertext, plaintext_length, s);
    
    aes_encrypt(s, s, &ctx->expanced_enc_key);
    xor_block(tag, s, iv);
}

// ---------------------------------------------------------------------

static int constant_time_memcmp(const void* av, const void* bv, size_t len) {
    const uint8_t* a = (const uint8_t*) av;
    const uint8_t* b = (const uint8_t*) bv;
    uint8_t result = 0;
    size_t i;

    for (i = 0; i < len; i++) {
        result |= *a ^ *b;
        a++;
        b++;
    }
    
    return (int)result;
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
    uint8_t iv[TAGLEN];
    uint8_t iv_prime[TAGLEN];
    memset(iv, 0x00, TAGLEN);
    memset(iv_prime, 0x00, TAGLEN);

    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_1, ciphertext, ciphertext_length, iv);
    aes_encrypt(iv, iv, &ctx->expanced_enc_key);
    xor_block(iv, iv, tag);
    
    counter_mode(ctx, iv, ciphertext, ciphertext_length, plaintext);
    clhash(&(ctx->prf_context), header, header_length, 
        DOMAIN_0, plaintext, ciphertext_length, iv_prime);
    aes_encrypt(iv_prime, iv_prime, &ctx->expanced_enc_key);
    return constant_time_memcmp(iv, iv_prime, TAGLEN);
}
