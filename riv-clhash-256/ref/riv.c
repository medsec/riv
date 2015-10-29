/**
// RIV[CLHASH-256, Counter[Deoxys-128-128]] reference code.
// Note: This uses CLHASH with four Toeplitz iterations.
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
// Print
// ---------------------------------------------------------------------

#ifdef DEBUG
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
// Deoxys constants
// ---------------------------------------------------------------------

static const block H_PERMUTATION = { 
    7,0,13,10, 11,4,1,14, 15,8,5,2, 3,12,9,6 
};
static const block DOMAIN_MASK = { 
    0x1f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff 
};
static const block DOMAIN_ENC = { 
    0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};
static const unsigned char RCON[17] = {
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 
    0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91, 0x39, 
    0x72
};

// ---------------------------------------------------------------------
// Utils
// ---------------------------------------------------------------------

static inline void to_be_array(uint8_t* result, const uint64_t src)
{
    for (size_t i = 0; i < 8; ++i) {
        result[7-i] = (src >> (8*i)) & 0xFF;
    }
}

// ---------------------------------------------------------------------

static inline void vand(block out, const block a, const block b)
{
    for(size_t i = 0; i < BLOCKLEN; ++i) {
        out[i] = a[i] & b[i];
    }
}

// ---------------------------------------------------------------------

static inline void vor(block out, const block a, const block b)
{
    for(size_t i = 0; i < BLOCKLEN; ++i) {
        out[i] = a[i] | b[i];
    }
}

// ---------------------------------------------------------------------

static inline void vxor(uint8_t* out, const uint8_t* a, const uint8_t* b, 
                        const size_t len)
{
    for(size_t i = 0; i < len; ++i) {
        out[i] = a[i] ^ b[i];
    }
}

// ---------------------------------------------------------------------

static inline void xor_block(block out, const block a, const block b)
{
    vxor(out, a, b, BLOCKLEN);
}

// ---------------------------------------------------------------------
// Permutations
// ---------------------------------------------------------------------

static inline void permute(block x, const block mask) 
{
    block tmp;

    for (size_t i = 0; i < BLOCKLEN; ++i) {
        tmp[i] = x[mask[i]];
    }

    memcpy(x, tmp, BLOCKLEN);
}

// ---------------------------------------------------------------------

static inline void permute_tweak(block x)
{
    permute(x, H_PERMUTATION);
}

// ---------------------------------------------------------------------
// Deoxys key/tweak setup
// ---------------------------------------------------------------------

static inline void gf2_8_double_bytes(block out, const block in)
{
    for (int i = 0; i < BLOCKLEN; ++i) {
        const uint8_t msb = (in[i] & 0x80) >> 7;
        out[i] = (in[i] << 1) ^ (msb * 0x1b);
    }
}

// ---------------------------------------------------------------------

static inline void tweakey_round(block out, const block in) 
{
    gf2_8_double_bytes(out, in);
    permute_tweak(out);
}

// ---------------------------------------------------------------------

static void deoxys_keysetup(uint8_t* subkeys, const uint8_t key[BLOCKLEN])
{
    memcpy(subkeys, key, BLOCKLEN);

    for (size_t i = 0; i < DEOXYS_ROUNDS; ++i) {
        tweakey_round(subkeys + (i+1)*BLOCKLEN, subkeys + i*BLOCKLEN);
    }

    for (size_t i = 0; i <= DEOXYS_ROUNDS; ++i) {
        const block rcon = {
            1,2,4,8,RCON[i],RCON[i],RCON[i],RCON[i],0,0,0,0, 0,0,0,0
        };
        xor_block(subkeys + i*BLOCKLEN, subkeys + i*BLOCKLEN, rcon);
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
    deoxys_keysetup(ctx->expanded_key, ctx->enc_key);

    uint8_t prf_key[CLHASH_KEYLEN];
    const size_t num_prf_key_blocks = CLHASH_KEYLEN / BLOCKLEN;
    
    for (size_t i = 1, j = 0; i <= num_prf_key_blocks; ++i, j += BLOCKLEN) {
        ctr[0] = i;
        aes_encrypt(ctr, &(prf_key[j]), &aes_enc);
    }
    
    clhash_keysetup(&(ctx->prf_context), prf_key);
}

// ---------------------------------------------------------------------
// Encryption
// ---------------------------------------------------------------------

static void deoxys_encrypt(block output, 
                           const block input, 
                           const block tweak, 
                           const DEOXYS_KEY key)
{
    block state;
    xor_block(state, input, key);
    xor_block(state, state, tweak);

    // #ifdef DEBUG
    // print_hex("Round 0 state", state);
    // print_hex("Round 0 tweak", tweak);
    // print_hex("Round 0 key", key);
    // #endif

    block tmp_tweak;
    memcpy(tmp_tweak, tweak, BLOCKLEN);
    permute_tweak(tmp_tweak);

    for(size_t i = 1; i < DEOXYS_ROUNDS; ++i) {
        aes_encrypt_round(state, state, tmp_tweak, key + i*BLOCKLEN);

        // #ifdef DEBUG
        // printf("\nRound %zu\n", i);
        // print_hex("state", state);
        // print_hex("tweak", tmp_tweak);
        // print_hex("key", key + i*BLOCKLEN);
        // #endif

        permute_tweak(tmp_tweak);
    }

    aes_encrypt_round(state, output, tmp_tweak, key + DEOXYS_ROUNDS*BLOCKLEN);

    // #ifdef DEBUG
    // printf("\nRound 14\n");
    // print_hex("state", state);
    // print_hex("tweak", tweak);
    // print_hex("key", key + DEOXYS_ROUNDS*BLOCKLEN);
    // #endif   
}

// ---------------------------------------------------------------------

/**
 * The \Psi_3 3-round tweakable Feistel construction from CDMS'09
 * Coron, Dodis, Mandal, and Seurin: "A Domain Extender for the Ideal Cipher".
 */
static void cdms(uint8_t out[32], const uint8_t in[32], const DEOXYS_KEY key)
{
    const uint8_t* l = in;
    const uint8_t* r = in + BLOCKLEN;
    uint8_t* x = out;
    uint8_t* y = out + BLOCKLEN;

    block s;
    // output, input, tweak, key
    deoxys_encrypt(s, l, r, key); // S <- E_K^{R}(L)
    deoxys_encrypt(y, r, s, key); // Y <- E_K^{S}(R)
    deoxys_encrypt(x, s, y, key); // X <- E_K^{Y}(S)
}

// ---------------------------------------------------------------------

static inline void set_domain_in_tweak(block x, const block mask)
{
    vand(x, x, DOMAIN_MASK);
    vor(x, x, mask);
}

// ---------------------------------------------------------------------

static inline void xor_with_counter(block out, 
                                    const block in, 
                                    const uint64_t counter)
{
    memcpy(out, in, 8);
    to_be_array(out+8, counter);
    vxor(out+8, out+8, in+8, 8);
}

// ---------------------------------------------------------------------

static inline void sct_mode(riv_context_t* ctx, 
                            const uint8_t iv[DEOXYS_IVLEN],
                            const uint8_t* input,
                            const uint64_t length, 
                            uint8_t* output)
{
    uint64_t len = length;
    uint8_t* in = (uint8_t*)input;
    uint8_t* out = (uint8_t*)output;
    uint64_t counter = 0L;

    uint8_t* key = ctx->expanded_key;
    block nonce;
    block tweak;
    memcpy(nonce, iv, BLOCKLEN);
    memcpy(tweak, iv+BLOCKLEN, BLOCKLEN);

    block current_tweak;
    block cipher_output;

    while(len >= BLOCKLEN) {
        xor_with_counter(current_tweak, tweak, counter);
        set_domain_in_tweak(current_tweak, DOMAIN_ENC);
        deoxys_encrypt(cipher_output, nonce, current_tweak, key);
        xor_block(out, cipher_output, in);
        in += BLOCKLEN;
        out += BLOCKLEN;
        len -= BLOCKLEN;
        counter++;
    }
    
    if (len > 0) {
        xor_with_counter(current_tweak, tweak, counter);
        set_domain_in_tweak(current_tweak, DOMAIN_ENC);
        deoxys_encrypt(cipher_output, nonce, current_tweak, key);
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

    #ifdef DEBUG
    print_hex_var("clhash result 1", iv, TAGLEN);
    #endif

    cdms(iv, iv, ctx->expanded_key);

    #ifdef DEBUG
    print_hex_var("iv", iv, TAGLEN);
    #endif

    uint8_t iv_[TAGLEN];
    cdms(iv_, iv, ctx->expanded_key);
    sct_mode(ctx, iv_, plaintext, plaintext_length, ciphertext);
    
    clhash(&(ctx->prf_context), 
        header, header_length, DOMAIN_1, ciphertext, plaintext_length, s);

    #ifdef DEBUG
    print_hex_var("clhash result 2", s, TAGLEN);
    #endif

    cdms(s, s, ctx->expanded_key);

    #ifdef DEBUG
    print_hex_var("s", s, TAGLEN);
    #endif

    vxor(tag, s, iv, TAGLEN);
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
    cdms(iv, iv, ctx->expanded_key);
    vxor(iv, iv, tag, TAGLEN);

    cdms(iv_prime, iv, ctx->expanded_key);
    sct_mode(ctx, iv_prime, ciphertext, ciphertext_length, plaintext);
    
    clhash(&(ctx->prf_context), header, header_length, 
        DOMAIN_0, plaintext, ciphertext_length, iv_prime);
    cdms(iv_prime, iv_prime, ctx->expanded_key);
    return constant_time_memcmp(iv, iv_prime, TAGLEN);
}
