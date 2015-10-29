#include <stdint.h>
#include <stdlib.h>
#ifdef DEBUG
    #include <stdio.h>
#endif
#include <string.h>
#include "clhash.h"

// ---------------------------------------------------------------------

#ifndef BLOCKLEN 
#define BLOCKLEN 16 
typedef uint8_t block[BLOCKLEN];
#endif

static const uint8_t TABLE[16] = {
    0, 27, 54, 45,108,119, 90, 65,216,195,238,245,180,175,130,153
};

// ---------------------------------------------------------------------

#ifdef DEBUG
static void print_hex_var(const char* message, const uint8_t* x, const size_t len)
{
    puts(message);
    
    for (size_t i = 0; i < len; i++) {
        if ((i != 0) && (i % BLOCKLEN == 0)) {
            puts("");
        }

        printf("%02x ", x[i]);
    }

    puts("\n");
}

static void print_hex(const char* message, const uint8_t x[BLOCKLEN])
{
    print_hex_var(message, x, BLOCKLEN);
}
#endif

// ---------------------------------------------------------------------
// Utils
// ---------------------------------------------------------------------

static inline uint64_t ceil(const uint64_t x, const uint64_t y)
{
    if (x % y == 0) {
        return x;
    }

    return x + (y - (x % y));
}

// ---------------------------------------------------------------------

static inline void to_array(uint8_t* result, const uint64_t src)
{
    for (size_t i = 0; i < 8; ++i) {
        result[7-i] = (src >> (8*i)) & 0xFF;
    }
}

// ---------------------------------------------------------------------

static inline uint64_t to_uint64(const uint8_t* x)
{
    uint64_t result = 0L;

    for (size_t i = 0; i < 8; ++i) {
        result |= (uint64_t)x[i] << (i*8);
    }

    return result;
}

// ---------------------------------------------------------------------

static inline void vxor(block result, const block a, const block b) 
{
    for (size_t i = 0; i < BLOCKLEN; ++i) {
        result[i] = a[i] ^ b[i];
    }
}

// ---------------------------------------------------------------------

static inline void zeroize_block(block in) 
{
    memset(in, 0x00, BLOCKLEN);
}

// ---------------------------------------------------------------------

static inline void zeroize(uint8_t* in, const uint64_t num_bytes) 
{
    memset(in, 0x00, num_bytes);
}

// ---------------------------------------------------------------------

static inline void store(uint8_t* result, const uint64_t value) 
{
    for (size_t i = 0; i < 8; ++i) {
        result[i] = (value >> (i*8)) & 0xFF;
    }
}

// ---------------------------------------------------------------------
// Multiplications
// ---------------------------------------------------------------------

static void gf64_multiply(const uint64_t a, const uint64_t b, uint64_t r[2])
{
    uint8_t w = 64;
    uint8_t s = 4; // Window size
    uint64_t two_s = 1 << s; // 2^s
    uint64_t smask = two_s-1; // s-1 bits
    uint64_t u[two_s];
    uint64_t tmp;
    uint64_t ifmask;

    // -----------------------------------------------------------------
    // Precomputation
    // -----------------------------------------------------------------
    uint64_t i;
    u[0] = 0;
    u[1] = b;

    for(i = 2 ; i < two_s; i += 2) {
        u[i] = u[i >> 1] << 1; // Get even numbers by left shift.
        u[i + 1] = u[i] ^ b;   // Get odd numbers by xoring with b.
    }

    // -----------------------------------------------------------------
    // Multiplication
    // -----------------------------------------------------------------
    
    // The first window affects only the least-significant word.
    r[0] = u[a & smask] ^ (u[(a >> s) & smask] << s); 
    r[1] = 0;

    for(i = 2*s; i < w; i += 2*s) {
        tmp = u[(a >> i) & smask] ^ (u[(a >> (i+s)) & smask] << s);
        r[0] ^= tmp << i;
        r[1] ^= tmp >> (w - i);
    }

    // -----------------------------------------------------------------
    // Reparation
    // -----------------------------------------------------------------
    
    uint64_t m = 0xFEFEFEFEFEFEFEFEL; // s = 2*4
    
    for(i = 1 ; i < 2*s ; i++) {
        tmp = ((a & m) >> i);
        m &= m << 1;
        ifmask = -((b >> (w-i)) & 1); // If the (w-i)-th bit of b is 1.
        r[1] ^= (tmp & ifmask);
    }
}

// ---------------------------------------------------------------------

static inline void reduce_mod_2_127(const uint64_t x01[2], 
                                    const uint64_t x23[2], 
                                    uint64_t r[2])
{
    const uint64_t shifted_x23_lo = ((x23[1] & 0xFF) << 56) | (x23[0] >> 8);
    const uint64_t shifted_x23_hi = x23[1] >> 8;

    uint64_t sd1_lo = shifted_x23_lo << 1;
    uint64_t sd1_hi = shifted_x23_hi << 1;
    sd1_hi = (sd1_hi << 8) | ((sd1_lo >> 56) & 0xFF);
    sd1_lo = (sd1_lo << 8);

    uint64_t sd2_lo = shifted_x23_lo << 2;
    uint64_t sd2_hi = shifted_x23_hi << 2;
    sd2_hi = (sd2_hi << 8) | ((sd2_lo >> 56) & 0xFF);
    sd2_lo = (sd2_lo << 8);

    const uint64_t s1_lo = (x23[0] << 1) | sd1_lo;
    const uint64_t s1_hi = (x23[1] << 1) | sd1_hi;

    const uint64_t s2_lo = (x23[0] << 2) | sd2_lo;
    const uint64_t s2_hi = (x23[1] << 2) | sd2_hi;

    r[0] = s1_lo ^ s2_lo ^ x01[0];
    r[1] = s1_hi ^ s2_hi ^ x01[1];
}

// ---------------------------------------------------------------------

static inline void gf128_multiply(const uint64_t a[2], 
                                  const uint64_t b[2], 
                                  uint64_t r[2])
{
    uint64_t e[2];
    uint64_t x[4];

    gf64_multiply(a[0], b[0], x); // lower = d
    gf64_multiply(a[1], b[1], x+2); // higher = c
    gf64_multiply(a[0] ^ a[1], b[0] ^ b[1], e); // mid = e

    const uint64_t d1Xc0 = x[1] ^ x[2]; // d1 ^ c0
    x[1] = x[0] ^ e[0] ^ d1Xc0; // d0 ^ e0 ^ d1 ^ c0
    x[2] = x[3] ^ e[1] ^ d1Xc0; // d1 ^ e1 ^ d1 ^ c0

    reduce_mod_2_127(x, x+2, r);
}

// ---------------------------------------------------------------------

static inline uint64_t reduce(block a)
{
    uint64_t z[2];
    gf64_multiply(to_uint64(a+8), 27L, z);

    uint8_t index;
    uint64_t y = 0L;

    for (size_t i = 0; i < 8; ++i) {
        index = (z[1] >> (i*8)) & 0x0F;
        y |= (uint64_t)(TABLE[index]) << (i*8);
    }

    return to_uint64(a) ^ y ^ z[0];
}

// ---------------------------------------------------------------------

static inline void hash_asu(const uint64_t* k, 
                            const block poly_sum, 
                            block result)
{
    const uint64_t x0 = k[0] ^ to_uint64(poly_sum);
    const uint64_t x1 = k[1] ^ to_uint64(poly_sum+8);
    gf64_multiply(x0, x1, (uint64_t*)result);
}

// ---------------------------------------------------------------------

static inline void hash_size(const uint64_t k_pprime, 
                             const uint64_t len, 
                             block result) 
{
    gf64_multiply(k_pprime, len, (uint64_t*)result);
}

// ---------------------------------------------------------------------

static inline void poly_step(const uint64_t* k, block sum, block term, block result)
{
    block x;
    vxor(x, sum, term);
    gf128_multiply(k, (const uint64_t*)x, (uint64_t*)result);
}

// ---------------------------------------------------------------------
// CLNH
// ---------------------------------------------------------------------

static inline void xor_and_multiply_word(const block m, 
                                         const block k, 
                                         block sum) 
{
    block m_xor_k;
    vxor(m_xor_k, m, k);

    block product;
    gf64_multiply(to_uint64(m_xor_k), to_uint64(m_xor_k+8), (uint64_t*)product);

    vxor(sum, sum, product);
}

// ---------------------------------------------------------------------

static void clnh_partial_block(const uint8_t* key, 
                               const uint8_t* message, 
                               const size_t len, 
                               block sum1, 
                               block sum2)
{
    zeroize_block(sum1);
    zeroize_block(sum2);

    uint8_t* k1 = (uint8_t*)key;            // Key for Toeplitz iteration #1
    uint8_t* k2 = (uint8_t*)key + BLOCKLEN; // Key for Toeplitz iteration #2

    uint8_t* m1 = (uint8_t*)message;
    uint8_t* m2 = (uint8_t*)message + BLOCKLEN;

    size_t num_bytes = len;

    while (num_bytes >= 2*BLOCKLEN) {
        xor_and_multiply_word(m1, k1, sum1);
        xor_and_multiply_word(m1, k2, sum2);

        xor_and_multiply_word(m2, k1+BLOCKLEN, sum1);
        xor_and_multiply_word(m2, k2+BLOCKLEN, sum2);

        m1 += 2*BLOCKLEN;
        m2 += 2*BLOCKLEN;
        k1 += 2*BLOCKLEN;
        k2 += 2*BLOCKLEN;
        num_bytes -= 2*BLOCKLEN;
    }

    if (num_bytes >= BLOCKLEN) {
        xor_and_multiply_word(m1, k1, sum1);
        xor_and_multiply_word(m1, k2, sum2);

        num_bytes -= BLOCKLEN;
        m1 += BLOCKLEN;
        k1 += BLOCKLEN;
        k2 += BLOCKLEN;
    }

    if (num_bytes > 0) {
        block last_block;
        zeroize_block(last_block);
        memcpy(last_block, m1, num_bytes);
        xor_and_multiply_word(last_block, k1, sum1);
        xor_and_multiply_word(last_block, k2, sum2);
    }
}

// ---------------------------------------------------------------------

static void clnh_full_block(const uint8_t* key, 
                            const uint8_t* message, 
                            block sum1, 
                            block sum2)
{
    zeroize_block(sum1);
    zeroize_block(sum2);

    uint8_t* k1 = (uint8_t*)key;            // Key for Toeplitz iteration #1
    uint8_t* k2 = (uint8_t*)key + BLOCKLEN; // Key for Toeplitz iteration #2

    uint8_t* m1 = (uint8_t*)message;
    uint8_t* m2 = (uint8_t*)message + BLOCKLEN;

    for (size_t i = 0; i < CL_NH_NUM_WORDS; i += 2) {
        xor_and_multiply_word(m1, k1, sum1);
        xor_and_multiply_word(m1, k2, sum2);

        xor_and_multiply_word(m2, k1+BLOCKLEN, sum1);
        xor_and_multiply_word(m2, k2+BLOCKLEN, sum2);

        m1 += 2*BLOCKLEN;
        m2 += 2*BLOCKLEN;
        k1 += 2*BLOCKLEN;
        k2 += 2*BLOCKLEN;
    }
}

// ---------------------------------------------------------------------
// Encoding and chunk-splitting
// ---------------------------------------------------------------------

static void encode_domain_and_size(const uint8_t domain, 
                                   const uint64_t hlen, 
                                   const uint64_t mlen, 
                                   uint8_t* result)
{
    to_array(result,   hlen);
    to_array(result+8, mlen);
    result[0] = ((domain << 4) & 0xF0) | (result[0] & 0x0F);
}

// ---------------------------------------------------------------------

static void copy_header_rest(uint64_t* h_pos, 
                             const uint64_t hlen, 
                             const uint8_t* h, 
                             uint8_t* chunk, 
                             uint64_t* chunklen)
{
    *chunklen = 0L;

    if (*h_pos == hlen) {
        return;
    }

    const uint64_t num_remaining_h_bytes = hlen - *h_pos;
    memcpy(chunk, h + *h_pos, num_remaining_h_bytes);
    *chunklen += num_remaining_h_bytes;

    // -----------------------------------------------------------------
    // Pad header
    // -----------------------------------------------------------------

    const uint64_t num_pad_bytes = BLOCKLEN - (num_remaining_h_bytes % BLOCKLEN); 

    if ((num_remaining_h_bytes % BLOCKLEN) != 0) { 
        zeroize(chunk + *chunklen, num_pad_bytes);
        *chunklen += num_pad_bytes;
    }

    *h_pos = hlen;
}

// ---------------------------------------------------------------------

static void copy_message_rest(uint64_t* m_pos, 
                              const uint64_t mlen, 
                              const uint8_t* m, 
                              const uint64_t hlen, 
                              const uint8_t domain, 
                              uint8_t* chunk, 
                              uint64_t* chunklen, 
                              int* processed_length)
{
    uint64_t num_free_bytes_in_chunk = CL_NH_BLOCKLEN - *chunklen;
    *processed_length = 0;
    
    // -----------------------------------------------------------------
    // Full message block
    // -----------------------------------------------------------------

    if (*m_pos + num_free_bytes_in_chunk <= mlen) {
        memcpy(chunk + *chunklen, m + *m_pos, 
            num_free_bytes_in_chunk);
        *m_pos += num_free_bytes_in_chunk;
        *chunklen += num_free_bytes_in_chunk;
        return;
    }

    // -----------------------------------------------------------------
    // Pad message
    // -----------------------------------------------------------------

    const uint64_t num_remaining_m_bytes = mlen - *m_pos;
    memcpy(chunk + *chunklen, m + *m_pos, num_remaining_m_bytes);
    *chunklen += num_remaining_m_bytes;

    const uint64_t num_pad_bytes = BLOCKLEN - (num_remaining_m_bytes % BLOCKLEN); 

    if ((num_remaining_m_bytes % BLOCKLEN) != 0) { // pad
        zeroize(chunk + *chunklen, num_pad_bytes);
        *chunklen += num_pad_bytes;
    }

    *m_pos = mlen;

    // -----------------------------------------------------------------
    // If was final message block and the current block still has space,
    // encode domain and the lengths.
    // -----------------------------------------------------------------

    if ((*chunklen + BLOCKLEN) <= CL_NH_BLOCKLEN) {
        encode_domain_and_size(domain, hlen, mlen, chunk + *chunklen);
        *chunklen += BLOCKLEN;
        *processed_length = 1;
    }
}

// ---------------------------------------------------------------------
// CLHASH main
// ---------------------------------------------------------------------

static void clhash_short_message(const uint8_t* nh_key, 
                                 const uint64_t* k_pprime, 
                                 const uint8_t* h, 
                                 const uint64_t hlen, 
                                 const uint8_t domain, 
                                 const uint8_t* m, 
                                 const uint64_t mlen, 
                                 const uint64_t num_bytes, 
                                 uint8_t target[CLHASH_HASHLEN])
{
    uint8_t* chunk = (uint8_t*)malloc(CL_NH_BLOCKLEN);
    zeroize(chunk, CL_NH_BLOCKLEN);

    uint64_t h_pos = 0L;
    uint64_t m_pos = 0L;
    uint64_t chunklen = 0L;
    int processed_length = 0;

    copy_header_rest(&h_pos, hlen, h, chunk, &chunklen);
    copy_message_rest(&m_pos, mlen, m, hlen, domain, chunk, &chunklen, 
        &processed_length);

    block nh_sum1;
    block nh_sum2;
    zeroize_block(nh_sum1);
    zeroize_block(nh_sum2);
    clnh_partial_block(nh_key, chunk, chunklen, nh_sum1, nh_sum2);
    
    block size_hash1;
    block size_hash2;

    hash_size(k_pprime[0], num_bytes, size_hash1);
    hash_size(k_pprime[1], num_bytes, size_hash2);

    vxor(nh_sum1, nh_sum1, size_hash1);
    vxor(nh_sum2, nh_sum2, size_hash2);

    store(target,   reduce(nh_sum1));
    store(target+8, reduce(nh_sum2));

    free(chunk);
}

// ---------------------------------------------------------------------

static void clhash_long_message(const uint8_t* nh_key, 
                                const uint64_t* k, 
                                const uint64_t* k_prime, 
                                const uint64_t* k_pprime, 
                                const uint8_t* h, 
                                const uint64_t hlen, 
                                const uint8_t domain, 
                                const uint8_t* m, 
                                const uint64_t mlen, 
                                const uint64_t num_bytes, 
                                uint8_t* target)
{
    uint64_t h_pos = 0L;
    uint64_t m_pos = 0L;
    int processed_length = 0;
    
    block poly_sum1;
    block poly_sum2;
    zeroize_block(poly_sum1);
    zeroize_block(poly_sum2);

    block poly_term1;
    block poly_term2;
    
    // -----------------------------------------------------------------
    // Compute a_i = CLNH(H_i) for each header block 1 <= i <= h-1.
    // -----------------------------------------------------------------

    while (h_pos + CL_NH_BLOCKLEN <= hlen) {
        clnh_full_block(nh_key, h+h_pos, poly_term1, poly_term2);
        poly_step(k,   poly_sum1, poly_term1, poly_sum1);
        poly_step(k+2, poly_sum2, poly_term2, poly_sum2);
        h_pos += CL_NH_BLOCKLEN;
    }

    // -----------------------------------------------------------------
    // Junction block: Final header and/or first message block. If the 
    // message is short, this block may also contain domain and lengths.
    // -----------------------------------------------------------------

    uint8_t* chunk = (uint8_t*)malloc(CL_NH_BLOCKLEN);
    uint64_t chunklen;
    copy_header_rest(&h_pos, hlen, h, chunk, &chunklen);
    copy_message_rest(&m_pos, mlen, m, hlen, domain, chunk, &chunklen, 
        &processed_length);

    if (chunklen == CL_NH_BLOCKLEN) {
        clnh_full_block(nh_key, chunk, poly_term1, poly_term2);
    } else {
        clnh_partial_block(nh_key, chunk, chunklen, 
            poly_term1, poly_term2); 
    }
    
    // -----------------------------------------------------------------
    // If the message fills (at least one) more blocks: 
    // Compute a_i = CLNH(M_i) for each message block 1 <= i <= n-1.
    // -----------------------------------------------------------------

    if (!processed_length) {
        // -----------------------------------------------------------------
        // The poly steps are here since we did not know before whether
        // the message was short; if so, then the final block is not 
        // processed in a poly step. We need not manuelly increment m_pos 
        // after the poly step since this happened inside copy_message_rest.
        // -----------------------------------------------------------------

        poly_step(k,   poly_sum1, poly_term1, poly_sum1);
        poly_step(k+2, poly_sum2, poly_term2, poly_sum2);

        while (m_pos + CL_NH_BLOCKLEN <= mlen) {
            clnh_full_block(nh_key, m+m_pos, poly_term1, poly_term2);
            poly_step(k,   poly_sum1, poly_term1, poly_sum1);
            poly_step(k+2, poly_sum2, poly_term2, poly_sum2);
            m_pos += CL_NH_BLOCKLEN;
        }

        // -----------------------------------------------------------------
        // Compute a_n = CLNH(M_n) for the final block.
        // -----------------------------------------------------------------

        chunklen = 0L;
        copy_message_rest(&m_pos, mlen, m, hlen, domain, chunk, &chunklen, 
            &processed_length);


        if (chunklen == CL_NH_BLOCKLEN) {
            clnh_full_block(nh_key, chunk, poly_term1, poly_term2);

            // -----------------------------------------------------------------
            // If the final message block was full, we still must encode the 
            // domain and header/message lengths.
            // -----------------------------------------------------------------

            if (!processed_length) {
                poly_step(k,   poly_sum1, poly_term1, poly_sum1);
                poly_step(k+2, poly_sum2, poly_term2, poly_sum2);

                chunklen = 0L;
                encode_domain_and_size(domain, hlen, mlen, chunk);
                clnh_partial_block(nh_key, chunk, BLOCKLEN, 
                    poly_term1, poly_term2);
                m_pos = mlen;
            }
        } else {
            clnh_partial_block(nh_key, chunk, chunklen, 
                poly_term1, poly_term2); 
        }
    }
    
    // -----------------------------------------------------------------
    // Post-process the final hash value.
    // -----------------------------------------------------------------
    
    vxor(poly_sum1, poly_sum1, poly_term1);
    vxor(poly_sum2, poly_sum2, poly_term2);

    block size_hash1;
    block size_hash2;
    hash_size(k_pprime[0], num_bytes, size_hash1);
    hash_size(k_pprime[1], num_bytes, size_hash2);
    
    block inner_product_hash1;
    block inner_product_hash2;
    hash_asu(k_prime,   poly_sum1, inner_product_hash1);
    hash_asu(k_prime+2, poly_sum2, inner_product_hash2);

    vxor(size_hash1, size_hash1, inner_product_hash1);
    vxor(size_hash2, size_hash2, inner_product_hash2);
    
    store(target,   reduce(size_hash1));
    store(target+8, reduce(size_hash2));

    free(chunk);
}

// ---------------------------------------------------------------------

void clhash(const clhash_ctx_t* ctx, 
            const uint8_t* header, 
            const uint64_t num_header_bytes, 
            const uint8_t domain, 
            const uint8_t* message, 
            const uint64_t num_message_bytes, 
            uint8_t target[CLHASH_HASHLEN])
{
    const uint8_t* h = (uint8_t*)header;
    const uint8_t* m = (uint8_t*)message;
    const uint64_t num_bytes = ceil(num_header_bytes, BLOCKLEN)
        + ceil(num_message_bytes, BLOCKLEN) 
        + BLOCKLEN;

    // -----------------------------------------------------------------
    // Short or long message.
    // -----------------------------------------------------------------

    if (num_bytes <= CL_NH_BLOCKLEN) {
        clhash_short_message(
            ctx->nh_key, ctx->k_pprime, 
            h, num_header_bytes, domain, 
            m, num_message_bytes, num_bytes, target);
        clhash_short_message(
            ctx->nh_key+2*BLOCKLEN, ctx->k_pprime+2, 
            h, num_header_bytes, domain, 
            m, num_message_bytes, num_bytes, target + BLOCKLEN);
    } else {
        clhash_long_message(
            ctx->nh_key, ctx->k, ctx->k_prime, ctx->k_pprime, 
            h, num_header_bytes, domain, 
            m, num_message_bytes, num_bytes, target);
        clhash_long_message(
            ctx->nh_key+2*BLOCKLEN, ctx->k+4, ctx->k_prime+4, ctx->k_pprime+2, 
            h, num_header_bytes, domain, 
            m, num_message_bytes, num_bytes, target + BLOCKLEN);
    }
}

// ---------------------------------------------------------------------

void clhash_keysetup(clhash_ctx_t* ctx, 
                     const unsigned char key[CLHASH_KEYLEN])
{
    uint8_t* k = (uint8_t*)key;
    memcpy(ctx->nh_key, k, CL_NH_KEYLEN);
    
    // -----------------------------------------------------------------
    // The poly-hash keys
    // -----------------------------------------------------------------
    
    k += CL_NH_KEYLEN;

    // Each k is a 126-bit integer, and k[15] the highest byte.
    const uint64_t MASK_126 = 0x3FFFFFFFFFFFFFFFL;
    const size_t word_size = sizeof(uint64_t);
    size_t i;

    for (i = 0; i < 2 * NUM_TOEPLITZ_ITERATIONS; ++i) {
        ctx->k[i] = to_uint64(k+(i*word_size));

        if (i & 1) { // For odd i, the higher part of the 128-bit keys
            ctx->k[i] &= MASK_126;
        }
    }

    k += CL_K_LEN;

    for (i = 0; i < 2 * NUM_TOEPLITZ_ITERATIONS; ++i) {
        ctx->k_prime[i] = to_uint64(k+(i*word_size));
    }

    k += CL_K_PRIME_LEN;

    for (i = 0; i < NUM_TOEPLITZ_ITERATIONS; ++i) {
        ctx->k_pprime[i] = to_uint64(k+(i*word_size));
    }
}
