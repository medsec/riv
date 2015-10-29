#include <emmintrin.h>
#include <wmmintrin.h>
#include <smmintrin.h>
#include <stdint.h>
#ifdef DEBUG
    #include <stdio.h>
#endif
#include <string.h>
#include "clhash.h"

// ---------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------

#ifndef BLOCKLEN 
#define BLOCKLEN 16 
typedef uint8_t block[BLOCKLEN];
#endif

static const uint8_t MASK_126 = 0x3F;

#define zero          _mm_setzero_si128()
#define R             _mm_cvtsi64_si128(27)
#define TABLE         _mm_setr_epi8(0, 27, 54, 45,108,119, 90, 65,216,195,238,245,180,175,130,153)
#define BSWAP_MASK    _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

// ---------------------------------------------------------------------
// Load and store
// ---------------------------------------------------------------------

static inline __m128i loadu(const void *p) { return _mm_loadu_si128((__m128i*)p); }
static inline __m128i load(const void *p)  { return _mm_load_si128((__m128i*)p);  }
static inline void store(const void *p, __m128i x)  {_mm_store_si128((__m128i*)p, x); }

#define vxor(x,y)     _mm_xor_si128(x, y)
#define vxor3(x,y,z)  _mm_xor_si128(x, _mm_xor_si128(y, z))

// ---------------------------------------------------------------------
// Utils
// ---------------------------------------------------------------------

static inline uint64_t min(const uint64_t a, const uint64_t b) 
{
    return a <= b ? a : b;
}

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

#ifdef DEBUG
static inline void storeu(const void *p, __m128i x)  {_mm_storeu_si128((__m128i*)p, x); }
static void store_partial(uint8_t* p, const uint64_t* x)
{
    uint8_t* x_ = (uint8_t*)x;

    for (size_t i = 0; i < 8; ++i) {
        p[i] = x_[i];
    }
}

// ---------------------------------------------------------------------

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

// ---------------------------------------------------------------------

static void print_hex(const char* message, const __m128i var)
{
    size_t i;
    puts(message);
    uint8_t x[BLOCKLEN];
    storeu(&x, var);

    for (i = 0; i < BLOCKLEN; i++) {
        if ((i != 0) && (i % BLOCKLEN == 0)) {
            puts("");
        }

        printf("%02x ", x[i]);
    }

    puts("\n");
}
#endif

// ---------------------------------------------------------------------

static __m128i load_partial(const void *p, size_t n)
{
    if (n == 0) {
        return zero;
    } else if (n % BLOCKLEN == 0) {
        return loadu(p);
    }

    __m128i tmp = zero;

    for (size_t i = 0; i < n; ++i) {
        ((char*)&tmp)[i] = ((char*)p)[i];
    }

    return tmp;
}

// ---------------------------------------------------------------------

static inline uint64_t reduce(__m128i a)
{
    const __m128i z = _mm_clmulepi64_si128(a, R, 0x01);
    const __m128i y = _mm_shuffle_epi8(TABLE, _mm_srli_si128(z, 8));
    const __m128i x = vxor(a, vxor(y, z));
    return (uint64_t)_mm_cvtsi128_si64(x);
}

// ---------------------------------------------------------------------

static inline __m128i hash_asu(const uint64_t* k, const __m128i poly_sum)
{
    const __m128i x = vxor(load(k), poly_sum);
    return _mm_clmulepi64_si128(x, x, 0x01);
}

// ---------------------------------------------------------------------

static inline __m128i hash_size(const uint64_t k_pprime, const uint64_t len) 
{
    const __m128i x = _mm_set_epi64((__m64)k_pprime, (__m64)len);
    return _mm_clmulepi64_si128(x, x, 0x01);
}

// ---------------------------------------------------------------------

static inline __m128i reduce_mod_2_127(__m128i x01, __m128i x23)
{
    const __m128i shifted_x23 = _mm_srli_si128(x23, 1);
    __m128i s1 = _mm_slli_epi64(x23, 1);
    __m128i s2 = _mm_slli_epi64(x23, 2);
    const __m128i sd1 = _mm_slli_si128(_mm_slli_epi64(shifted_x23, 1), 1);
    const __m128i sd2 = _mm_slli_si128(_mm_slli_epi64(shifted_x23, 2), 1);
    s1 = _mm_or_si128(s1, sd1);
    s2 = _mm_or_si128(s2, sd2);
    return vxor3(x01, s1, s2);
}

// ---------------------------------------------------------------------

static inline __m128i gf128_multiply(__m128i a, __m128i b)
{
    __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5;

    tmp0 = _mm_clmulepi64_si128(a, b, 0x11); // C = [C1:C0] = A1 * B1
    tmp5 = _mm_clmulepi64_si128(a, b, 0x00); // D = [D1:D0] = A0 * B0

    /* Computing E */
    tmp1 = _mm_shuffle_epi32(a, 78); // [A0:A1]
    tmp2 = _mm_shuffle_epi32(b, 78); // [B0:B1]
    tmp1 = _mm_xor_si128(tmp1, a);   // [A1 ^ A0:A1 ^ A0]
    tmp2 = _mm_xor_si128(tmp2, b);   // [B1 ^ B0:B1 ^ B0]
    tmp2 = _mm_clmulepi64_si128(tmp1, tmp2, 0x00); // [E:E1] = (A1 xor A0) * (B1 xor B0)
    
    /* Computing X */
    tmp3 = _mm_xor_si128(tmp0, tmp5);// [          D1 ^ C1:          D0 ^ C0]
    tmp4 = _mm_slli_si128(tmp3, 8);  // [          D0 ^ C0:                0]
    tmp3 = _mm_srli_si128(tmp3, 8);  // [                0:          D1 ^ C1]
    tmp4 = _mm_xor_si128(tmp4, tmp5);// [     D1 ^ D0 ^ C0:               D0]
    tmp3 = _mm_xor_si128(tmp3, tmp0);// [               C1:     D1 ^ C1 ^ C0]
    tmp1 = _mm_slli_si128(tmp2, 8);  // [               E0:                0]
    tmp2 = _mm_srli_si128(tmp2, 8);  // [                0:               E1]

    tmp1 = _mm_xor_si128(tmp4, tmp1);// [E0 ^ D1 ^ D0 ^ C0:               D0] = [X1:X0]
    tmp2 = _mm_xor_si128(tmp3, tmp2);// [               C1:E1 ^ D1 ^ C1 ^ C0] = [X3:X2]

    return reduce_mod_2_127(tmp1, tmp2);
}

// ---------------------------------------------------------------------

static inline __m128i poly_step(__m128i* k, __m128i sum, __m128i term)
{
    return gf128_multiply(*k, vxor(sum, term));
}

// ---------------------------------------------------------------------

// @author Daniel Lemire, Owen Kaser
// @return A1 * B1 xor A2 * B2 xor A3 * B3 xor A4 * B4 mod 2^{127} + 2 + 1.
static inline __m128i gf128_multiply_four(__m128i* key, __m128i* terms) 
{
    __m128i A1 = load(key+3);
    __m128i A2 = load(key+2);
    __m128i A3 = load(key+1);
    __m128i A4 = load(key  );

    __m128i B1 = terms[0];
    __m128i B2 = terms[1];
    __m128i B3 = terms[2];
    __m128i B4 = terms[3];

    __m128i Amix11 = _mm_clmulepi64_si128(A1,B1,0x01);
    __m128i Amix21 = _mm_clmulepi64_si128(A1,B1,0x10);
    __m128i Amix12 = _mm_clmulepi64_si128(A2,B2,0x01);
    __m128i Amix22 = _mm_clmulepi64_si128(A2,B2,0x10);
    __m128i Amix13 = _mm_clmulepi64_si128(A3,B3,0x01);
    __m128i Amix23 = _mm_clmulepi64_si128(A3,B3,0x10);
    __m128i Amix14 = _mm_clmulepi64_si128(A4,B4,0x01);
    __m128i Amix24 = _mm_clmulepi64_si128(A4,B4,0x10);

    __m128i Alow1  = _mm_clmulepi64_si128(A1,B1,0x00);
    __m128i Ahigh1 = _mm_clmulepi64_si128(A1,B1,0x11);
    __m128i Alow2  = _mm_clmulepi64_si128(A2,B2,0x00);
    __m128i Ahigh2 = _mm_clmulepi64_si128(A2,B2,0x11);
    __m128i Alow3  = _mm_clmulepi64_si128(A3,B3,0x00);
    __m128i Ahigh3 = _mm_clmulepi64_si128(A3,B3,0x11);
    __m128i Alow4  = _mm_clmulepi64_si128(A4,B4,0x00);
    __m128i Ahigh4 = _mm_clmulepi64_si128(A4,B4,0x11);

    __m128i Amix1 = _mm_xor_si128(Amix11,Amix21);
    __m128i Amix2 = _mm_xor_si128(Amix12,Amix22);
    __m128i Amix3 = _mm_xor_si128(Amix13,Amix23);
    __m128i Amix4 = _mm_xor_si128(Amix14,Amix24);

    Amix12 = _mm_xor_si128(Amix1,Amix2);
    Amix23 = _mm_xor_si128(Amix3,Amix4);
    __m128i Amix = _mm_xor_si128(Amix12,Amix23);

    Amix1 = _mm_slli_si128(Amix,8);
    Amix2 = _mm_srli_si128(Amix,8);

    __m128i Alow12 = _mm_xor_si128(Alow1,Alow2);
    __m128i Alow34 = _mm_xor_si128(Alow3,Alow4);
    __m128i Alow = _mm_xor_si128(Alow12,Alow34);

    __m128i Ahigh12 = _mm_xor_si128(Ahigh1,Ahigh2);
    __m128i Ahigh34 = _mm_xor_si128(Ahigh3,Ahigh4);
    __m128i Ahigh = _mm_xor_si128(Ahigh12,Ahigh34);

    Alow = _mm_xor_si128(Alow,Amix1);
    Ahigh = _mm_xor_si128(Ahigh,Amix2);

    return reduce_mod_2_127(Alow, Ahigh);
}

// ---------------------------------------------------------------------
// Encoding and chunk-splitting
// ---------------------------------------------------------------------

static inline void encode_domain_and_size(const uint8_t domain, 
                                          const uint64_t hlen, 
                                          const uint64_t mlen, 
                                          uint8_t* result)
{
    to_array(result,   hlen);
    to_array(result+8, mlen);
    result[0] = ((domain << 4) & 0xF0) | (result[0] & 0x0F);
}

// ---------------------------------------------------------------------
// CLNH
// ---------------------------------------------------------------------

static inline void xor_and_multiply_word(const __m128i* m, 
                                         const __m128i* k, 
                                         __m128i* sum) 
{
    const __m128i m_xor_k = vxor(*m, *k);
    const __m128i product = _mm_clmulepi64_si128(m_xor_k, m_xor_k, 0x01);
    *sum = vxor(*sum, product);
}

// ---------------------------------------------------------------------

#define clnh(k, m, num_bytes, sum1, sum2) \
    while (num_bytes >= 32) { \
        xor_and_multiply_word(m, k, sum1); \
        xor_and_multiply_word(m, k+1, sum2); \
        xor_and_multiply_word(m+1, k+1, sum1); \
        xor_and_multiply_word(m+1, k+2, sum2); \
        num_bytes -= 32; \
        m += 2; \
        k += 2; \
    }

// ---------------------------------------------------------------------

static inline void clnh_partial_block_non_zero(__m128i* k, 
                                               __m128i* m, 
                                               size_t len, 
                                               __m128i* sum1, 
                                               __m128i* sum2)
{
    clnh(k, m, len, sum1, sum2);
    const size_t rest = len % 32;
    const size_t less = len - rest;

    k += less / BLOCKLEN;
    m += less / BLOCKLEN;
    len -= less;

    if (len >= BLOCKLEN) {
        xor_and_multiply_word(m, k,   sum1);
        xor_and_multiply_word(m, k+1, sum2);
        
        len -= BLOCKLEN;
        m++;
        k++;
    }

    if (len > 0) {
        const __m128i m1 = load_partial(m, len % BLOCKLEN);
        xor_and_multiply_word(&m1, k,   sum1);
        xor_and_multiply_word(&m1, k+1, sum2);
        
        m++;
        k++;
    }
}

// ---------------------------------------------------------------------

static inline void clnh_partial_block(__m128i* k, 
                                      __m128i* m, 
                                      size_t len, 
                                      __m128i* sum1, 
                                      __m128i* sum2)
{
    *sum1 = zero;
    *sum2 = zero;
    clnh_partial_block_non_zero(k, m, len, sum1, sum2);
}

// ---------------------------------------------------------------------

static inline void clnh_full_block(__m128i* k, 
                                   __m128i* m, 
                                   __m128i* sum1, 
                                   __m128i* sum2)
{
    *sum1 = zero;
    *sum2 = zero;
    size_t len = CL_NH_BLOCKLEN;
    clnh(k, m, len, sum1, sum2);
}

// ---------------------------------------------------------------------

#define clnh_full_block_four(nh_key, m, terms1, terms2) \
    clnh_full_block(nh_key, m, terms1, terms2);     m += CL_NH_NUM_WORDS; \
    clnh_full_block(nh_key, m, terms1+1, terms2+1); m += CL_NH_NUM_WORDS; \
    clnh_full_block(nh_key, m, terms1+2, terms2+2); m += CL_NH_NUM_WORDS; \
    clnh_full_block(nh_key, m, terms1+3, terms2+3); m += CL_NH_NUM_WORDS

// ---------------------------------------------------------------------
// CLHASH main
// ---------------------------------------------------------------------

static inline void short_tail(const uint64_t k_pprime, __m128i nh_sum, 
                              const uint64_t size, uint8_t* target)
{
    const __m128i size_hash = hash_size(k_pprime, size);
    const uint64_t result = reduce(vxor(nh_sum, size_hash));
    memcpy(target, (uint64_t*)(&result), 8);
}

// ---------------------------------------------------------------------

static inline void clhash_short_message(__m128i* nh_key, 
                                        const uint64_t* k_pprime, 
                                        __m128i* h, 
                                        const uint64_t hlen, 
                                        const uint8_t domain, 
                                        __m128i* m, 
                                        const uint64_t mlen, 
                                        const uint64_t num_bytes, 
                                        uint8_t target[CLHASH_HASHLEN])
{
    
    __m128i nh_sum1;
    __m128i nh_sum2;
    clnh_partial_block(nh_key, h, hlen, &nh_sum1, &nh_sum2);

    const uint64_t h_words = ceil(hlen, BLOCKLEN) / BLOCKLEN;
    nh_key += h_words;
    clnh_partial_block_non_zero(nh_key, m, mlen, &nh_sum1, &nh_sum2);
    
    const uint64_t m_words = ceil(mlen, BLOCKLEN) / BLOCKLEN;
    nh_key += m_words;

    ALIGN(BLOCKLEN) 
    uint8_t len_block[BLOCKLEN];
    encode_domain_and_size(domain, hlen, mlen, len_block);
    clnh_partial_block_non_zero(nh_key, (__m128i*)len_block, BLOCKLEN, 
        &nh_sum1, &nh_sum2);

    short_tail(k_pprime[0], nh_sum1, num_bytes, target);
    short_tail(k_pprime[1], nh_sum2, num_bytes, target+8);
}

// ---------------------------------------------------------------------

static inline void long_tail(const uint64_t* k_prime, 
                             const uint64_t k_pprime, 
                             __m128i poly_sum,  
                             const uint64_t num_bytes, 
                             uint8_t* target)
{
    const __m128i size_hash = hash_size(k_pprime, num_bytes);
    const __m128i inner_product_hash = hash_asu(k_prime, poly_sum);
    const uint64_t result = reduce(vxor(size_hash, inner_product_hash));
    memcpy(target, (uint64_t*)(&result), 8);
}

// ---------------------------------------------------------------------

static inline void clhash_long_message(__m128i* orig_nh_key, 
                                       __m128i* poly_key, 
                                       const uint64_t* k_prime, 
                                       const uint64_t* k_pprime, 
                                       __m128i* h, 
                                       const uint64_t hlen, 
                                       const uint8_t domain, 
                                       __m128i* m, 
                                       const uint64_t mlen, 
                                       const uint64_t num_bytes, 
                                       uint8_t* target)
{
    uint64_t h_pos = 0L;
    uint64_t m_pos = 0L;
    
    __m128i poly_sum1 = zero;
    __m128i poly_sum2 = zero;
    
    __m128i poly_term1;
    __m128i poly_term2;
    __m128i* nh_key = orig_nh_key;

    // -----------------------------------------------------------------
    // Compute a_i = CLNH(H_i) for each header block 1 <= i <= h-1.
    // -----------------------------------------------------------------

    while (h_pos + CL_NH_BLOCKLEN <= hlen) {
        clnh_full_block(nh_key, h, &poly_term1, &poly_term2);

        poly_sum1 = poly_step(poly_key, poly_sum1, poly_term1);
        poly_sum2 = poly_step(poly_key+CL_NUM_PRECOMPUTED_POWERS_OF_K, 
            poly_sum2, poly_term2);
        
        h_pos += CL_NH_BLOCKLEN;
        h += CL_NH_QWORDLEN;
    }

    // -----------------------------------------------------------------
    // Junction block: Final header and/or first message block. If the 
    // message is short, this block may also contain domain and lengths.
    // -----------------------------------------------------------------

    const uint64_t h_rest = hlen - h_pos;
    const uint64_t h_junction = ceil(h_rest, BLOCKLEN);
    clnh_partial_block(nh_key, h, h_rest, &poly_term1, &poly_term2);
    nh_key += h_junction / BLOCKLEN;

    uint64_t m_junction = min(mlen, CL_NH_BLOCKLEN - h_junction);

    clnh_partial_block_non_zero(nh_key, m, m_junction, 
        &poly_term1, &poly_term2);
    
    m_pos = m_junction;
    m += m_junction / BLOCKLEN;
    
    m_junction = ceil(m_junction, BLOCKLEN);
    nh_key += m_junction / BLOCKLEN;

    if (h_junction + m_junction == CL_NH_BLOCKLEN) {
        poly_sum1 = poly_step(poly_key, poly_sum1, poly_term1);
        poly_sum2 = poly_step(poly_key+CL_NUM_PRECOMPUTED_POWERS_OF_K, 
            poly_sum2, poly_term2);

        nh_key = orig_nh_key;
        poly_term1 = zero;
        poly_term2 = zero;
    }

    // -----------------------------------------------------------------
    // If the message fills at least 4 more blocks: Compute 4 
    // blocks to process their poly results thereafter in parallel.
    // -----------------------------------------------------------------

    while (m_pos + 4 * CL_NH_BLOCKLEN <= mlen) {
        __m128i poly_terms1[CL_NUM_PRECOMPUTED_POWERS_OF_K];
        __m128i poly_terms2[CL_NUM_PRECOMPUTED_POWERS_OF_K];

        clnh_full_block_four(nh_key, m, poly_terms1, poly_terms2);
        
        poly_terms1[0] = vxor(poly_terms1[0], poly_sum1);
        poly_terms2[0] = vxor(poly_terms2[0], poly_sum2);
        
        poly_sum1 = gf128_multiply_four(poly_key, poly_terms1);
        poly_sum2 = gf128_multiply_four(poly_key+CL_NUM_PRECOMPUTED_POWERS_OF_K, 
            poly_terms2
        );

        m_pos += 4 * CL_NH_BLOCKLEN;
    }

    // -----------------------------------------------------------------
    // If the message fills (at least one) more blocks: 
    // Compute a_i = CLNH(M_i) for each message block 1 <= i <= n-1.
    // -----------------------------------------------------------------

    while (m_pos + CL_NH_BLOCKLEN <= mlen) {
        clnh_full_block(nh_key, m, &poly_term1, &poly_term2);

        poly_sum1 = poly_step(poly_key, poly_sum1, poly_term1);
        poly_sum2 = poly_step(poly_key+CL_NUM_PRECOMPUTED_POWERS_OF_K, 
            poly_sum2, poly_term2);
        
        m_pos += CL_NH_BLOCKLEN;
        m += CL_NH_QWORDLEN;
    }
    
    // -----------------------------------------------------------------
    // Compute a_n = CLNH(M_n) for the final block.
    // -----------------------------------------------------------------

    if (mlen - m_pos != 0) {
        // Must not be called if messeage is fully processed since it would 
        // zeroize the poly terms, which yields wrong results for processing
        // the length later.
        clnh_partial_block(nh_key, m, mlen - m_pos, &poly_term1, &poly_term2);
        const uint64_t m_rest = ceil(mlen - m_pos, BLOCKLEN);

        m += m_rest / BLOCKLEN;
        nh_key += m_rest / BLOCKLEN;

        if (m_rest == CL_NH_BLOCKLEN) {
            poly_sum1 = poly_step(poly_key, poly_sum1, poly_term1);
            poly_sum2 = poly_step(poly_key+CL_NUM_PRECOMPUTED_POWERS_OF_K, 
                poly_sum2, poly_term2);

            nh_key = orig_nh_key;
            poly_term1 = zero;
            poly_term2 = zero;
        }
    }

    // -----------------------------------------------------------------
    // Process length
    // -----------------------------------------------------------------

    ALIGN(BLOCKLEN) 
    uint8_t len_block[BLOCKLEN];
    encode_domain_and_size(domain, hlen, mlen, len_block);

    clnh_partial_block_non_zero(nh_key, (__m128i*)len_block, BLOCKLEN, 
        &poly_term1, &poly_term2);

    // -----------------------------------------------------------------
    // Post-process the final hash value.
    // -----------------------------------------------------------------
    
    long_tail(k_prime,   k_pprime[0], 
        vxor(poly_sum1, poly_term1), num_bytes, target);
    long_tail(k_prime+2, k_pprime[1], 
        vxor(poly_sum2, poly_term2), num_bytes, target+8);
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
    __m128i* h = (__m128i*)header;
    __m128i* m = (__m128i*)message;

    const uint64_t num_bytes = ceil(num_header_bytes, BLOCKLEN)
        + ceil(num_message_bytes, BLOCKLEN) 
        + BLOCKLEN;

    // -----------------------------------------------------------------
    // Short or long message.
    // -----------------------------------------------------------------

    if (num_bytes <= CL_NH_BLOCKLEN) {
        clhash_short_message(
            (__m128i*)(ctx->nh_key), ctx->k_pprime, 
            h, num_header_bytes, domain, 
            m, num_message_bytes, num_bytes, target);
        clhash_short_message(
            (__m128i*)(ctx->nh_key)+2, ctx->k_pprime+2, 
            h, num_header_bytes, domain, 
            m, num_message_bytes, num_bytes, target + BLOCKLEN);
    } else {
        clhash_long_message(
            (__m128i*)(ctx->nh_key), (__m128i*)(ctx->k), ctx->k_prime, ctx->k_pprime, 
            h, num_header_bytes, domain, 
            m, num_message_bytes, num_bytes, target);
        clhash_long_message(
            (__m128i*)(ctx->nh_key)+2, 
            (__m128i*)(ctx->k)+2*CL_NUM_PRECOMPUTED_POWERS_OF_K, 
            ctx->k_prime+4, ctx->k_pprime+2, 
            h, num_header_bytes, domain, 
            m, num_message_bytes, num_bytes, target + BLOCKLEN);
    }
}

// ---------------------------------------------------------------------

static void precompute_poly_keys(uint8_t* key)
{
    const __m128i k1 = load(key);
    __m128i k_power = k1;

    for (size_t i = 1; i < CL_NUM_PRECOMPUTED_POWERS_OF_K; ++i) {
        k_power = gf128_multiply(k1, k_power);
        store(key+(i*CL_NH_QWORDLEN), k_power);
    }
}

// ---------------------------------------------------------------------

void clhash_keysetup(clhash_ctx_t* ctx, const uint8_t key[CLHASH_KEYLEN])
{
    memcpy(ctx->nh_key, key, CL_NH_KEYLEN);


    // Due to Toeplitz, we have several keys K, K', K''. E.g., for 2 Toeplitz 
    // iterations: K_1, K_2, K'_1, K'_2, K''_1, K''_2.
    // 
    // For improving the performance of polynomial hashing, we precompute the 
    // CL_NUM_PRECOMPUTED_POWERS_OF_K powers of K and store it in, 
    // e.g., for (CL_NUM_PRECOMPUTED_POWERS_OF_K = 4) as
    // K_1, K_1^2, K_1^3, K_1^4, K_2, K_2^2, K_2^3, K_2^4. 

    for(size_t i = 0; i < NUM_TOEPLITZ_ITERATIONS; ++i) {
        memcpy(
            ctx->k + i*CL_NUM_BYTES_BETWEEN_K, 
            key + CL_NH_KEYLEN + i*BLOCKLEN, 
            BLOCKLEN
        );

        // Each k is a 126-bit integer, and k[15] the highest byte.
        ctx->k[i*CL_NUM_BYTES_BETWEEN_K + BLOCKLEN-1] &= MASK_126;
        precompute_poly_keys(ctx->k + i*CL_NUM_BYTES_BETWEEN_K);
    }

    memcpy((uint8_t*)(ctx->k_prime), 
        key + CL_NH_KEYLEN + CL_K_LEN, 
        CL_K_PRIME_LEN);
    memcpy((uint8_t*)(ctx->k_pprime), 
        key + CL_NH_KEYLEN + CL_K_LEN + CL_K_PRIME_LEN, 
        CL_K_PPRIME_LEN);
}

