#include <stdio.h>
#include <string.h>
#include "riv.h"
#include "api.h"

// ---------------------------------------------------------------------

static void print_hex(const char *message, 
                      const unsigned char *x, 
                      const long long unsigned len)
{
    unsigned long long i;
    puts(message);

    for (i = 0; i < len; i++) {
        if ((i != 0) && (i % 16 == 0)) {
            puts("");
        }

        printf("%02x ", x[i]);
    }

    printf("     %llu (octets)\n\n", len);
}

// ---------------------------------------------------------------------

static void print_context(const riv_context_t *ctx)
{
    print_hex("K: ", ctx->enc_key, BLOCKLEN);

    const clhash_ctx_t prf_ctx = ctx->prf_context;
    // print_hex("K_N: ", prf_ctx.nh_key, CL_NH_KEYLEN);
    print_hex("K_P: ", (uint8_t*)prf_ctx.k, CL_K_LEN);
    print_hex("K_A:", (uint8_t*)prf_ctx.k_prime, CL_K_PRIME_LEN);
    print_hex("K_L: ", (uint8_t*)prf_ctx.k_pprime, CL_K_PPRIME_LEN);
}

// ---------------------------------------------------------------------

static void test_output(const riv_context_t *ctx,
                        const unsigned char *k, const unsigned long long klen,
                        const unsigned char *h, const unsigned long long hlen,
                        const unsigned char *m, const unsigned long long mlen,
                        const unsigned char *c, const unsigned long long clen,
                        const unsigned char *t, const unsigned long long tlen)

{
    print_hex("SK: ", k, klen);
    print_context(ctx);
    print_hex("Header/Nonce: ", h, hlen);
    print_hex("Plaintext:", m, mlen);
    print_hex("Ciphertext:", c, clen);
    print_hex("Tag:", t, tlen);
    puts("\n\n");
}

// ---------------------------------------------------------------------

static void free_up(unsigned char* c, 
                    unsigned char* m, 
                    const unsigned long long mlen)
{
    if (mlen != 0) {
        free(m);
        free(c);
    }
}

// ---------------------------------------------------------------------

static int run_test(const unsigned char *k,
                    const unsigned char *h,
                    const unsigned long long hlen,
                    const unsigned char *expected_m,
                    unsigned long long mlen,
                    const unsigned char *expected_c,
                    const unsigned char *expected_t)
{
    riv_context_t ctx;
    unsigned char* c = NULL;
    unsigned char* m = NULL;

    if (mlen != 0) {
        c = (unsigned char*)calloc((size_t)mlen, 1);
        m = (unsigned char*)calloc((size_t)mlen, 1);
    }

    unsigned long long clen = mlen;
    unsigned char t[TAGLEN];

    keysetup(&ctx, k);
    encrypt_final(&ctx, expected_m, mlen, h, hlen, c, t);

    if (((mlen != 0) && memcmp(expected_c, c, clen))
        || memcmp(expected_t, t, TAGLEN)) {
        test_output(&ctx, k, KEYLEN, h, hlen, expected_m, mlen, c, clen, t, TAGLEN);
        puts("Encryption produced incorrect result");
        free_up(c, m, mlen);
        return -1;
    }

    keysetup(&ctx, k);

    const int result = decrypt_final(&ctx, c, clen, h, hlen, t, m);
    test_output(&ctx, k, KEYLEN, h, hlen, m, mlen, c, clen, t, TAGLEN);

    if ((mlen != 0) && memcmp(expected_m, m, mlen)) {
        puts("Decryption produced incorrect result");
        free_up(c, m, mlen);
        return -1;
    }
    
    free_up(c, m, mlen);
    return result;
}

// ---------------------------------------------------------------------

/*static int flexible_size_tests(const uint8_t* kat)
{
    static const size_t MAX_H_LEN = 1024;
    static const size_t MAX_M_LEN = 4096;
    
    static const size_t H_LENS[11] = {0, 1, 2, 8, 16, 32, 64, 128, 256, 512, 1024};
    static const size_t M_LENS[15] = {0, 1, 2, 8, 16, 32, 64, 128, 256, 512, 1023, 1024, 1025, 2048, 4096};

    unsigned char* w = (uint8_t*)malloc(MAX_M_LEN);
    unsigned char* h = (uint8_t*)malloc(MAX_H_LEN);
    unsigned char* m = (unsigned char*)malloc(MAX_M_LEN);
    unsigned char k[KEYLEN];

    size_t i, j;
    int result = 0;

    for(i = 0; i < MAX_H_LEN; ++i) {
        h[i] = 0xFF & (i*193 + 123);
    }

    for(i = 0; i < KEYLEN; ++i) {
        k[i] = 0xFF & (i*191 + 123);
    }

    for(i = 0; i < MAX_M_LEN; ++i) {
        w[i] = 0xFF & (i*181 + 123);
    }

    for(i = 0; i < sizeof H_LENS; ++i) {
        const size_t hlen = H_LENS[i];

        for(j = 0; j < sizeof M_LENS; ++j) {
            const size_t mlen = M_LENS[j];
            const size_t clen = mlen + TAGLEN;
            memcpy(m, w, mlen);

            result |= run_test(k, h, hlen, m, mlen, kat, kat + mlen);
            kat += clen;
        }
    }

    free(w);
    free(h);
    free(m);

    return result;
}*/

// ---------------------------------------------------------------------

static int test1()
{
    unsigned long long mlen = BLOCKLEN;
    const unsigned long long hlen = 0;
    const unsigned char k[KEYLEN] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
    };
    const unsigned char m[BLOCKLEN] = {
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
        0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff
    };
    const unsigned char c[BLOCKLEN] = {
        0x2c, 0x0f, 0xf0, 0x3a, 0x4b, 0x7d, 0x1f, 0x63, 
        0x9c, 0x62, 0xa4, 0x35, 0xb5, 0xac, 0x05, 0xb9
    };
    const unsigned char t[TAGLEN] = {
        0x94, 0x4c, 0xbb, 0xe1, 0x20, 0x41, 0xf7, 0xd5, 
        0xe6, 0x44, 0xc3, 0xc1, 0x81, 0xfb, 0xcc, 0x6f
    };
    return run_test(k, NULL, hlen, m, mlen, c, t);
}

// ---------------------------------------------------------------------

static int test2()
{
    unsigned long long mlen = 56;
    const unsigned long long hlen = BLOCKLEN;
    const unsigned char k[KEYLEN] = {
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
        0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff
    };
    const unsigned char h[BLOCKLEN] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
    };
    const unsigned char m[56] = {
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
        0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
        0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
        0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
        0xde, 0xad, 0xbe, 0xef, 0xde, 0xaf, 0xba, 0xbe
    };
    const unsigned char c[56] = {
        0x26, 0xef, 0xa5, 0xe0, 0x2c, 0x66, 0x08, 0x5c, 
        0x35, 0x26, 0x0a, 0x9a, 0x10, 0x63, 0x2e, 0x49, 
        0xe0, 0xd0, 0x29, 0x42, 0x42, 0xba, 0xdf, 0x3c, 
        0xfd, 0xb0, 0x29, 0x1c, 0x1d, 0xea, 0xce, 0xf4, 
        0xee, 0x6d, 0x19, 0x08, 0x9c, 0x09, 0xd7, 0x55, 
        0xd0, 0x8a, 0x4b, 0xa8, 0xe4, 0x13, 0xac, 0x0f, 
        0xc7, 0x5c, 0x7c, 0x12, 0x94, 0x32, 0x2a, 0x2c
    };
    const unsigned char t[TAGLEN] = {
        0x96, 0x16, 0xa9, 0x5f, 0x63, 0x18, 0x38, 0x1c, 
        0x58, 0x09, 0xd1, 0x05, 0xb1, 0x8d, 0xe1, 0x24
    };
    return run_test(k, h, hlen, m, mlen, c, t);
}

// ---------------------------------------------------------------------

static int test3()
{
    unsigned long long mlen = 0;
    const unsigned long long hlen = 24;
    const unsigned char k[KEYLEN] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
    };
    const unsigned char h[24] = {
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
        0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff,
        0xde, 0xad, 0xbe, 0xef, 0xde, 0xaf, 0xba, 0xbe
    };
    const unsigned char t[TAGLEN] = {
        0x2b, 0xd3, 0x5f, 0x61, 0x02, 0x8c, 0xae, 0xbc, 
        0x55, 0x09, 0xcd, 0x69, 0xcd, 0xcc, 0x92, 0x1f
    };
    return run_test(k, h, hlen, NULL, mlen, NULL, t);
}

// ---------------------------------------------------------------------

static int test4()
{
    unsigned long long mlen = 52;
    const unsigned long long hlen = 24;
    const unsigned char k[KEYLEN] = {
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
        0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff
    };
    const unsigned char h[24] = {
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
        0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff,
        0xde, 0xad, 0xbe, 0xef, 0xde, 0xaf, 0xba, 0xbe
    };
    const unsigned char m[52] = {
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
        0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
        0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
        0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
        0xfe, 0xfe, 0xba, 0xbe
    };
    const unsigned char c[52] = {
        0x04, 0xb2, 0xee, 0xb8, 0x18, 0x7c, 0x52, 0xac, 
        0xda, 0xad, 0x68, 0x00, 0x9c, 0x1e, 0x7c, 0xbf, 
        0x5e, 0xeb, 0x20, 0x67, 0x95, 0x99, 0x0a, 0x0c, 
        0x00, 0xd8, 0xe6, 0xc0, 0x24, 0xda, 0xce, 0x70, 
        0x7c, 0x1f, 0x31, 0x0b, 0xb7, 0x2f, 0x84, 0x2b, 
        0xc2, 0xf7, 0x12, 0x57, 0xd4, 0xd9, 0xfb, 0x02, 
        0xd3, 0x59, 0xfe, 0xbf
    };
    const unsigned char t[TAGLEN] = {
        0xdb, 0x5d, 0x66, 0x2d, 0x24, 0xb7, 0xa7, 0xea, 
        0x1a, 0xf7, 0x65, 0x67, 0x04, 0x0f, 0x20, 0x1c
    };
    return run_test(k, h, hlen, m, mlen, c, t);
}

// ---------------------------------------------------------------------

static int test5()
{
    const char k[] = "Edgar Allan Poe.";
    const char h[] = "\"Seldom we find,\" says Solomon Don Dunce,\n\"Half an idea in the profoundest sonnet.\nThrough all the flimsy things we see at once\nAs easily as through a Naples bonnet-\nTrash of all trash!- how can a lady don it?\nYet heavier far than your Petrarchan stuff-\nOwl-downy nonsense that the faintest puff\nTwirls into trunk-paper the while you con it.\"\nAnd, veritably, Sol is right enough.\nThe general tuckermanities are arrant\nBubbles- ephemeral and so transparent-\nBut this is, now- you may depend upon it-\nStable, opaque, immortal- all by dint\nOf the dear names that he concealed within 't.";
    char m[] = "The noblest name in Allegory's page,\nThe hand that traced inexorable rage;\nA pleasing moralist whose page refined,\nDisplays the deepest knowledge of the mind;\nA tender poet of a foreign tongue,\n(Indited in the language that he sung.)\nA bard of brilliant but unlicensed page\nAt once the shame and glory of our age,\nThe prince of harmony and stirling sense,\nThe ancient dramatist of eminence,\nThe bard that paints imagination's powers,\nAnd him whose song revives departed hours,\nOnce more an ancient tragic bard recall,\nIn boldness of design surpassing all.\nThese names when rightly read, a name [make] known\nWhich gathers all their glories in its own.";

    const unsigned long long hlen = (unsigned long long)strlen(h);
    unsigned long long mlen = (unsigned long long)strlen(m);
    const unsigned char c[650] = {
        0x9e, 0xe2, 0x45, 0xa7, 0x90, 0xfc, 0xbf, 0xfb, 
        0x9c, 0xe0, 0x77, 0xa3, 0x47, 0x2a, 0x21, 0xd4, 
        0x03, 0xe3, 0x8b, 0x88, 0x98, 0x6a, 0xe7, 0x1a, 
        0x92, 0x36, 0xf3, 0x9a, 0x67, 0x96, 0x28, 0x9d, 
        0x38, 0xf9, 0x69, 0xf3, 0xbc, 0xa9, 0x76, 0x94, 
        0xb8, 0x8d, 0x5b, 0x8f, 0xff, 0x27, 0x5f, 0xde, 
        0xf8, 0xef, 0x6d, 0x32, 0x07, 0x2e, 0x36, 0x29, 
        0x94, 0x64, 0x54, 0xd3, 0x62, 0xf0, 0xee, 0x3d, 
        0x24, 0x59, 0xc9, 0x89, 0x8f, 0x41, 0x76, 0xe5, 
        0x62, 0xea, 0x45, 0xe1, 0xe7, 0x21, 0xfd, 0x48, 
        0x58, 0x2b, 0xba, 0x00, 0xd3, 0xb4, 0x46, 0x93, 
        0xc4, 0x85, 0x41, 0x91, 0x9c, 0x8a, 0x1a, 0x4c, 
        0x74, 0x9b, 0x70, 0xc1, 0x26, 0x2c, 0xb8, 0x2e, 
        0xc8, 0x39, 0x62, 0x79, 0x88, 0xd0, 0x73, 0xdf, 
        0x2a, 0x6c, 0x3c, 0x9b, 0xcb, 0xf4, 0x12, 0x79, 
        0xce, 0xe8, 0x74, 0xe7, 0xb7, 0x03, 0x5c, 0x7a, 
        0xc0, 0x63, 0x50, 0x27, 0xf2, 0xa7, 0x3e, 0xa7, 
        0xaf, 0x4c, 0xfa, 0x53, 0x8f, 0xcc, 0x83, 0x12, 
        0x3f, 0x03, 0x89, 0x85, 0x2a, 0xbb, 0x3d, 0xac, 
        0x31, 0xf4, 0x9f, 0x4a, 0xb7, 0x4c, 0xe8, 0x27, 
        0x80, 0xc8, 0x24, 0x4d, 0x5e, 0xa8, 0xe9, 0xc6, 
        0x03, 0x4c, 0xb0, 0xaa, 0x84, 0xe3, 0x79, 0x74, 
        0x3e, 0xf4, 0x7f, 0x8f, 0x16, 0xd7, 0x5c, 0x4d, 
        0xb1, 0xe1, 0xfe, 0x01, 0xfa, 0x6a, 0x9c, 0x49, 
        0x3f, 0x56, 0x82, 0x4c, 0x0b, 0x8b, 0x8d, 0x4f, 
        0x52, 0xeb, 0x8d, 0x44, 0xde, 0xb5, 0x7f, 0xcd, 
        0x55, 0x45, 0x44, 0x35, 0x2b, 0x66, 0x55, 0xb9, 
        0x73, 0x6b, 0x1f, 0x3f, 0xdf, 0x55, 0x03, 0x99, 
        0x10, 0xe0, 0xac, 0x94, 0x9e, 0xeb, 0x5e, 0x25, 
        0x96, 0xf5, 0xe1, 0x3e, 0xbe, 0x6a, 0x45, 0xb4, 
        0x6e, 0xb7, 0x20, 0xd9, 0xf4, 0x9b, 0x4b, 0x40, 
        0x80, 0xf6, 0x49, 0xd6, 0xc9, 0x30, 0x27, 0xa5, 
        0x95, 0x64, 0x0e, 0x5c, 0x61, 0xc7, 0x6b, 0xb8, 
        0x18, 0x65, 0x19, 0x18, 0xde, 0x77, 0x81, 0x3f, 
        0x2a, 0xc2, 0x2f, 0x72, 0xf1, 0xc7, 0x6b, 0x34, 
        0x32, 0xcf, 0x32, 0xf7, 0x8c, 0xb8, 0xb6, 0x7d, 
        0xaa, 0x17, 0x7c, 0x96, 0xb5, 0x99, 0xa3, 0x2f, 
        0xd1, 0xea, 0xc7, 0xcd, 0x53, 0xdf, 0x4a, 0x06, 
        0x30, 0x35, 0xc4, 0x15, 0x16, 0x95, 0xa2, 0x7a, 
        0x2d, 0xb9, 0xab, 0x74, 0x25, 0xe5, 0x62, 0x1b, 
        0x28, 0x62, 0xb2, 0x9d, 0x9c, 0x94, 0xa8, 0x7b, 
        0xeb, 0xb4, 0x08, 0xa6, 0x07, 0x1f, 0x29, 0xb2, 
        0xde, 0x71, 0x25, 0x0f, 0x99, 0x7f, 0x2e, 0xaf, 
        0xff, 0xac, 0xd9, 0x02, 0xdd, 0xac, 0xc6, 0x74, 
        0xb4, 0x67, 0xe3, 0x99, 0x7d, 0x2d, 0xd4, 0xce, 
        0xe0, 0x31, 0xfe, 0xea, 0x03, 0x8f, 0x3a, 0x51, 
        0x89, 0x5b, 0x44, 0xb4, 0x86, 0x03, 0xcd, 0x5c, 
        0xf0, 0x4e, 0x12, 0x39, 0x5c, 0x57, 0xd6, 0xa5, 
        0x9c, 0xd9, 0x8f, 0x3a, 0x44, 0xb6, 0x97, 0xeb, 
        0xaf, 0x1d, 0xfe, 0x96, 0xae, 0x27, 0x6c, 0x8c, 
        0xa6, 0x57, 0x4d, 0x84, 0xce, 0x93, 0xe5, 0x5a, 
        0x59, 0xe6, 0xe1, 0x10, 0x66, 0xfd, 0x6b, 0x4e, 
        0x6a, 0xd0, 0x81, 0x7b, 0xda, 0x8e, 0x81, 0x06, 
        0x37, 0x38, 0x5f, 0x87, 0xb7, 0x66, 0x43, 0xdf, 
        0x0f, 0x76, 0x1e, 0x2d, 0x73, 0xe6, 0x06, 0xa2, 
        0x05, 0x8f, 0x7c, 0x26, 0x17, 0x02, 0x2f, 0x4f, 
        0xf3, 0x18, 0x4c, 0x32, 0xa8, 0x05, 0xd6, 0xda, 
        0x9c, 0xc5, 0x80, 0x6d, 0x52, 0xf1, 0x0e, 0x9e, 
        0x3f, 0x29, 0xd7, 0x9b, 0xf5, 0xb9, 0x60, 0xbe, 
        0xdb, 0x83, 0xf2, 0x25, 0x42, 0x2c, 0x2e, 0x08, 
        0x4b, 0x72, 0x57, 0xfd, 0xd9, 0xf6, 0xcb, 0x4f, 
        0x14, 0xe3, 0xb3, 0x25, 0x21, 0x50, 0x64, 0x6b, 
        0xc5, 0xcc, 0x27, 0x21, 0x47, 0x30, 0xf9, 0xab, 
        0x10, 0x46, 0xc4, 0x59, 0x16, 0xae, 0xc3, 0x43, 
        0x18, 0x23, 0x7c, 0x0c, 0x56, 0x84, 0x61, 0x1d, 
        0xc9, 0xbf, 0xa2, 0xba, 0x50, 0x4e, 0x4b, 0x46, 
        0xd6, 0x3e, 0x21, 0x64, 0x74, 0xd8, 0x2b, 0x94, 
        0x48, 0x4e, 0x7e, 0xf3, 0xdd, 0x9f, 0x83, 0xe3, 
        0x30, 0xe7, 0xa5, 0xad, 0x73, 0x70, 0xe4, 0x6c, 
        0xca, 0x16, 0x7b, 0x09, 0x53, 0x51, 0xf7, 0x59, 
        0xea, 0x0e, 0x77, 0xf6, 0x75, 0x0c, 0x8d, 0xa8, 
        0xa4, 0x82, 0xa8, 0x9d, 0x41, 0x3a, 0x7a, 0xf2, 
        0xd8, 0xf5, 0x5b, 0xb6, 0xab, 0xea, 0x64, 0xba, 
        0xf0, 0x07, 0x54, 0x1b, 0xfa, 0xba, 0x1f, 0xd2, 
        0xc7, 0x32, 0x9f, 0x4e, 0x23, 0x0b, 0x5b, 0xe7, 
        0x9f, 0x86, 0x4d, 0xec, 0x31, 0x02, 0x91, 0xef, 
        0x1f, 0x8f, 0x41, 0x4b, 0x40, 0xa4, 0x3b, 0x44, 
        0x0f, 0x93, 0xae, 0xcd, 0xb6, 0x5e, 0xd2, 0x30, 
        0x16, 0x86, 0x42, 0x75, 0x33, 0x00, 0x9e, 0x24, 
        0xfc, 0x55, 0xea, 0xca, 0xfb, 0x39, 0xda, 0x1b, 
        0x39, 0xff, 0x9c, 0x3b, 0xa9, 0x63, 0xcf, 0x24, 
        0x5f, 0xe4
    };
    const unsigned char t[TAGLEN] = {
        0x6d, 0xd8, 0x86, 0xf8, 0x8b, 0x9f, 0xc6, 0xa0, 
        0xc3, 0x3d, 0x2e, 0x3a, 0x84, 0x81, 0xf4, 0x99
    };
    return run_test((const unsigned char*)k, (const unsigned char*)h, hlen,
                    (unsigned char*)m, mlen, c, t);
}

// ---------------------------------------------------------------------

int main()
{
    int result = 0;

    result |= test1();
    result |= test2();
    result |= test3();
    result |= test4();
    result |= test5();
    // result |= flexible_size_tests(kat);

    if (result) {
        puts("Test result:  FAILED");
    } else {
        puts("Tests result: SUCCESS");
    }

    return 0;
}
