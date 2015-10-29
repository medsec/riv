#include "api.h"
#include "riv.h"

// ---------------------------------------------------------------------

void crypto_aead_encrypt(unsigned char *c, size_t *clen,
                         const unsigned char *h, size_t hlen,
                         const unsigned char *m, size_t mlen,
                         const unsigned char *n,
                         const unsigned char *k)
{
    riv_context_t ctx;
    (void)n;
    keysetup(&ctx, k);
    *clen = mlen + TAGLEN;
    encrypt_final(&ctx, m, mlen, h, hlen, c, c+mlen);
}

// ---------------------------------------------------------------------

int crypto_aead_decrypt(unsigned char *m, size_t *mlen,
                        const unsigned char *h, size_t hlen,
                        const unsigned char *c, size_t clen,
                        const unsigned char *n,
                        const unsigned char *k)
{
    riv_context_t ctx;
    (void)n;
    keysetup(&ctx, k);
    *mlen = clen - TAGLEN;
    return decrypt_final(&ctx, c, clen, h, hlen, c+(*mlen), m);
}
