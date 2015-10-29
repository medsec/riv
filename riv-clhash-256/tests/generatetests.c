#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "riv.h"
#include "api.h"

// ---------------------------------------------------------------------

extern 
void crypto_aead_encrypt(unsigned char *c, size_t *clen,
                         const unsigned char *h, size_t hlen,
                         const unsigned char *m, size_t mlen,
                         const unsigned char *n,
                         const unsigned char *k);

// ---------------------------------------------------------------------

static void print_hex(const unsigned char *x, 
                      const long long unsigned len)
{
    for (unsigned long long i = 0; i < len; i++) {
        if ((i != 0) && (i % 8 == 0)) {
            puts("");
        }

        printf("0x%02x, ", x[i]);
    }
}

// ---------------------------------------------------------------------

static void genkat()
{
	static const size_t MAX_H_LEN = 1024;
	static const size_t MAX_M_LEN = 4096;

	static const uint32_t H_LENS[11] = {0, 1, 2, 8, 16, 32, 64, 128, 256, 512, 1024};
	static const uint32_t M_LENS[15] = {0, 1, 2, 8, 16, 32, 64, 128, 256, 512, 1023, 1024, 1025, 2048, 4096};

	unsigned char* w = (uint8_t*)malloc(MAX_M_LEN);
	unsigned char* h = (uint8_t*)malloc(MAX_H_LEN);
	unsigned char* m = (unsigned char*)malloc(MAX_M_LEN);
	unsigned char* c = (unsigned char*)malloc(MAX_M_LEN + TAGLEN);
	unsigned char k[KEYLEN];

	size_t i, j;

	for(i = 0; i < MAX_H_LEN; ++i) {
		h[i] = 0xFF & (i*193 + 123);
	}

	for(i = 0; i < KEYLEN; ++i) {
		k[i] = 0xFF & (i*191 + 123);
	}

	for(i = 0; i < MAX_M_LEN; ++i) {
		w[i] = 0xFF & (i*181 + 123);
	}

	puts("#pragma once\n");
	puts("static const unsigned char kat[] = {");

	for(i = 0; i < 11; ++i) {
		size_t hlen = H_LENS[i];

		for(j = 0; j < 15; ++j) {
			size_t mlen = M_LENS[j];
			size_t clen = 0;

			memcpy(m, w, mlen);
			memset(c, 0x00, mlen+TAGLEN);

			crypto_aead_encrypt(c, &clen, h, hlen, m, mlen, NULL, k);
			print_hex(c, clen);
			puts("\n");
		}
	}

	puts("};");

	free(w);
	free(h);
	free(m);
	free(c);
}

// ---------------------------------------------------------------------

int main()
{
	genkat();
	return 0;
}

