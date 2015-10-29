#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "api.h"
#include "riv.h"

// ---------------------------------------------------------------------

#define _Is_X86_            1
#define HI_RES_CLK_OK
#define TIMER_SAMPLE_CNT    10000

#define MAX_BUFFER_LEN      16384
#define MIN_CONT_MSG_LEN    4
#define MAX_CONT_MSG_LEN    256
#define NUM_INTERVALS       16

// ---------------------------------------------------------------------

// Default
#ifndef FN_TO_BENCH
#define FN_TO_BENCH encrypt_final
#endif

// ---------------------------------------------------------------------

/**
 * Quicksort comparison
 */
static int compare_doubles(const void *aPtr, const void *bPtr)
{
    double a = *((double*) aPtr);
    double b = *((double*) bPtr);
    if (a > b) return  1;
    if (a < b) return -1;
    return 0;
}

// ---------------------------------------------------------------------

static uint64_t get_time()
{
    uint64_t x[2];
    __asm__ volatile("rdtsc" : "=a"(x[0]), "=d"(x[1]));
    return x[0];
}

// ---------------------------------------------------------------------

static uint64_t calibrate_timer()
{
    // big number to start
    uint64_t dtMin = 0xFFFFFFFFL; 
    uint64_t t0, t1;
    int i;

    for (i = 0; i < TIMER_SAMPLE_CNT; ++i) {
        t0 = get_time();
        t1 = get_time();

        if (dtMin > t1 - t0) {
            dtMin = t1 - t0;
        }
    }

    return dtMin;
}

// ---------------------------------------------------------------------

static void prepare_k(uint8_t* k, const size_t klen)
{
    for (size_t i = 0; i < klen; ++i) {
        k[i] = (23 * i) & 0xFF;
    }
}

// ---------------------------------------------------------------------

static void prepare_m(uint8_t* m, const size_t mlen)
{
    for (size_t i = 0; i < mlen; ++i) {
        m[i] = (123 * i) & 0xFF;
    }
}

// ---------------------------------------------------------------------

static int benchmark(const uint32_t num_iterations)
{
    ALIGN(16) uint8_t k[KEYLEN];
    ALIGN(16) uint8_t* m = (uint8_t*)malloc(MAX_BUFFER_LEN);
    ALIGN(16) uint8_t* c = (uint8_t*)malloc(MAX_BUFFER_LEN);
    ALIGN(16) uint8_t t[TAGLEN];

    uint64_t mlen;
    const uint64_t hlen = 0L;

    const uint64_t calibration = calibrate_timer();
    uint64_t t0, t1;
    uint32_t i, j;

    prepare_k(k, KEYLEN);
    prepare_m(m, MAX_BUFFER_LEN);

    riv_context_t ctx;
    keysetup(&ctx, k);

    puts("#mlen cpb");

    double timings[num_iterations];
    const uint32_t median = num_iterations / 2;
    // To load the timing code into the instruction cache
    get_time(); 

    for (j = MIN_CONT_MSG_LEN; j < MAX_CONT_MSG_LEN; ++j) {
        mlen = j;
        
        t0 = get_time();
        t1 = get_time();

        for (i = 0; i < num_iterations; ++i) {
            t0 = get_time();
            
            FN_TO_BENCH(&ctx, m, mlen, NULL, hlen, c, t);

            t1 = get_time();
            timings[i] = (double)(t1 - t0 - calibration) / j;
        }

        // Sort the measurements and print the median
        qsort(timings, num_iterations, sizeof(double), compare_doubles);
        printf("%5d %4.3lf \n", j, timings[median]);
    }

    uint32_t interval_len;
    uint32_t j_outer;

    for (j_outer = MAX_CONT_MSG_LEN; j_outer < MAX_BUFFER_LEN; 
        j_outer = 2*j_outer) {
        
        interval_len = j_outer / NUM_INTERVALS;

        for (j = j_outer; j < 2*j_outer; j += interval_len) {
            mlen = j;
            
            t0 = get_time();
            t1 = get_time();

            for (i = 0; i < num_iterations; ++i) {
                t0 = get_time();
                
                FN_TO_BENCH(&ctx, m, mlen, NULL, hlen, c, t);

                t1 = get_time();
                timings[i] = (double)(t1 - t0 - calibration) / j;
            }

            // Sort the measurements and print the median
            qsort(timings, num_iterations, sizeof(double), compare_doubles);
            printf("%5d %4.3lf \n", j, timings[median]);
        }
    }


    t0 = get_time();
    t1 = get_time();
    mlen = MAX_BUFFER_LEN;

    for (i = 0; i < num_iterations; ++i) {
        t0 = get_time();
        
        FN_TO_BENCH(&ctx, m, mlen, NULL, hlen, c, t);

        t1 = get_time();
        timings[i] = (double)(t1 - t0 - calibration) / mlen;
    }

    // Sort the measurements and print the median
    qsort(timings, num_iterations, sizeof(double), compare_doubles);
    printf("%5lu %4.3lf \n", mlen, timings[median]);

    free(c);
    free(m);
    return 0;
}

// ---------------------------------------------------------------------

int main()
{
    return benchmark(TIMER_SAMPLE_CNT);
}
