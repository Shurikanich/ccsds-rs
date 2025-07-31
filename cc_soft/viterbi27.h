/* Copyright 1994 Phil Karn, KA9Q
 * May be used under the terms of the GNU Public License
 */

#ifndef __VITERBI27_H__
#define __VITERBI27_H__

#undef DEBUG

#ifdef DEBUG
    #define DEBUG_PRINT(...) printf(__VA_ARGS__)
#else
    #define DEBUG_PRINT(...)
#endif

/* The two generator polynomials for the NASA Standard K=7 rate 1/2 code. */
#define	POLYA	0x6d
#define	POLYB	0x4f


/* This parameter sizes the path memory in bits, which is organized as a
 * circular buffer through which we periodically "trace back" to
 * produce the decoded data. PATHMEM must be greater than
 * MERGEDIST+TRACECHUNK, and for efficiency it should also be a power of 2.
 * Don't make it *too* large, or it will spill out of the CPU's on-chip cache
 * and decrease performance. Each bit of path memory costs 8 bytes for the
 * K=7 code.
 */
#define PATHMEM	256

/* In theory, a Viterbi decoder is true maximum likelihood only if
 * the path memory is as long as the entire message and a single traceback
 * is made from the terminal state (usually zero) after the entire message
 * is received.
 *
 * In practice, performance is essentially optimum as long as decoding
 * decisions are deferred by at least 4-5 constraint lengths (28-35 bits
 * for K=7) from the most recently received symbols. MERGEDIST sets this
 * parameter. We give ourselves some margin here in case the code is
 * punctured (which slows merging) and also to let us start each traceback
 * from an arbitrary current state instead of taking the time to find the
 * path with the highest current metric.
 */
#define	MERGEDIST	128	/* Distance to trace back before decoding */

/* Since each traceback is costly (thanks to the overhead of having to
 * go back MERGEDIST bits before we produce our first decoded bit) we'd like
 * to decode as many bits as possible per traceback at the expense of
 * increased decoding delay. TRACECHUNK sets how many bits to
 * decode on each traceback. Since output is produced in 8-bit bytes,
 * TRACECHUNK MUST be a multiple of 8.
 */
#define	TRACECHUNK	8	/* How many bits to decode on each traceback */

/* The path metrics need to be periodicially adjusted downward
 * to prevent an integer overflow that could cause the signed comparisons
 * in the butterfly macros to fail.
 *
 * It's possible to code the comparisons to work in modulo fashion, e.g.,
 * as 'if((a-b) > 0)' rather than 'if(a >b)'. A good optimizer would generate
 * code like 'cmp a,b;js foo' for this, but GCC doesn't.
 *
 * This constant should be larger than the maximum path metric spread.
 * Experimentally this seems to be 2040, which is probably related to the
 * free distance of the code (10) and the symbol metric scale (0-255).
 */
#define	RENORMALIZE	10000

#if (TRACECHUNK + MERGEDIST > PATHMEM)
#error "TRACECHUNK + MERGEDIST > PATHMEM"
#endif

#if ((TRACECHUNK % 8) != 0)
#error "TRACECHUNK not multiple of 8"
#endif


// code rate = 1/2
#define CODE_RATE_12 (1.0/2.0)
#define PUNCTURE_PATTERN_LEN_12 1  // Length of the puncturing pattern
static const int puncture_C1_12[PUNCTURE_PATTERN_LEN_12] = {1};  // C1 puncturing
static const int puncture_C2_12[PUNCTURE_PATTERN_LEN_12] = {1};  // C2 puncturing

// code rate = 3/4
#define CODE_RATE_34 (3.0/4.0)
#define PUNCTURE_PATTERN_LEN_34 3  // Length of the puncturing pattern
static const int puncture_C1_34[PUNCTURE_PATTERN_LEN_34] = {1, 0, 1};  // C1 puncturing
static const int puncture_C2_34[PUNCTURE_PATTERN_LEN_34] = {1, 1, 0};  // C2 puncturing

// code rate = 7/8
#define CODE_RATE_78 (7.0/8.0)
#define PUNCTURE_PATTERN_LEN_78 7  // Length of the puncturing pattern
static const int puncture_C1_78[PUNCTURE_PATTERN_LEN_78] = {1, 0, 0, 0, 1, 0, 1};  // C1 puncturing
static const int puncture_C2_78[PUNCTURE_PATTERN_LEN_78] = {1, 1, 1, 1, 0, 1, 0};  // C2 puncturing

// code rate = 2/3
#define CODE_RATE_23 (2.0/3.0)
#define PUNCTURE_PATTERN_LEN_23 2  // Length of the puncturing pattern
static const int puncture_C1_23[PUNCTURE_PATTERN_LEN_23] = {1, 0};  // C1 puncturing
static const int puncture_C2_23[PUNCTURE_PATTERN_LEN_23] = {1, 1};  // C2 puncturing

// code rate = 5/6
#define CODE_RATE_56 (5.0/6.0)
#define PUNCTURE_PATTERN_LEN_56 5  // Length of the puncturing pattern
static const int puncture_C1_56[PUNCTURE_PATTERN_LEN_56] = {1, 0, 1, 0, 1};  // C1 puncturing
static const int puncture_C2_56[PUNCTURE_PATTERN_LEN_56] = {1, 1, 0, 1, 0};  // C2 puncturing





typedef struct v27
{
    long cmetric[64];
    long nmetric[64];
    unsigned long paths[2*PATHMEM];
    unsigned int pi;
    unsigned long dec;
    int mets[4];
} v27;

#ifdef __cplusplus
extern "C" {
#endif

unsigned int encode27(unsigned char *encstate,
    unsigned char *symbols,
    unsigned char *data,
    unsigned int nbytes,
    const int* puncture_C1_ptr,
    const int* puncture_C2_ptr,
    int puncture_pattern_len);
void vitfilt27_init(v27 *vi);
void vitfilt27_decode(v27 *vi, unsigned char *syms, unsigned char *data, unsigned int nbits);
void encode27_bit(unsigned char *encstate, unsigned char *symbols, unsigned char *data);

#ifdef __cplusplus
}
#endif

#endif

