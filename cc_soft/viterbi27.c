/* Viterbi decoder for K=7 rate=1/2 convolutional code
 * continuous traceback version
 * Copyright 1996 Phil Karn, KA9Q
 *
 * This version of the Viterbi decoder reads a continous stream of
 * 8-bit soft decision samples from standard input in offset-binary
 * form, i.e., a 255 sample is the strongest possible "1" symbol and a
 * 0 is the strongest possible "0" symbol. 128 is an erasure (unknown).
 *
 * The decoded output is written to stdout in big-endian form (the first
 * decoded bit appears in the high order bit of the first output byte).
 *
 * The metric table is fixed, and no attempt is made (yet) to find proper
 * symbol synchronization. These are likely future enhancements.
 */
//#include <stdio.h>
#include <limits.h>
#include "viterbi27.h"
#include <stdio.h>




/* The basic Viterbi decoder operation, called a "butterfly"
 * operation because of the way it looks on a trellis diagram. Each
 * butterfly involves an Add-Compare-Select (ACS) operation on the two nodes
 * where the 0 and 1 paths from the current node merge at the next step of
 * the trellis.
 *
 * The code polynomials are assumed to have 1's on both ends. Given a
 * function encode_state() that returns the two symbols for a given
 * encoder state in the low two bits, such a code will have the following
 * identities for even 'n' < 64:
 *
 * 	encode_state(n) = encode_state(n+65)
 *	encode_state(n+1) = encode_state(n+64) = (3 ^ encode_state(n))
 *
 * Any convolutional code you would actually want to use will have
 * these properties, so these assumptions aren't too limiting.
 *
 * Doing this as a macro lets the compiler evaluate at compile time the
 * many expressions that depend on the loop index and encoder state and
 * emit them as immediate arguments.
 * This makes an enormous difference on register-starved machines such
 * as the Intel x86 family where evaluating these expressions at runtime
 * would spill over into memory.
 *
 * Two versions of the butterfly are defined. The first reads cmetric[]
 * and writes nmetric[], while the other does the reverse. This allows the
 * main decoding loop to be unrolled to two bits per loop, avoiding the
 * need to reference the metrics through pointers that are swapped at the
 * end of each bit. This was another performance win on the register-starved
 * Intel CPU architecture.
 */

#define	BUTTERFLY(i,sym) { \
	long m0,m1;\
	/* ACS for 0 branch */\
    DEBUG_PRINT("i:%d sym:%d", i, sym);\
	m0 = vi->cmetric[i] + vi->mets[sym];	/* 2*i */\
	m1 = vi->cmetric[i+32] + vi->mets[3^sym];	/* 2*i + 64 */\
    DEBUG_PRINT(" m0:%ld m1:%ld", m0, m1);\
	vi->nmetric[2*i] = m0;\
	if(m1 > m0){\
		vi->nmetric[2*i] = m1;\
		vi->dec |= 1 << ((2*i) & 31);\
	}\
    DEBUG_PRINT(" cmetric[%d]:%ld dec:%lx val:%ld @%d",2*i,vi->cmetric[2*i], vi->dec, (vi->dec >> ((2*i) & 31)), ((2*i) & 31));\
    DEBUG_PRINT("\n");\
	/* ACS for 1 branch */\
	m0 -= (vi->mets[sym] - vi->mets[3^sym]);\
	m1 += (vi->mets[sym] - vi->mets[3^sym]);\
    DEBUG_PRINT("m0:%ld m1:%ld", m0, m1);\
	vi->nmetric[2*i+1] = m0;\
	if(m1 > m0){\
		vi->nmetric[2*i+1] = m1;\
		vi->dec |= 1 << ((2*i+1) & 31);\
	}\
    DEBUG_PRINT(" cmetric[%d]:%ld dec:%lx val:%ld @%d",2*i+1,vi->cmetric[2*i+1], vi->dec, (vi->dec >> ((2*i+1) & 31)), ((2*i+1) & 31));\
    DEBUG_PRINT("\n");\
}

#define	BUTTERFLY2(i,sym) { \
	long m0,m1;\
	/* ACS for 0 branch */\
    DEBUG_PRINT("(2)i:%d sym:%d", i, sym);\
	m0 = vi->nmetric[i] + vi->mets[sym];	/* 2*i */\
	m1 = vi->nmetric[i+32] + vi->mets[3^sym]; /* 2*i + 64 */\
    DEBUG_PRINT(" m0:%ld m1:%ld", m0, m1);\
	vi->cmetric[2*i] = m0;\
	if(m1 > m0){\
		vi->cmetric[2*i] = m1;\
		vi->dec |= 1 << ((2*i) & 31);\
	}\
    DEBUG_PRINT(" cmetric[%d]:%ld dec:%lx val:%ld @%d",2*i,vi->cmetric[2*i], vi->dec, (vi->dec >> ((2*i) & 31)), ((2*i) & 31));\
    DEBUG_PRINT("\n");\
	/* ACS for 1 branch */\
	m0 -= (vi->mets[sym] - vi->mets[3^sym]);\
	m1 += (vi->mets[sym] - vi->mets[3^sym]);\
    DEBUG_PRINT("m0:%ld m1:%ld", m0, m1);\
	vi->cmetric[2*i+1] = m0;\
	if(m1 > m0){\
		vi->cmetric[2*i+1] = m1;\
		vi->dec |= 1 << ((2*i+1) & 31);\
	}\
    DEBUG_PRINT(" cmetric[%d]:%ld dec:%lx val:%ld @%d",2*i+1,vi->cmetric[2*i+1], vi->dec, (vi->dec >> ((2*i+1) & 31)), ((2*i+1) & 31));\
    DEBUG_PRINT("\n");\
}



int mettab[2][256];

void vitfilt27_init(v27 *vi)
{
    int i;

    /* Initialize metric table (make this an option)
     * This table assumes a symbol of 0 is the
     * strongest possible '0', and a symbol
     * of 255 is the strongest possible '1'. A symbol
     * of 128 is an erasure
     */
    for(i=0; i<256; i++)
    {
        mettab[0][i] = 128 - i;
        mettab[1][255-i] = 127 - i;
    }

    vi->cmetric[0] = 0;
    for(i=1; i<64; i++)
        vi->cmetric[i] = -99999;

    vi->pi = 0;
}

/* Periodic traceback to produce decoded data */
static void
traceback(unsigned long paths[],unsigned int pi, unsigned char *dst)
{
    int beststate,i,j;
    unsigned char data[TRACECHUNK/8];

    /* Start on an arbitrary path and trace it back until it's almost
     * certain we've merged onto the best path
     */
    beststate = 0;	/* arbitrary */
    pi = (pi - 1) % PATHMEM;	/* Undo last increment of pi */
    for(i=0; i < MERGEDIST-6; i++)
    {
        if(paths[2*pi + (beststate >> 5)] & (1 << (beststate & 31)))
        {
            beststate |= 64;	/* 2^(K-1) */
        }
        beststate >>= 1;
        pi = (pi - 1) % PATHMEM;
    }
    /* bestpath is now the encoder state on the best path, MERGEDIST
     * bits back. We continue to chain back until we accumulate
     * TRACECHUNK bits of decoded data
     */
    for(j=sizeof(data)-1; j >= 0; j--)
    {
        data[j] = 0;
        for(i=0; i<8; i++)
        {
            if(paths[2*pi + (beststate >> 5)] & (1 << (beststate & 31)))
            {
                beststate |= 64;	/* 2^(K-1) */
                data[j] |= 1 << i;
            }
            beststate >>= 1;
            pi = (pi - 1) % PATHMEM;
        }
    }

    for(i=0; i<(int)sizeof(data); i++)
    {
        DEBUG_PRINT("data[%d]:%02X\n", i, data[i]);
        dst[i] = data[i];
    } 
}

void vitfilt27_decode(v27 *vi, unsigned char *syms, unsigned char *data, unsigned int nbits)
{
    int i;
    unsigned char symbols[2];

#if ((nbits % (2*TRACECHUNK) ) != 0)
#error "nbits not multiple of 2*TRACECHUNK"
#endif

    /* Main loop -- read input symbols and run ACS butterflies,
     * periodically tracing back to produce decoded output data.
     * The loop is unrolled to process two bits per iteration.
     */
    while(nbits)
    {
        /* Renormalize metrics to prevent overflow */
        if(vi->cmetric[0] > (LONG_MAX - RENORMALIZE))
        {
            for(i=0; i<64; i++)
                vi->cmetric[i] -= LONG_MAX;
        }
        else if(vi->cmetric[0] < LONG_MIN+RENORMALIZE)
        {
            for(i=0; i<64; i++)
                vi->cmetric[i] += LONG_MAX;
        }
        /* Read input symbol pair and compute branch metrics */
        symbols[0] = *(syms++);
        symbols[1] = *(syms++);
        DEBUG_PRINT("symbols[0]:%d symbols[1]:%d\n", symbols[0], symbols[1]);
        nbits-=2;

        vi->mets[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
        vi->mets[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
        vi->mets[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
        vi->mets[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];

        DEBUG_PRINT("mets[0]:%d mets[1]:%d mets[2]:%d mets[3]:%d\n", vi->mets[0], vi->mets[1], vi->mets[2], vi->mets[3]);   

        /* On even numbered bits, the butterflies read from cmetrics[]
         * and write to nmetrics[]. On odd numbered bits, the reverse
         * is done
         */
        vi->dec = 0;
        BUTTERFLY(0,1);
        BUTTERFLY(1,3);
        BUTTERFLY(2,2);
        BUTTERFLY(3,0);
        BUTTERFLY(4,2);
        BUTTERFLY(5,0);
        BUTTERFLY(6,1);
        BUTTERFLY(7,3);
        BUTTERFLY(8,1);
        BUTTERFLY(9,3);
        BUTTERFLY(10,2);
        BUTTERFLY(11,0);
        BUTTERFLY(12,2);
        BUTTERFLY(13,0);
        BUTTERFLY(14,1);
        BUTTERFLY(15,3);
        vi->paths[2*vi->pi] = vi->dec;
        DEBUG_PRINT("dec:%lx\n", vi->dec);
        vi->dec = 0;
        BUTTERFLY(16,0);
        BUTTERFLY(17,2);
        BUTTERFLY(18,3);
        BUTTERFLY(19,1);
        BUTTERFLY(20,3);
        BUTTERFLY(21,1);
        BUTTERFLY(22,0);
        BUTTERFLY(23,2);
        BUTTERFLY(24,0);
        BUTTERFLY(25,2);
        BUTTERFLY(26,3);
        BUTTERFLY(27,1);
        BUTTERFLY(28,3);
        BUTTERFLY(29,1);
        BUTTERFLY(30,0);
        BUTTERFLY(31,2);
        vi->paths[2*vi->pi+1] = vi->dec;
        DEBUG_PRINT("dec:%lx\n", vi->dec);
        vi->pi++;

        /* Read input symbol pair and compute branch metrics */
        symbols[0] = *(syms++);
        symbols[1] = *(syms++);
        nbits-=2;

        vi->mets[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
        vi->mets[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
        vi->mets[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
        vi->mets[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];

        vi->dec = 0;
        BUTTERFLY2(0,1);
        BUTTERFLY2(1,3);
        BUTTERFLY2(2,2);
        BUTTERFLY2(3,0);
        BUTTERFLY2(4,2);
        BUTTERFLY2(5,0);
        BUTTERFLY2(6,1);
        BUTTERFLY2(7,3);
        BUTTERFLY2(8,1);
        BUTTERFLY2(9,3);
        BUTTERFLY2(10,2);
        BUTTERFLY2(11,0);
        BUTTERFLY2(12,2);
        BUTTERFLY2(13,0);
        BUTTERFLY2(14,1);
        BUTTERFLY2(15,3);
        vi->paths[2*vi->pi] = vi->dec;
        DEBUG_PRINT("dec:%lx\n", vi->dec);
        vi->dec = 0;
        BUTTERFLY2(16,0);
        BUTTERFLY2(17,2);
        BUTTERFLY2(18,3);
        BUTTERFLY2(19,1);
        BUTTERFLY2(20,3);
        BUTTERFLY2(21,1);
        BUTTERFLY2(22,0);
        BUTTERFLY2(23,2);
        BUTTERFLY2(24,0);
        BUTTERFLY2(25,2);
        BUTTERFLY2(26,3);
        BUTTERFLY2(27,1);
        BUTTERFLY2(28,3);
        BUTTERFLY2(29,1);
        BUTTERFLY2(30,0);
        BUTTERFLY2(31,2);
        vi->paths[2*vi->pi+1] = vi->dec;
        DEBUG_PRINT("dec:%lx\n", vi->dec);
        vi->pi = (vi->pi + 1) % PATHMEM;
        if((vi->pi % TRACECHUNK) == 0)
        {
            traceback(vi->paths, vi->pi, data);

            //printf("data: %d", data);
            data += TRACECHUNK/8;
        }
    }
}


extern unsigned char Partab[];	/* Parity lookup table */

unsigned int encode27(unsigned char *encstate,
                     unsigned char *symbols,
                     unsigned char *data,
                     unsigned int nbytes,
                     const int* puncture_C1_ptr,
                     const int* puncture_C2_ptr,
                     int puncture_pattern_len)
{
    unsigned char c;
    int i;
    // variable to track puncturing pattern
    int pattern_index = 0;
    unsigned int number_of_coded_symbols = 0;

    while(nbytes--)
    {
        c = *(data++);

        for(i=7; i>=0; i--)
        {
            DEBUG_PRINT("s%d", (*encstate) & ((1 << 6) - 1));
            (*encstate) = ((*encstate) << 1) | ((c >> 7) & 1);
            DEBUG_PRINT("->s%d", (*encstate) & ((1 << 6) - 1));
            DEBUG_PRINT(" :%d", ((c >> 7) & 1));
            c <<= 1;

            unsigned char s1 = Partab[(*encstate) & POLYB];  // First bit from C1
            unsigned char s2 = !Partab[(*encstate) & POLYA]; // Second bit from C2

            DEBUG_PRINT("%d%d", s1, s2);
            DEBUG_PRINT("\n");

            // Apply puncturing pattern
            if (puncture_C1_ptr[pattern_index])
            {
                *(symbols++) = s1;
                number_of_coded_symbols++;
            }
            
            if (puncture_C2_ptr[pattern_index])
            {
                *(symbols++) = s2;
                number_of_coded_symbols++;
            }
            
            // Cycle through puncturing pattern
            pattern_index = (pattern_index + 1) % puncture_pattern_len; 

            /* 1-sym -> 255, 0-sym -> 0 */
            //*(symbols++) = 0 - Partab[vi->encstate & POLYB];
            //*(symbols++) = 0 - !Partab[vi->encstate & POLYA];
        }
    }

#if 0
    // Append 8 tail bits
    for(i=0; i<8; i++)
    {
        (*encstate) <<= 1;
        *(symbols++) = Partab[(*encstate) & POLYB];
        *(symbols++) = !Partab[(*encstate) & POLYA];
    }
#endif

    return number_of_coded_symbols;
}

void encode27_bit(unsigned char *encstate, unsigned char *symbols, unsigned char *data)
{
    unsigned char c;
    c = data[0];
    (*encstate) = ((*encstate) << 1) | (c & 1);
   
    /* 1-sym -> 1, 0-sym -> 0 */
    *(symbols) = Partab[(*encstate) & POLYB];
    symbols++;
    *(symbols) = !Partab[(*encstate) & POLYA];
    symbols++;
    (*encstate) &= (1 << 6) - 1;
    return;
}


