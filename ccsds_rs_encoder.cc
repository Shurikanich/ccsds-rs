// Standalone CCSDS Encoder (GNU Radio dependencies removed)

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "ccsds.h"
#include "reed_solomon.h"
#include "ccsds_rs_encoder.h"


ccsds_rs_encoder::ccsds_rs_encoder(bool rs_encode, bool interleave, bool scramble,
                             bool printing, bool verbose, int n_interleave, bool dual_basis)
    : d_rs_encode(rs_encode), d_interleave(interleave), d_scramble(scramble),
      d_printing(printing), d_verbose(verbose),
      d_n_interleave(n_interleave), d_dual_basis(dual_basis)
{
    memcpy(d_pkt.sync_word, SYNC_WORD, SYNC_WORD_LEN);
}

int ccsds_rs_encoder::encode(const uint8_t* in, uint8_t* out)
{
    if (!in || !out) return 0;

    uint8_t rs_block[RS_BLOCK_LEN];
    for (uint8_t i = 0; i < d_n_interleave; i++)
    {
        if (d_interleave)
        {
            for (uint8_t j = 0; j < RS_BLOCK_LEN; j++)
            {
                rs_block[j] = in[i + (d_n_interleave * j)];
            }
        }
        else
        {
            memcpy(rs_block, &in[i * RS_DATA_LEN], RS_DATA_LEN);
        }

        if (d_rs_encode)
        {
            d_rs.encode(rs_block, d_dual_basis);
        }
        else
        {
            memset(&rs_block[RS_DATA_LEN], 0, RS_PARITY_LEN);
        }

        if (d_interleave)
        {
            for (uint8_t j = 0; j < RS_BLOCK_LEN; j++)
                d_pkt.codeword[i + (d_n_interleave * j)] = rs_block[j];
        }
        else
        {
            memcpy(&d_pkt.codeword[i * RS_BLOCK_LEN], rs_block, RS_BLOCK_LEN);
        }
    }

    if (d_scramble)
    {
        scramble(d_pkt.codeword, codeword_len());
    }

    d_num_frames++;
    if (d_verbose)
    {
        printf("sending %i bytes of data\n", total_frame_len());
        printf("number of frames transmitted: %i\n", d_num_frames);
    }

    if (d_printing)
    {
        print_bytes(d_pkt.codeword, codeword_len());
    }

    memcpy(out, &d_pkt, total_frame_len());
    return total_frame_len();
}