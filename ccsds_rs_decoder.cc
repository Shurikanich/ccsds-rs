// Standalone CCSDS Reed-Solomon Decoder (GNU Radio dependencies removed)

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "reed_solomon.h"
#include "ccsds.h"
#include "ccsds_rs_decoder.h"

#define STATE_SYNC_SEARCH 0
#define STATE_CODEWORD 1



ccsds_rs_decoder::ccsds_rs_decoder(int threshold,
                                  bool rs_decode,
                                  bool deinterleave,
                                  bool descramble,
                                  bool verbose,
                                  bool printing,
                                  int n_interleave,
                                  bool dual_basis)
    : d_threshold(threshold), d_rs_decode(rs_decode), d_deinterleave(deinterleave), d_descramble(descramble),
      d_verbose(verbose), d_printing(printing), d_n_interleave(n_interleave), d_dual_basis(dual_basis)
{
    for (uint8_t i = 0; i < SYNC_WORD_LEN; i++)
    {
        d_sync_word = (d_sync_word << 8) | (SYNC_WORD[i] & 0xff);
    }
    enter_sync_search();
}

int ccsds_rs_decoder::process(const uint8_t* in, int ninput_items, const uint8_t *out, int *noutput_items)
{
    uint16_t count = 0;
    while (count < ninput_items)
    {
        switch (d_decoder_state)
        {
            case STATE_SYNC_SEARCH:
            {
                d_data_reg = (d_data_reg << 1) | (in[count++] & 0x01);
                if (compare_sync_word())
                {
                    if (d_verbose) printf("\tsync word detected\n");
                    d_num_frames_received++;
                    enter_codeword();
                }
                break;
            }
            case STATE_CODEWORD:
            {
                d_data_reg = (d_data_reg << 1) | (in[count++] & 0x01);
                d_bit_counter++;
                if (d_bit_counter == 8)
                {
                    d_codeword[d_byte_counter] = d_data_reg;
                    d_byte_counter++;
                    d_bit_counter = 0;
                }
                if (d_byte_counter == codeword_len())
                {
                    if (d_verbose) printf("\tloaded codeword of length %i\n", codeword_len());
                    if (d_printing) print_bytes(d_codeword, codeword_len());

                    bool success = decode_frame();

                    if (success)
                    {
                        out = d_payload;
                        *noutput_items = data_len();
                        if (d_verbose)
                        {
                            printf("\tframes received: %i\n\tframes decoded: %i\n\tsubframes decoded: %i\n", d_num_frames_received, d_num_frames_decoded, d_num_subframes_decoded);
                        }
                    }
                    enter_sync_search();
                }
                break;
            }
        }
    }
    return ninput_items;
}

void ccsds_rs_decoder::enter_sync_search()
{
    if (d_verbose) printf("enter sync search\n");
    d_decoder_state = STATE_SYNC_SEARCH;
    d_data_reg = 0;
}

void ccsds_rs_decoder::enter_codeword()
{
    if (d_verbose) printf("enter codeword\n");
    d_decoder_state = STATE_CODEWORD;
    d_byte_counter = 0;
    d_bit_counter = 0;
}

bool ccsds_rs_decoder::compare_sync_word()
{
    uint32_t wrong_bits = d_data_reg ^ d_sync_word;
    int nwrong = __builtin_popcount(wrong_bits);
    return nwrong <= d_threshold;
}

bool ccsds_rs_decoder::decode_frame()
{
    bool success = true;

    if (d_descramble)
    {
        descramble(d_codeword, codeword_len());
    }

    uint8_t rs_block[RS_BLOCK_LEN];
    int8_t nerrors;
    for (uint8_t i = 0; i < d_n_interleave; i++)
    {
        for (uint8_t j = 0; j < RS_BLOCK_LEN; j++)
        {
            if (d_deinterleave)
            {
                rs_block[j] = d_codeword[i + (j * d_n_interleave)];
            }
            else
            {
                rs_block[j] = d_codeword[i * RS_BLOCK_LEN + j];
            }
        }
        if (d_rs_decode)
        {
            nerrors = d_rs.decode(rs_block, d_dual_basis);
            if (nerrors == -1)
            {
                if (d_verbose) printf("\tcould not decode rs block #%i\n", i);
                success = false;
            }
            else
            {
                if (d_verbose) printf("\tdecoded rs block #%i with %i errors\n", i, nerrors);
                d_num_subframes_decoded++;
            }
        }
        if (d_deinterleave)
        {
            for (uint8_t j = 0; j < RS_DATA_LEN; j++)
            {
                d_payload[i + (j * d_n_interleave)] = rs_block[j];
            }
        }
        else
        {
            memcpy(&d_payload[i * RS_DATA_LEN], rs_block, RS_DATA_LEN);
        }
    }

    if (success) d_num_frames_decoded++;

    return success;
}