#ifndef CCSDS_RS_DECODER_H
#define CCSDS_RS_DECODER_H

#include <stdint.h>
#include "reed_solomon.h"
#include "ccsds.h"

class ccsds_rs_decoder {
public:
    ccsds_rs_decoder(int threshold,
                     bool rs_decode,
                     bool deinterleave,
                     bool descramble,
                     bool verbose,
                     bool printing,
                     int n_interleave,
                     bool dual_basis);
    ~ccsds_rs_decoder() = default;

    int process(const uint8_t* in, int ninput_items, const uint8_t* out, int* noutput_items);

    uint32_t num_frames_received() const { return d_num_frames_received; }
    uint32_t num_frames_decoded()  const { return d_num_frames_decoded; }
    uint32_t num_subframes_decoded() const { return d_num_subframes_decoded; }

private:
    void enter_sync_search();
    void enter_codeword();
    bool compare_sync_word();
    bool decode_frame();

    inline int data_len() const { return RS_DATA_LEN * d_n_interleave; }
    inline int codeword_len() const { return RS_BLOCK_LEN * d_n_interleave; }

    int d_threshold;
    bool d_rs_decode;
    bool d_deinterleave;
    bool d_descramble;
    bool d_verbose;
    bool d_printing;
    int d_n_interleave;
    bool d_dual_basis;

    uint32_t d_data_reg = 0;
    uint32_t d_sync_word = 0;
    int d_decoder_state = 0;
    int d_bit_counter = 0;
    int d_byte_counter = 0;


    uint8_t d_codeword[CODEWORD_MAX_LEN] = {0};
    uint8_t d_payload[DATA_MAX_LEN] = {0};

    uint32_t d_num_frames_received = 0;
    uint32_t d_num_frames_decoded = 0;
    uint32_t d_num_subframes_decoded = 0;

    reed_solomon d_rs;
};

#endif // CCSDS_RS_DECODER_H