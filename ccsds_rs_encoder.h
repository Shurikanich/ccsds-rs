#ifndef ccsds_rs_encoder_H
#define ccsds_rs_encoder_H

#include <stdint.h>
#include "reed_solomon.h"
#include "ccsds.h"

class ccsds_rs_encoder {
public:
    /**
     * Constructor for CCSDS encoder
     */
    ccsds_rs_encoder(bool rs_encode,
                  bool interleave,
                  bool scramble,
                  bool printing,
                  bool verbose,
                  int n_interleave,
                  bool dual_basis);

    ~ccsds_rs_encoder() = default;

    /**
     * Encodes an input payload into a CCSDS frame.
     *
     * @param in_payload   Pointer to input data (size = RS_DATA_LEN * n_interleave)
     * @param out_frame    Output buffer to hold the encoded frame (must be at least SYNC_WORD_LEN + RS_BLOCK_LEN * n_interleave)
     * @return             Total number of bytes written to the output buffer
     */
    int encode(const uint8_t* in_payload, uint8_t* out_frame);

    /**
     * @return Number of frames transmitted
     */
    uint32_t num_frames() const { return d_num_frames; }

private:
    inline int data_len() const { return RS_DATA_LEN * d_n_interleave; }
    inline int codeword_len() const { return RS_BLOCK_LEN * d_n_interleave; }
    inline int total_frame_len() const { return SYNC_WORD_LEN + codeword_len(); }

    bool d_rs_encode;
    bool d_interleave;
    bool d_scramble;
    bool d_printing;
    bool d_verbose;
    int d_n_interleave;
    bool d_dual_basis;

    reed_solomon d_rs;
    struct ccsds_tx_pkt d_pkt;


    uint32_t d_num_frames = 0;
};

#endif // ccsds_rs_encoder_H