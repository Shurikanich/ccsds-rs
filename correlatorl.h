#ifndef CCSDS_CORRELATOR_H
#define CCSDS_CORRELATOR_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// Standalone CCSDS Correlator (no GNU Radio)

class ccsds_correlator {
public:
    enum state_t { SEARCH, LOCK };
    enum ambiguity_t { NONE, INVERTED };

    ccsds_correlator(uint64_t asm_word, uint64_t asm_mask, uint8_t threshold, size_t frame_len);
    ~ccsds_correlator();

    /**
     * Process input bits to detect ASM and extract CCSDS frames.
     * 
     * @param in            Pointer to input bitstream (1-bit packed as uint8_t)
     * @param ninput_items  Number of input bits
     * @param out_frame     Pointer to buffer where decoded frame will be written
     * @param frame_ready   Set to true if a complete frame is detected and copied to out_frame
     * @return              Number of bits consumed
     */
    int process(const uint8_t* in, int ninput_items, uint8_t* out_frame, bool* frame_ready);

    /**
     * Get total number of successfully extracted frames.
     */
    uint64_t frame_count() const;

private:
    bool check_asm(uint64_t asm_buf);
    void enter_state(state_t state);

    // Config parameters
    uint64_t d_asm;
    uint64_t d_asm_mask;
    uint8_t  d_threshold;
    size_t   d_frame_len;

    // State variables
    uint64_t d_asm_buf;
    uint8_t  d_byte_buf;
    uint8_t  d_bit_ctr;
    size_t   d_frame_buffer_len;
    uint8_t* d_frame_buffer;
    uint64_t d_frame_count;
    state_t d_state;
    ambiguity_t d_ambiguity;
};

#endif // CCSDS_CORRELATOR_H