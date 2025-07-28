// Standalone CCSDS Correlator (GNU Radio dependencies removed)

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


class ccsds_correlator
{
public:
    enum state_t { SEARCH, LOCK };
    enum ambiguity_t { NONE, INVERTED };

    ccsds_correlator(uint64_t asm_word, uint64_t asm_mask, uint8_t threshold, size_t frame_len);
    ~ccsds_correlator();

    int process(const uint8_t* in, int ninput_items, uint8_t* out_frame, bool* frame_ready);
    uint64_t frame_count() const { return d_frame_count; }

private:
    bool check_asm(uint64_t asm_buf);
    void enter_state(state_t state);

    uint64_t d_asm;
    uint64_t d_asm_mask;
    uint8_t d_threshold;
    size_t d_frame_len;

    uint64_t d_asm_buf = 0;
    uint8_t d_byte_buf = 0;
    uint8_t d_bit_ctr = 0;
    uint8_t* d_frame_buffer = nullptr;
    size_t d_frame_buffer_len = 0;
    uint64_t d_frame_count = 0;

    state_t d_state = SEARCH;
    ambiguity_t d_ambiguity = NONE;
};

ccsds_correlator::ccsds_correlator(uint64_t asm_word, uint64_t asm_mask, uint8_t threshold, size_t frame_len)
    : d_asm(asm_word), d_asm_mask(asm_mask), d_threshold(threshold), d_frame_len(frame_len)
{
    d_frame_buffer = (uint8_t*) malloc(frame_len * sizeof(uint8_t));
    enter_state(SEARCH);
}

ccsds_correlator::~ccsds_correlator()
{
    free(d_frame_buffer);
}

int ccsds_correlator::process(const uint8_t* in, int ninput_items, uint8_t* out_frame, bool* frame_ready)
{
    *frame_ready = false;
    uint64_t count = 0;
    while (count < ninput_items)
    {
        switch (d_state)
        {
            case SEARCH:
            {
                d_asm_buf = (d_asm_buf << 1) | (in[count++] & 0x01);
                if (check_asm(d_asm_buf))
                {
                    d_ambiguity = NONE;
                    enter_state(LOCK);
                } else if (check_asm(d_asm_buf ^ 0xffffffffffffffffULL))
                {
                    d_ambiguity = INVERTED;
                    enter_state(LOCK);
                }
                break;
            }
            case LOCK:
            {
                d_byte_buf = (d_byte_buf << 1) | (in[count++] & 0x01);
                d_bit_ctr++;
                if (d_bit_ctr == 8)
                {
                    if (d_ambiguity == NONE)
                    {
                        d_frame_buffer[d_frame_buffer_len] = d_byte_buf;
                    } 
                    else
                    {
                        d_frame_buffer[d_frame_buffer_len] = d_byte_buf ^ 0xff;
                    }
                    d_frame_buffer_len++;
                    d_bit_ctr = 0;
                }
                if (d_frame_buffer_len == d_frame_len)
                {
                    memcpy(out_frame, d_frame_buffer, d_frame_len);
                    *frame_ready = true;
                    d_frame_count++;
                    enter_state(SEARCH);
                    return count; // return after one complete frame
                }
                break;
            }
        }
    }
    return count;
}

bool ccsds_correlator::check_asm(uint64_t asm_buf)
{
    uint64_t syndrome = (asm_buf ^ d_asm) & d_asm_mask;
    uint64_t nerrors = 0;
    //volk_64u_popcnt(&nerrors, syndrome);
    nerrors = __builtin_popcountll(syndrome);
    return nerrors <= d_threshold;
}

void ccsds_correlator::enter_state(state_t state)
{
    switch (state)
    {
        case SEARCH:
        {
            d_asm_buf = 0;
            break;
        }
        case LOCK:
        {
            d_byte_buf = 0;
            d_bit_ctr = 0;
            d_frame_buffer_len = 0;
            break;
        }
    }
    d_state = state;
}