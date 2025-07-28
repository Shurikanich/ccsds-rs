

#include "reed_solomon.h"

#include <stdint.h>
#include <assert.h>
#include <stdio.h>

extern "C" {
#include "fec-3.0.1/fec.h"
}
#include "ccsds.h"

extern unsigned char CCSDS_alpha_to[];
extern unsigned char CCSDS_index_of[];
extern unsigned char CCSDS_poly[];
extern unsigned char Taltab[];
extern unsigned char Tal1tab[];



reed_solomon::reed_solomon() {}
reed_solomon::~reed_solomon() {}

void reed_solomon::encode(uint8_t *data, bool use_dual_basis)
{
    if (use_dual_basis)
    {
        encode_rs_ccsds(data, &data[RS_DATA_LEN],0);
    }
    else
    {
        encode_rs_8(data, &data[RS_DATA_LEN], 0);
    }
}

int16_t reed_solomon::decode(uint8_t *data, bool use_dual_basis)
{
    if (use_dual_basis)
    {
        return decode_rs_ccsds(data, 0, 0, 0);
    }
    else
    {
        return decode_rs_8(data, 0, 0, 0);
    }

}
