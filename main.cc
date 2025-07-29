#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include "ccsds_rs_encoder.h"
#include "ccsds_rs_decoder.h"
#include "ccsds.h"
#include <iomanip>

// Configuration parameters
constexpr bool rs_encode = true;
constexpr bool interleave = true;
constexpr bool scramble_val = true;
constexpr bool printing = false;
constexpr bool verbose = false;
constexpr int n_interleave = 8;
constexpr bool dual_basis = true;

int main() {
    srand(time(nullptr));

    const int payload_len = RS_DATA_LEN * n_interleave;
    const int frame_len = SYNC_WORD_LEN + RS_BLOCK_LEN * n_interleave;

    uint8_t input_payload[payload_len];
    uint8_t encoded_frame[frame_len];
    uint8_t decoded_output[payload_len];

    // Fill input with pseudo-random data
    //int val = 1;
    for (int i = 0; i < payload_len; ++i)
    {
        input_payload[i] = rand() % 256;
        //input_payload[i] = val; //i % RS_DATA_LEN; //rand() % 256;
        //if ((i % RS_DATA_LEN == 0) && (i != 0)) val++; // Increment every RS_DATA_LEN bytes
    }
    if (verbose)
    {
      std::cout << "Original input payload:\n";
      print_bytes(input_payload, payload_len);
    }

    // Encode
    ccsds_rs_encoder encoder(rs_encode, interleave, scramble_val, printing, verbose, n_interleave, dual_basis);
    int encoded_len = encoder.encode(input_payload, encoded_frame);
    std::cout << "\nEncoded frame length: " << encoded_len << "\n";

#if 1
    if (verbose)
    {
      std::cout << "\n--- Encoded Frame Breakdown ---\n";
      constexpr bool interleave_output = true;

      // Print SYNC_WORD
      std::cout << "[SYNC_WORD] (" << SYNC_WORD_LEN << " bytes):\n";
      print_bytes(encoded_frame, SYNC_WORD_LEN);

      // Print each RS block
      for (int i = 0; i < n_interleave; ++i)
      {
        std::cout << "[RS_BLOCK #" << i << "] (" << RS_BLOCK_LEN << " bytes):\n";

        if (interleave_output)
        {
          // Deinterleave-style layout (transpose)
          for (int j = 0; j < RS_BLOCK_LEN; ++j)
          {
              std::cout << std::hex << std::uppercase
                        << "0x" << std::setw(2) << std::setfill('0')
                        << static_cast<int>(encoded_frame[SYNC_WORD_LEN + i + j * n_interleave]) << " ";
              if ((j + 1) % 16 == 0) std::cout << "\n";
          }
        } 
        else
        {
          // Block-by-block layout
          print_bytes(&encoded_frame[SYNC_WORD_LEN + i * RS_BLOCK_LEN], RS_BLOCK_LEN);
        }

        std::cout << "\n";
      }
    }
#endif

    // Convert to bitstream for decoder
    uint8_t bitstream[frame_len * 8];
    for (int i = 0; i < encoded_len; ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            bitstream[i * 8 + j] = (encoded_frame[i] >> (7 - j)) & 1;
        }
    }

    // Decode
    ccsds_rs_decoder decoder(0, rs_encode, interleave, scramble_val, verbose, printing, n_interleave, dual_basis);
    int noutput_items = 0;
    decoder.process(bitstream, encoded_len * 8, decoded_output, &noutput_items);

    if (verbose)
    {
        std::cout << "\n--- Decoded Frame Breakdown ---\n";
        for (int i = 0; i < n_interleave; ++i)
        {
            std::cout << "[RS_BLOCK #" << i << "] (" << RS_BLOCK_LEN << " bytes):\n";
            print_bytes(&decoded_output[i * RS_BLOCK_LEN], RS_BLOCK_LEN);
        }
        std::cout << "\nDecoded output payload (" << noutput_items << " bytes):\n";
        print_bytes(decoded_output, noutput_items);
    }


    // Validate
    bool match = (noutput_items == payload_len) && (memcmp(input_payload, decoded_output, payload_len) == 0);
    std::cout << "\n✅ Decode " << (match ? "SUCCESS" : "FAILED") << "\n";

    return match ? 0 : 1;
}