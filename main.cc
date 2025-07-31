// Modified RS-only simulation to include convolutional encoding and decoding
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "ccsds_rs_encoder.h"
#include "ccsds_rs_decoder.h"
#include "ccsds.h"
#include "viterbi27.h"

using namespace std;

// Function to generate AWGN noise
static double gaussian_noise(double stddev)
{
    static bool has_saved = false;
    static double saved;

    if (has_saved)
    {
        has_saved = false;
        return stddev * saved;
    }

    double u1, u2;
    do {
        u1 = rand() / (double)RAND_MAX;
        u2 = rand() / (double)RAND_MAX;
    } while (u1 <= 1e-10);

    double mag = sqrt(-2.0 * log(u1));
    double z0 = mag * cos(2.0 * M_PI * u2);
    double z1 = mag * sin(2.0 * M_PI * u2);

    saved = z1;
    has_saved = true;

    return stddev * z0;
}

// Soft-decision mapping for BPSK with optional puncturing
static unsigned char soft_decision(double x, bool is_punctured = false)
{
    if (is_punctured) return 128;
    double normalized = (x + 1.0) / 2.0;
    int value = int(normalized * 255.0 + 0.5);
    value = std::min(255, std::max(0, value));
    return 0xFF - value;
}

// Configurable parameters
constexpr bool rs_encode = true;
constexpr bool interleave = true;
constexpr bool scramble_val = true;
constexpr bool printing = false;
constexpr bool verbose = false;
constexpr int n_interleave = 8;
constexpr bool dual_basis = true;

int main(int argc, char* argv[])
{
    // Read convolutional code config file
    string config_filename = argc > 1 ? argv[1] : "config.txt";
    double snr_db = 100; // Default SNR in dB
    string puncturing_type = "1/2";
    ifstream config(config_filename);
    if (config)
    {
        string line;
        while (getline(config, line))
        {
            if (line.empty() || line[0] == '#') continue;
            istringstream iss(line);
            string key, value;
            if (getline(iss, key, '=') && getline(iss, value))
            {
                if (key == "start_snr") snr_db = stod(value);
                else if (key == "puncturing_type") puncturing_type = value;
            }
        }
    }

    // Select puncturing
    double code_rate;
    const int* puncture_C1_ptr = nullptr;
    const int* puncture_C2_ptr = nullptr;
    int puncture_pattern_len = 0;
    if (puncturing_type == "1/2")
    {
        code_rate = CODE_RATE_12;
        puncture_pattern_len = PUNCTURE_PATTERN_LEN_12;
        puncture_C1_ptr = puncture_C1_12;
        puncture_C2_ptr = puncture_C2_12;
    }
    else
    {
        cerr << "Unsupported puncturing type." << endl;
        return 1;
    }

    srand(time(nullptr));
    // RS Encode
    ccsds_rs_encoder encoder(rs_encode, interleave, scramble_val, printing, verbose, n_interleave, dual_basis);
    // RS Decode
    ccsds_rs_decoder decoder(0, rs_encode, interleave, scramble_val, verbose, printing, n_interleave, dual_basis);

    const int payload_len = RS_DATA_LEN * n_interleave;
    const int frame_len = SYNC_WORD_LEN + RS_BLOCK_LEN * n_interleave;
    uint8_t input_payload[payload_len];
    uint8_t encoded_frame[frame_len];
    uint8_t decoded_output[payload_len];

    // Convolutional decode initialization
    v27 vi;
    memset(&vi, 0, sizeof(v27));
    vitfilt27_init(&vi);


    // BPSK modulation + AWGN
    double EbN0 = pow(10.0, snr_db / 10.0);
    // TODO: adjust code rate for reed solomon
    double N0 = 1.0 / (2.0 * EbN0 * code_rate);
    double noise_std = sqrt(N0);

    // Generate input data
    for (int i = 0; i < payload_len; ++i)
    {
        input_payload[i] = rand() % 256;
    }

    
    int encoded_len = encoder.encode(input_payload, encoded_frame);

    if (verbose)
    {
        std::cout << "\n--- encoded_frame ---\n";
        print_bytes(encoded_frame, encoded_len);
    }

    // Convolutional encode
    // TODO: here check the size of conv_encoded
    // encode produces 2 bits for every input bit. therefore, 8 bits pro byte * 2 bits = 16 times the size of the input
    unsigned char conv_encoded[frame_len * 16] = {0};
    unsigned char state = 0;
    unsigned int conv_len = encode27(&state, conv_encoded, encoded_frame, encoded_len,
                                     puncture_C1_ptr, puncture_C2_ptr, puncture_pattern_len);
    if (verbose)
    {
        std::cout << "\n--- conv_encoded ---\n";
        print_bytes(conv_encoded, conv_len);
    }



    double received[conv_len];
    unsigned char soft[conv_len];
    for (unsigned i = 0; i < conv_len; ++i)
    {
        double bpsk = (conv_encoded[i] == 0) ? 1.0 : -1.0;
        received[i] = bpsk + gaussian_noise(noise_std);
    }

    // Convert to soft decisions

    int idx = 0, c1 = 0, c2 = 0;
    for (int i = 0; i < (int)conv_len; ++i)
    {
        bool is_c1 = (i % 2 == 0);
        bool is_punct = is_c1 ? !puncture_C1_ptr[c1++ % puncture_pattern_len] : !puncture_C2_ptr[c2++ % puncture_pattern_len];
        soft[i] = soft_decision(received[i], is_punct);
    }


    unsigned char conv_decoded[frame_len + 16] = {0}; // +16 for traceback bias
    vitfilt27_decode(&vi, soft, conv_decoded, conv_len + 256);

    if (verbose)
    {
        std::cout << "\n--- conv_decoded ---\n";
        print_bytes(conv_decoded, frame_len + 16);

        // Compare original encoded_frame and conv_decoded
        std::cout << "\n--- Comparing conv_decoded with encoded_frame ---\n";

        int bit_errors = 0;
        int bit_count = 0;
        for (int i = 0; i < frame_len; ++i)
        {
            uint8_t original = encoded_frame[i];
            uint8_t decoded  = conv_decoded[i + 16];  // skip traceback bias
            uint8_t diff = original ^ decoded;

          if (diff != 0)
          {
            std::cout << "Byte " << std::setw(4) << i << ": "
                      << "original = 0x" << std::hex << std::setw(2) << std::setfill('0') << (int)original
                      << ", decoded = 0x" << std::hex << std::setw(2) << std::setfill('0') << (int)decoded
                      << ", diff bits = " << std::dec << std::bitset<8>(diff) << "\n";
          }

          for (int b = 0; b < 8; ++b)
          {
            if (diff & (1 << b)) ++bit_errors;
            ++bit_count;
          }
        }
  
      std::cout << "\nBit error count: " << bit_errors << " / " << bit_count
                << "  => BER = " << (double)bit_errors / bit_count << "\n";
    }

    int noutput_items = 0;
    
#if 0
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
#endif

    decoder.decode_aligned_bytes(&conv_decoded[16], frame_len, decoded_output, &noutput_items);

    // Validate
    bool match = (noutput_items == payload_len && memcmp(input_payload, decoded_output, payload_len) == 0);
    if (match)
    std::cout << "\n\033[1;32m✅ Decode SUCCESS\033[0m\n";  // Green bold
    else
    std::cout << "\n\033[1;31m❌ Decode FAILED\033[0m\n";   // Red bold
    return match ? 0 : 1;
}
