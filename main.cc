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

#include <cctype>

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

// Compute bit errors for one frame and update running totals.
// ref         : original payload/frame to compare against
// dec         : decoded bytes buffer
// len         : number of bytes to compare
inline void calc_errors(const uint8_t* ref,
                        const uint8_t* dec,
                        int            len,
                        unsigned long& total_errors,
                        unsigned long& total_bits,
                        int            dec_offset = 0)
{
    if (!ref || !dec || len <= 0) return;

    unsigned long bit_errors = 0;

    for (int j = 0; j < len; ++j)
    {
        uint8_t diff = static_cast<uint8_t>(ref[j] ^ dec[j + dec_offset]);
#if defined(__GNUC__) || defined(__clang__)
        bit_errors += static_cast<unsigned long>(__builtin_popcount(static_cast<unsigned>(diff)));
#else
        // Portable fallback
        for (int b = 0; b < 8; ++b)
            bit_errors += (diff >> b) & 0x1u;
#endif
    }

    total_errors += bit_errors;
    total_bits   += static_cast<unsigned long>(len) * 8ul;
}


static std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    size_t e = s.find_last_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    return s.substr(b, e - b + 1);
}
static std::string lower(std::string s) {
    for (auto& c : s) c = std::tolower(static_cast<unsigned char>(c));
    return s;
}
static bool parse_bool(const std::string& v, bool defval=false) {
    std::string s = lower(trim(v));
    if (s=="1" || s=="true" || s=="yes" || s=="on")  return true;
    if (s=="0" || s=="false"|| s=="no"  || s=="off") return false;
    return defval;
}



int main(int argc, char* argv[])
{
// --------- Read configuration from file ------------
    double start_snr, end_snr, step_snr;
    string output_filename;
    string puncturing_type;
    int packet_size;
    int i, j;
    int num_bits = 1000000;  // total number of bits to simulate per Eb/N0 point

    // -------- Runtime-configurable parameters (with defaults) --------
    typedef enum { ONLY_RS, ONLY_CC, RS_AND_CC } ccsds_mode_t;

    bool rs_encode     = true;
    bool interleave    = true;
    bool scramble_val  = true;
    bool printing      = false;
    bool verbose       = false;
    int  n_interleave  = 8;
    bool dual_basis    = true;
    ccsds_mode_t mode  = RS_AND_CC;

    // Use command-line argument for config file name if provided.
    string config_filename = "config.txt";
    if (argc > 1)
    {
            config_filename = argv[1];
    }
    cout << "Using config file: " << config_filename << endl;
    ifstream config(config_filename);
    if (!config)
    {
        cerr << "Error: Could not open config file:" << config_filename << endl;
        exit(1);
    }
    string line;
    while(getline(config, line))
    {
        // Skip empty lines or comments
        if(line.empty() || line[0]=='#')
            continue;
        istringstream iss(line);
        string key, value;
        if(getline(iss, key, '=') && getline(iss, value))
        {
          key   = trim(key);
          value = trim(value);

          if      (key == "start_snr")       start_snr = stod(value);
          else if (key == "end_snr")         end_snr   = stod(value);
          else if (key == "step_snr")        step_snr  = stod(value);
          else if (key == "output_file")     output_filename = value;
          else if (key == "puncturing_type") puncturing_type = value;



          else if (key == "rs_encode")       rs_encode    = parse_bool(value, rs_encode);
          else if (key == "interleave")      interleave   = parse_bool(value, interleave);
          else if (key == "scramble")        scramble_val = parse_bool(value, scramble_val);
          else if (key == "printing")        printing     = parse_bool(value, printing);
          else if (key == "verbose")         verbose      = parse_bool(value, verbose);
          else if (key == "n_interleave")    n_interleave = std::max(1, stoi(value));
          else if (key == "dual_basis")      dual_basis   = parse_bool(value, dual_basis);
          else if (key == "mode")
          {
            std::string m = lower(value);
            if      (m == "only_rs")   mode = ONLY_RS;
            else if (m == "only_cc")   mode = ONLY_CC;
            else if (m == "rs_and_cc") mode = RS_AND_CC;
            else
            {
                std::cerr << "[WARN] Unknown mode '" << value
                          << "'. Using RS_AND_CC.\n";
                mode = RS_AND_CC;
            }
          }
        }
    }
    config.close();



    // Replace '/' with '_' in puncturing_type.
    string punc = puncturing_type;
    size_t pos = 0;
    while ((pos = punc.find('/', pos)) != string::npos)
    {
        punc.replace(pos, 1, "_");
        pos++; // Move past the replaced character
    }

    ostringstream oss;
    oss << "res/" << output_filename
    << "_start" << start_snr
    << "_end" << end_snr
    << "_step" << step_snr
    << "_" << punc
    << "_mode";

    switch (mode)
    {
      case ONLY_RS:   oss << "ONLY_RS"; break;
      case ONLY_CC:   oss << "ONLY_CC"; break;
      case RS_AND_CC: oss << "RS_AND_CC"; break;
    }

    oss << "_intlv" << n_interleave << (dual_basis ? "_dualBasis" : "_noDualBasis") << ".txt";

    output_filename = oss.str();


    // Open the output file to write BER results
    ofstream results(output_filename);
    if(!results)
    {
        cerr << "Error: Could not open output file " << output_filename << endl;
        exit(1);
    }
    cout << "Using output file: " << output_filename << endl;


    // ------------------------------
    // Set up puncturing parameters based on configuration
    // ------------------------------
    double code_rate_cc;
    const int* puncture_C1_ptr = nullptr;
    const int* puncture_C2_ptr = nullptr;
    int puncture_pattern_len = 0;

    if(puncturing_type == "1/2")
    {
        code_rate_cc = CODE_RATE_12;
        puncture_pattern_len = PUNCTURE_PATTERN_LEN_12;
        puncture_C1_ptr = puncture_C1_12;
        puncture_C2_ptr = puncture_C2_12;
    }
    else if(puncturing_type == "3/4")
    {
        code_rate_cc = CODE_RATE_34;
        puncture_pattern_len = PUNCTURE_PATTERN_LEN_34;
        puncture_C1_ptr = puncture_C1_34;
        puncture_C2_ptr = puncture_C2_34;
    }
    else if(puncturing_type == "7/8")
    {
        code_rate_cc = CODE_RATE_78;
        puncture_pattern_len = PUNCTURE_PATTERN_LEN_78;
        puncture_C1_ptr = puncture_C1_78;
        puncture_C2_ptr = puncture_C2_78;
    }
    else if(puncturing_type == "2/3")
    {
        code_rate_cc = CODE_RATE_23;
        puncture_pattern_len = PUNCTURE_PATTERN_LEN_23;
        puncture_C1_ptr = puncture_C1_23;
        puncture_C2_ptr = puncture_C2_23;
    }
    else if(puncturing_type == "5/6")
    {
        code_rate_cc = CODE_RATE_56;
        puncture_pattern_len = PUNCTURE_PATTERN_LEN_56;
        puncture_C1_ptr = puncture_C1_56;
        puncture_C2_ptr = puncture_C2_56;
    }
    else
    {
        cerr << "Unknown puncturing type: " << puncturing_type << endl;
        exit(1);
    }

    if (mode == RS_AND_CC || mode == ONLY_CC)
    {
      std::cout << std::fixed << std::setprecision(3);
      std::cout << "code_rate_cc = " << code_rate_cc << std::endl;

      std::cout << "Puncture C1: [ ";
      for (int i = 0; i < puncture_pattern_len; ++i)
          std::cout << puncture_C1_ptr[i] << " ";
      std::cout << "]" << std::endl;

      std::cout << "Puncture C2: [ ";
      for (int i = 0; i < puncture_pattern_len; ++i)
          std::cout << puncture_C2_ptr[i] << " ";
      std::cout << "]" << std::endl;
    }

    
    switch (mode)
    {
      case ONLY_RS:
        std::cout << "\033[1;34m[Mode: ONLY_RS]\033[0m\n";  // Blue bold
        break;
      case ONLY_CC:
        std::cout << "\033[1;36m[Mode: ONLY_CC]\033[0m\n";  // Cyan bold
        break;
      case RS_AND_CC:
        std::cout << "\033[1;33m[Mode: RS_AND_CC]\033[0m\n"; // Yellow bold
        break;
    }

    srand(time(nullptr));
    // RS Encode
    ccsds_rs_encoder encoder(rs_encode, interleave, scramble_val, printing, verbose, n_interleave, dual_basis);
    // RS Decode
    ccsds_rs_decoder decoder(0, rs_encode, interleave, scramble_val, verbose, printing, n_interleave, dual_basis);

    int payload_len = RS_DATA_LEN * n_interleave;

    const int frame_len = SYNC_WORD_LEN + RS_BLOCK_LEN * n_interleave;
    if (mode == ONLY_CC)
    {
       payload_len = frame_len; // For ONLY_CC, payload is the same as frame length
    }
    uint8_t input_payload[payload_len];
    // add one since convolutional encoder not decodes the last byte make it always 0
    // 3 zeros at the end  and 5 at the beginning
    uint8_t encoded_frame[frame_len + 8];
    uint8_t decoded_output[payload_len];

    //cout << "Payload length: " << payload_len << " bytes" << endl;
    //cout << "Frame length: " << frame_len << " bytes" << endl;
    //cout << "Encoded frame length: " << frame_len + 1 << " bytes (+1 for CC padding)" << endl;

    // add one since convolutional encoder not decodes the last 3 bytes make it always 0
    // 3 zeros at the end and 5 at the beginning
    unsigned int conv_len = (frame_len + 8) * 16; // 8 bits per byte * 2 bits per input bit for convolutional encoding

    //cout << "Convolutional encoded length (conv_len): " << conv_len << " bits" << endl;

    // Convolutional decode initialization
    v27 vi;
    memset(&vi, 0, sizeof(v27));
    vitfilt27_init(&vi);

    // Generate SNR values from config (in dB)
    vector<double> EbN0_values;
    for (double snr = start_snr; snr <= end_snr; snr += step_snr)
    {
        EbN0_values.push_back(snr);
    }
    int num_ebn0_points = EbN0_values.size();

    for (i = 0; i < num_ebn0_points; i++)
    {
        double EbN0 = pow(10.0, EbN0_values[i] / 10.0);
        unsigned long total_errors = 0;
        unsigned long total_bits = 0;

        unsigned long nbits = (num_bits * pow(2.0, EbN0_values[i] / 2.0));
        //unsigned long nbits = (num_bits * 100);

        unsigned last_pct = 0; // for the progress print

      // Simulate until we process at least num_bits bits
      while (total_bits < nbits)
      {
        // BPSK modulation + AWGN
        double EbN0 = pow(10.0,  EbN0_values[i] / 10.0);
        double N0 = 0.0;
        // code rate for RS
        double code_rate_rs = static_cast<double>(RS_DATA_LEN) / static_cast<double>(RS_BLOCK_LEN);
        switch (mode)
        {
          case ONLY_RS:
            N0 = 1.0 / (2.0 * EbN0 * code_rate_rs);
            break;
          case ONLY_CC:
            N0 = 1.0 / (2.0 * EbN0 * code_rate_cc);
            break;
          case RS_AND_CC:
            N0 = 1.0 / (2.0 * EbN0 * code_rate_rs * code_rate_cc);
            break;
        }
        
        double noise_std = sqrt(N0);

        // variables for encoding and decoding
        int encoded_len = 0; 
        int noutput_items = 0;

        if (mode == RS_AND_CC || mode == ONLY_RS)
        {
            // Generate input data
            for (int i = 0; i < payload_len; ++i)
            {
                input_payload[i] = rand() % 256;
            }
            // make 3 byte  always 0 and 5 bytes at the beginning
            encoded_len = encoder.encode(input_payload, encoded_frame) + 8;
            // here is the last byte always 0
            encoded_frame[encoded_len - 1] = 0;
            encoded_frame[encoded_len - 2] = 0;
            encoded_frame[encoded_len - 3] = 0;




            if (verbose)
            {
              std::cout << "\n--- encoded_frame ---\n";
              print_bytes(encoded_frame, encoded_len);
            }
        }
        else
        {
            // mode == ONLY_CC
            // Generate random data for convolutional encoding
            for (int i = 0; i < frame_len; ++i)
            {
                encoded_frame[i] = rand() % 256;
            }
            encoded_len = frame_len + 8;
            // make 3 byte  always 0 and 5 bytes at the beginning
            encoded_frame[encoded_len - 1] = 0;
            encoded_frame[encoded_len - 2] = 0;
            encoded_frame[encoded_len - 3] = 0;

            // Add 5 zeros at the beginning
            encoded_frame[0] = 0;
            encoded_frame[1] = 0;
            encoded_frame[2] = 0;
            encoded_frame[3] = 0;
            encoded_frame[4] = 0;

            if (0)
            {
              std::cout << "\n--- (Initial) encoded_frame in the case of ONLY_CC---\n";
              print_bytes(encoded_frame, frame_len);
            }

        }

        

        if (mode == RS_AND_CC || mode == ONLY_CC)
        {
          // Convolutional encode
          // TODO: here check the size of conv_encoded
          // encode produces 2 bits for every input bit. therefore, 8 bits pro byte * 2 bits = 16 times the size of the input
          // 5 bytes at the beginning and 3 (totally 8) at the end are always 0
          unsigned char conv_encoded[(frame_len + 8) * 16];
          unsigned char state = 0;
          unsigned int conv_len_real = encode27(&state, conv_encoded, encoded_frame, encoded_len,
                                          puncture_C1_ptr, puncture_C2_ptr, puncture_pattern_len);

          //cout << "Convolutional encoded length (conv_len_real): " << conv_len_real << " bits" << endl;

          if (verbose)
          {
              std::cout << "\n--- conv_encoded ---\n";
              print_bytes(conv_encoded, conv_len_real);
          }

          



          double received[conv_len];
          unsigned char soft[conv_len];
          for (unsigned i = 0; i < conv_len_real; ++i)
          {
              double bpsk = (conv_encoded[i] == 0) ? 1.0 : -1.0;
              received[i] = bpsk + gaussian_noise(noise_std);
          }


#if 0
          // Convert to soft decisions

          int idx = 0, c1 = 0, c2 = 0;
          for (int i = 0; i < (int)conv_len; ++i)
          {
              bool is_c1 = (i % 2 == 0);
              bool is_punct = is_c1 ? !puncture_C1_ptr[c1++ % puncture_pattern_len] : !puncture_C2_ptr[c2++ % puncture_pattern_len];
              soft[i] = soft_decision(received[i], is_punct);
          }

#else
            int idx = 0;         // Index into received_signal array.
            int c1_count = 0;    // Counter for C1 bits.
            int c2_count = 0;    // Counter for C2 bits.
            
            for (int j = 0; j < (int)conv_len; j++)
            {
                if (j % 2 == 0)
                {
                    // Even index: C1 bit.
                    bool is_punctured = !puncture_C1_ptr[c1_count % puncture_pattern_len];
                    soft[j] = soft_decision(received[idx], is_punctured);
                    if (!is_punctured)
                    {
                        idx++;
                    }
                    c1_count++;
                } 
                else
                {
                    // Odd index: C2 bit.
                    bool is_punctured = !puncture_C2_ptr[c2_count % puncture_pattern_len];
                    soft[j] = soft_decision(received[idx], is_punctured);
                    if (!is_punctured)
                    {
                        idx++;
                    }
                    c2_count++;
                }
            }
#endif




          unsigned char conv_decoded[frame_len + 16]; // +16 for traceback bias

          //
          //cout << "5.5" << endl;
          //cout << "conv_len = " << conv_len << endl;
          //cout << "frame_len = " << frame_len << endl;
          //cout << "conv_decoded size = " << sizeof(conv_decoded) << endl;
          //

          vitfilt27_decode(&vi, soft, conv_decoded, conv_len + 256);

          // first 5 bytes at thhe beginning are set always to 0
          conv_decoded[0 + 16] = 0;
          conv_decoded[1 + 16] = 0;
          conv_decoded[2 + 16] = 0;
          conv_decoded[3 + 16] = 0;
          conv_decoded[4 + 16] = 0;

          if (0)
          {
              //std::cout << "\n--- conv_decoded ---\n";
              //print_bytes(conv_decoded, frame_len + 16);

              // Compare original encoded_frame and conv_decoded
              // std::cout << "\n--- Comparing conv_decoded with encoded_frame ---\n";

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
        
            // std::cout << "\nBit error count: " << bit_errors << " / " << bit_count << "  => BER = " << (double)bit_errors / bit_count << "\n";
          }

          if (mode == RS_AND_CC)
          {
            decoder.decode_aligned_bytes(&conv_decoded[16], frame_len, decoded_output, &noutput_items);
          }
          else
          {
            // mode == ONLY_CC
            for (int i = 0; i < frame_len; ++i)
            {
                decoded_output[i] = conv_decoded[16 + i]; // skip traceback bias
            }
            noutput_items = frame_len;

            if (0)
            {
              std::cout << "\n--- decoded_output in the case of ONLY_CC---\n";
              print_bytes(decoded_output, frame_len);
            }
          }
        }
        else
        {
          // RS Decode only mode == ONLY_RS
          // Convert to bitstream for decoder
          uint8_t bitstream[frame_len * 8];
          for (int i = 0; i < encoded_len; ++i)
          {
            for (int j = 0; j < 8; ++j)
            {
                bitstream[i * 8 + j] = (encoded_frame[i] >> (7 - j)) & 1;
            }
          }

          // BPSK modulation
          vector<double> modulated(encoded_len);
          for (size_t i = 0; i < encoded_len; ++i)
          {
            modulated[i] = bitstream[i] ? -1.0 : 1.0;
          }

          // Add AWGN noise
          vector<double> received(modulated.size());
          for (size_t i = 0; i < modulated.size(); ++i)
          {
            received[i] = modulated[i] + gaussian_noise(noise_std);
          }

          // Hard-decision demodulation
          for (size_t i = 0; i < received.size(); ++i)
          {
              bitstream[i] = received[i] >= 0.0 ? 0 : 1;
          }

          decoder.find_asm_and_decode(bitstream, encoded_len * 8, decoded_output, &noutput_items);
        }
        
        // Validate
        bool match = false;
        if (mode == ONLY_CC)
        {

            if (0)
            {
              std::cout << "\n (Second)--- encoded_frame in the case of ONLY_CC---\n";
              print_bytes(encoded_frame, frame_len);
            }

            //if (1)
            //{
            //  std::cout << "\n--- decoded_output in the case of ONLY_CC---\n";
            //  print_bytes(decoded_output, frame_len);
            //}

#if 0
            match = memcmp(encoded_frame, decoded_output, frame_len) == 0;
            if (!match)
            {
              for (int i = 0; i < frame_len; ++i)
              {
                if (encoded_frame[i] != decoded_output[i])
                {
                  match = false;
                  std::cout << "Mismatch at byte " << i
                            << ": encoded_frame = 0x" << std::hex << std::setw(2) << std::setfill('0') << (int)encoded_frame[i]
                            << ", decoded_output = 0x" << std::hex << std::setw(2) << std::setfill('0') << (int)decoded_output[i]
                            << std::dec << "\n";
                }
              }
            }

#endif
            //match = (noutput_items == frame_len && memcmp(encoded_frame, decoded_output, frame_len) == 0);
            calc_errors(encoded_frame, decoded_output, frame_len, total_errors, total_bits, 0);
        }
        else
        {
            //match = (noutput_items == payload_len && memcmp(input_payload, decoded_output, payload_len) == 0);
            calc_errors(input_payload, decoded_output, payload_len, total_errors, total_bits, 0);
        }
        //bool match = (noutput_items == payload_len && memcmp(input_payload, decoded_output, payload_len) == 0);
        //if (match)
        //std::cout << "\n\033[1;32m✅ Decode SUCCESS\033[0m\n";  // Green bold
        //else
        //std::cout << "\n\033[1;31m❌ Decode FAILED\033[0m\n";   // Red bold
        //return match ? 0 : 1;

        // --- Progress print ---
        unsigned pct = static_cast<unsigned>(100.0 * total_bits / nbits);
        if (pct != last_pct && pct % 2 == 0)
        { // print every 2%
          cout << "\r[" << setw(3) << pct << "%] Done..." << flush;
          last_pct = pct;
        }
      }
      double ber = (double)total_errors / total_bits;
      cout << fixed << setprecision(2) << "Eb/N0 (dB) = " << EbN0_values[i] << ", BER = " << scientific << setprecision(2) << ber << endl;
      //results << "Eb/N0 (dB) = " << EbN0_values[i] << ", BER = " << ber << endl;
      results <<  EbN0_values[i] << " " << ber << endl;
  }
      
  results.close();
  return 0;
}
