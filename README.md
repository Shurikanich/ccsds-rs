## 🧪 ccsds-rs

Simulation chain for the CCSDS Reed-Solomon (RS) encoder and decoder, including a convolutional code according to the CCSDS TM standard.
The convolutional code can be configured for different puncturing rates: 1/2, 2/3, 3/4, 5/6, 7/8.
The simulation chain supports three different modes:
	•	only_cc – only convolutional code
	•	only_rs – only Reed-Solomon code
	•	rs_and_cc – concatenated convolutional + Reed-Solomon code

## 🧪 How to compile
```bash
mkdir build
cd build
cmake ..
make
cd ..
./run_all.sh
```

## 🧪 How to clean the res directory
```bash
./run_all.sh clean
```

## ⚙️ Configuration file 
Configuration files are located in the ./conf directory.
Example configuration file:

```bash    
start_snr=-2          # Initial SNR value
end_snr=5             # Final SNR value; together with start_snr defines the simulation range
step_snr=0.25         # SNR step between start_snr and end_snr
output_file=results   # Output file prefix
puncturing_type=3/4   # Puncturing type for convolutional code

rs_encode=true        # Enable RS encoding
interleave=true       # Enable RS interleaving
scramble=true         # Enable scrambling
printing=false        # Print debug information
verbose=false         # Print detailed debug information
n_interleave=8        # Interleaver depth
dual_basis=true       # Use dual basis (as defined in CCSDS TM)

mode=rs_and_cc        # Mode: rs_and_cc, only_cc, or only_rs
```