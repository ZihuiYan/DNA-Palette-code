# DNA Palette Code

This project implements a single-threaded execution of the DNA Palette code, including encoding and decoding for DICOM files, as well as simulating random data.

## Example

- **Encoding for DICOM:** Encode a folder of DICOM files into oligos using the following command:
  ```bash
  python DICOM_encoder.py
- **Decoding for DICOM:** Decode a text file containing sequencing reads back into DICOM files using the following command:
  ```bash
  One need to first download a partially sampled dataset of sequencing reads (https://doi.org/10.6084/m9.figshare.25567545.v1) in this directory and run the demo.
  python DICOM_decoder.py
- **Testing for Random Data:** Simulate the DNA Palette code using random data with the following command:
  python example_test.py

- **Parallel Processes Testing with Random Data:**
  DNA Palette code_multiProcesses.zip

## Installation

Before running the code, make sure to install the required dependencies using the following command:

  numpy
  reedsolo
