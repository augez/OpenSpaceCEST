## Overview
This Python script processes Amide Proton Transfer (APT) Chemical Exchange Saturation Transfer (CEST) MRI data to generate APT-weighted (APTw) images. The implementation includes B0 inhomogeneity correction, interpolation of saturation offsets, and symmetric subtraction analysis.

## Requirements
- Python >= 3.9
- Required packages:
  - numpy
  - nibabel
  - scipy

## Usage
```python
main_apt(APT_dir_path, nOffsets, APT_file, AMP_file, PHA_file, APTw_out_file)
```

### Parameters:
- `APT_dir_path`: Directory path containing input NIfTI files
- `nOffsets`: Number of frequency offsets (currently supports 7 offsets)
- `APT_file`: Filename prefix for (motion corrected) APT CEST data (without extension)
- `AMP_file`: Filename prefix for magnitude image (without extension)
- `PHA_file`: Filename prefix for B0 field map in Hz (without extension)
- `APTw_out_file`: Output filename prefix for APT-weighted image (without extension)

### Example:
```python
main_apt("processed", 7, "apt_mcf", "b0_mag_brain", "fmap_hz", "13_APTw")
```

## Input Files
1. **APT CEST Data** (`{APT_file}.nii.gz`):
   - 4D NIfTI file containing CEST Z-spectra data
   - Dimensions: [X, Y, Slices, Offsets]
   - Offsets must include M0 image (first volume) and symmetric frequency offsets

2. **Magnitude Image** (`{AMP_file}.nii.gz`):
   - 3D NIfTI file used for masking
   - Dimensions: [X, Y, Slices]

3. **B0 Field Map** (`{PHA_file}.nii.gz`):
   - 3D NIfTI file containing B0 inhomogeneity map in Hz
   - Dimensions: [X, Y, Slices]

## Processing Steps
1. **Data Loading**:
   - Loads APT CEST data, magnitude image, and B0 field map
   - Extracts M0 image from first volume of APT data
   - Creates B0 null mask from magnitude image

2. **APT Calculation**:
   - Separates positive and negative frequency offsets
   - For each voxel:
     - Interpolates Z-spectra using linear interpolation
     - Applies B0 correction by shifting frequency offsets
     - Calculates APTw value: 100*(MTR_neg_447 - MTR_pos_447)/M0
   - Uses symmetric offsets at ±447 Hz for APT contrast

3. **Post-processing**:
   - Sets APTw to -5 where magnitude is 0
   - Clips values to range [-5, 5]
     
## Output
- **APTw Image** (`{APTw_out_file}.nii.gz`):
  - 3D NIfTI file containing APT-weighted values
  - Dimensions: [X, Y, Slices]
  - Values in percentage units (%)

## Technical Details
- **Frequency Offsets**: Currently configured for 7 offsets:
  - [∞, 383, -383, 447, -447, 510, -510] Hz
- **Interpolation**: Uses linear interpolation (scipy.interpolate.interp1d)
- **B0 Correction**: Shifts frequency offsets based on B0 field map
- **Symmetric Analysis**: Calculates APT contrast using symmetric offsets at ±447 Hz
- **Voxel Processing**: Processes each voxel individually with B0 correction

## How to Reference
When using this software in your research, please cite it as follows:

> Essed, R. (2025). OpenSpaceCEST: A processing pipeline for SPACE APT-CEST MRI data [Computer software]. GitHub repository. https://github.com/augez/OpenSpaceCEST

For example, in a bibliography:

Essed, R. (2025). OpenSpaceCEST: A processing pipeline for SPACE APT-CEST MRI data [Computer software]. GitHub repository. https://github.com/augez/OpenSpaceCEST

If you publish work using this software, please include a citation to the GitHub repository in your methods section and inform the maintainers by creating an issue in the repository.

## Notes
- The B0 field map must be in Hz (not radians)
- The code assumes the first volume of APT data is the M0 image
- For optimal results, ensure input images are properly co-registered
- The implementation is optimized for SPACE APT-CEST acquisition sequences
