import numpy as np
import nibabel as nib
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter

def main_apt(APT_dir_path, nOffsets, APT_file, AMP_file, PHA_file, APTw_out_file):
    # Parameters setup
    scanned_offset = {
        7: [np.inf, 383, -383, 447, -447, 510, -510],
    }[nOffsets]
    
    # Load NIfTI files
    apt_img = nib.load(f"{APT_dir_path}/{APT_file}.nii.gz")
    amp_img = nib.load(f"{APT_dir_path}/{AMP_file}.nii.gz")
    pha_img = nib.load(f"{APT_dir_path}/{PHA_file}.nii.gz")
    
    APT_data = apt_img.get_fdata()
    AMP_data = amp_img.get_fdata()
    PHA_data = pha_img.get_fdata()
    
    
    # Reshape data
    APT = APT_data
    print(APT.shape)
    M0 = APT[:, :, :, 0]

    B0_map = PHA_data  # Assuming PHA is already in Hz, otherwise calculate it with fsl_prepare_fieldmap
    
    # Create B0 null mask from magnitude image
    B0_null_masks = (AMP_data == 0)
    
    # APT processing parameters
    CESTfreq = 447
    input_params = {
        'APT1_rec': APT,
        'M0': M0,
        'b0_freq': B0_map,
        'scanned_offset': scanned_offset,
        'CESTfreq': CESTfreq,
        'B0_null_masks': B0_null_masks,
        'interpolation_choice': 'interp1'
    }
    
    # Calculate APTw
    APTw = apt_cal_linear_pn_separate(input_params)

    # Post processing
    APTw = np.where(AMP_data>0, APTw, -5)
    APTw = np.clip(APTw, -5, 5)
    # APTw = gaussian_filter(APTw, sigma=1)

    # Save result
    aptw_img = nib.Nifti1Image(APTw, apt_img.affine, apt_img.header)
    nib.save(aptw_img, f"{APT_dir_path}/{APTw_out_file}.nii.gz")

def apt_cal_linear_pn_separate(input):
    APT = input['APT1_rec']
    B0_map = input['b0_freq']
    scanned_offset = input['scanned_offset']
    B0_null_masks = input['B0_null_masks']
    M0 = input['M0']

    offsets = np.array(scanned_offset)
    valid_offsets = offsets[np.isfinite(offsets)]

    pos_offsets = valid_offsets[valid_offsets > 0]
    neg_offsets = valid_offsets[valid_offsets < 0]

    x, y, n_slices, _ = APT.shape
    APTw = np.zeros((x, y, n_slices))

    for slice_idx in range(n_slices):
        b0_slice = B0_map[:, :, slice_idx]
        mask = B0_null_masks[:, :, slice_idx]

        # Get indices after exclusion
        pos_indices = [np.where(offsets == o)[0][0] for o in pos_offsets] # [1, 3, 6]=[383, 447, 510]
        neg_indices = [np.where(offsets == o)[0][0] for o in neg_offsets] # [2, 4, 6]=[-383, -447, -510]

        pos_data = APT[:, :, slice_idx, pos_indices] # data sampled as positive frequences 
        neg_data = APT[:, :, slice_idx, neg_indices] # data samples at negative frequencies 

        # Create extended frequency ranges with margin extra margin (0.5ppm/64hz)
        if len(pos_offsets) > 1:
            pos_margin = pos_offsets[0] - pos_offsets[1]
            pos_start = pos_offsets[0] + pos_margin
            pos_end = pos_offsets[-1] - pos_margin
            pos_freq = pos_offsets
            pos_extended = np.arange(pos_start-1, pos_end+1, 1) #extend with 1hz step size, N_pos_extended = 257
        else:
            pos_extended = pos_offsets

        if len(neg_offsets) > 1:
            neg_margin = neg_offsets[0] - neg_offsets[1]
            neg_start = neg_offsets[0] + neg_margin
            neg_end = neg_offsets[-1] - neg_margin
            neg_freq = neg_offsets
            neg_extended = np.arange(neg_start+1, neg_end-1, -1)
        else:
            neg_extended = neg_offsets


        for i in range(x):
            for j in range(y):
                if not mask[i, j]:
                    shift = b0_slice[i,j]
                    interp_z_func_pos = interp1d(pos_freq, pos_data[i, j], bounds_error=False, fill_value='extrapolate')
                    z_new_pos = interp_z_func_pos(pos_extended)

                    interp_z_func_neg = interp1d(neg_freq, neg_data[i, j], bounds_error=False, fill_value='extrapolate')
                    z_new_neg = interp_z_func_neg(neg_extended)

                    shift = np.round(shift)
                    new_pos_pos = 447 + shift
                    new_neg_pos = -447 + shift

                    pos_index = np.where(pos_extended==new_pos_pos)[0]
                    neg_index = np.where(neg_extended==new_neg_pos)[0]

                    if len(pos_index>0) and len(neg_index>0):
                        mtr_n_447 = z_new_neg[neg_index]
                        mtr_p_447 = z_new_pos[pos_index]

                        print((mtr_n_447 - mtr_p_447))
                        print(M0[i, j, slice_idx])
                        APTw[i, j, slice_idx] = 100*(mtr_n_447 - mtr_p_447)/(M0[i, j, slice_idx] + 1e-6)

    return APTw

# Example usage:
main_apt("processed", 7, "apt_mcf", "b0_mag_brain", "fmap_hz", "13_APTw")
