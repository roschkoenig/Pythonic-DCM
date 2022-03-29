function dcm_path = pdcm_run_dcm(dcm_path)
load(dcm_path)
spm_dcm_csd(DCM)
disp('Inversion Complete')