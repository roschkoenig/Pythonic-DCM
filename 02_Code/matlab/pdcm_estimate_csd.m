function dcm_path = pdcm_estimate_csd(specfile)
load(specfile)
DCM             = pdcm_dcm_fix(DCM);
DCM.xY.Dfile    = pdcm_filemaker(specfile);
DCM             = spm_dcm_erp_data(DCM);
DCM             = spm_dcm_erp_dipfit(DCM, 1); 
DCM             = spm_dcm_csd_data(DCM);
DCM.options.DATA        = 0; 
save(DCM.name, 'DCM')
disp('CSD estimation complete')
dcm_path = DCM.name;