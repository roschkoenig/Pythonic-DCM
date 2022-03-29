function pdcm(cfg)

disp('SPM initialising')
spm('defaults', 'eeg')
disp('SPM initialised')

% If called with path to configure file
%--------------------------------------------------------------------------
if ischar(cfg), load(cfg); end 

% Check whether this is running in principle
%--------------------------------------------------------------------------
if strcmp(cfg.task, 'test'), disp('The App is working'), end

% Call cross-spectral density estimation (will update DCM path)
%--------------------------------------------------------------------------
if strcmp(cfg.task,'estimate_csd')
    dcm_path     = pdcm_estimate_csd(cfg.specfile);
    cfg.dcm_path = dcm_path;
end

% Call DCM inversion (will update DCM path)
%--------------------------------------------------------------------------
if strcmp(cfg.task,'run_dcm'), dcm_path = pdcm_run_dcm(cfg.dcm_path); end 