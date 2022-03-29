function fname = pdcm_filemaker(specfile)
    load(specfile);

    % Convert data to SPM file and add pointer to DCM structure 
    %--------------------------------------------------------------------------
    % Define trial structure 
    for d = 1:size(data,1)
        ftdata.trial{d} = squeeze(data(1,:,:));
        ftdata.time{d}  = time;
    end
    ftdata.label = cellstr(labels); 
    fname   = [Fbase filesep '01_Data' filesep 'matlab' filesep 'sub-' subject '_SEEG'];
    D       =  spm_eeg_ft2spm(ftdata, fname); 
    
    % Label conditions
    conds = cellstr(conds);
    for d = 1:size(D,3), D = conditions(D, d, conds{d}); end
    save(D)
