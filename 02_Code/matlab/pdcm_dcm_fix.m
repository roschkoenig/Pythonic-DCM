function DCM = pdcm_dcm_fix(DCM)

    % Reformat some DCM structure fields
    %--------------------------------------------------------------------------
    DCM.options.Fdcm = double(DCM.options.Fdcm); 
    if ischar(DCM.Sname), DCM.Sname = cellstr(DCM.Sname); end 
    if ~iscell(DCM.A)
        for a = 1:size(DCM.A,1)
            A{a} = squeeze(DCM.A(a,:,:));
        end
        DCM.A = A; 
    end
    if isempty(DCM.C), DCM.C = sparse(length(DCM.A{1}),0); end 
    DCM.M.nograph = 1; 