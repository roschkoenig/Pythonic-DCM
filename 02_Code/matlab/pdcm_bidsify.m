% BIDS Setup Code copied adapted from https://www.fieldtriptoolbox.org/example/bids_eeg/
%==========================================================================
% Github version of data2bids downloaded here: https://github.com/fieldtrip/fieldtrip/blob/release/data2bids.m
% 22-03-23 RER
%
%%%%%%%%% Manual Input Required After Running %%%%%%%%%%%%%
%    - need to manually edit JSON file to add dataset name 
%    - runs are labelled as '1' instead of '01' (mne_bids assumes 2 char) 


% Import relevant packages (optional)
%--------------------------------------------------------------------------
F.toolboxes = ['/Users/roschkoenig/Dropbox/Research/0002 Tools/Matlab'];
F.addpath = {[F.toolboxes filesep 'ft-extensions'], [F.toolboxes filesep 'jsonlab']};
for a = 1:length(F.addpath), addpath(F.addpath{a}); end

sub = {'72'};
age = [nan];
sex = {[]};
ses = {'01'};
run = [1]; 

F.base = '/Volumes/GoogleDrive/My Drive/Research/2201_TVP-integ/01_Data';
F.raw  = [F.base filesep 'raw'];
F.bids = [F.base filesep 'bids'];

for subindex = 1:numel(sub) 
for si = 1:numel(ses)
for ri = 1:numel(run)
    cfg = []; 
    cfg.method = 'copy';
    cfg.datatype = 'eeg';

    % Specify the input filename here 
    %----------------------------------------------------------------------
    edf_dat = dir([F.raw filesep '*sub-' num2str(sub{subindex}) '*' num2str(ri, '%02.f') '*edf']);
    cfg.dataset = [F.raw filesep edf_dat(1).name]; 

    cfg.bidsroot    = F.bids; 
    cfg.sub         = sub{subindex};
    cfg.ses         = ses{si};
    cfg.run         = run(ri); 

    % information for the tsv file (optional, extendable); 
    %----------------------------------------------------------------------
    cfg.participants.age = age(subindex); 
    cfg.participants.sex = sex(subindex); 

    % json metadata
    %----------------------------------------------------------------------
    cfg.InstitutionName = 'Great Ormond Street Hospital';
    cfg.InstitutionalDepartmentName = 'Department of Clinical Neurophysiology';

    % provide mnemonic and long description of the task
    %----------------------------------------------------------------------
    cfg.TaskName = 'spont';
    cfg.TaskDescription = 'Spontaneous activity during daytime hours';

    % EEG-specific settings 
    %----------------------------------------------------------------------
    cfg.eeg.PowerLineFrequency = 50; 
    cfg.eeg.EEGReference = 'average';

    data2bids(cfg)
end
end
end
