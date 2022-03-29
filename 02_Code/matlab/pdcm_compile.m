pdcm_path = which('pdcm');
build = compiler.build.standaloneApplication(pdcm_path);
compiler.package.installer(build);

%% Not currently working because of some filenameing issues
%--------------------------------------------------------------------------
[p,f] = fileparts(pdcm_path);
installfrom = [p filesep 'pdcminstaller' filesep 'MyAppInstaller.app'];
installto   = [p filesep 'pdcm_standalone'];
system(['sudo ' '"' installfrom '"' ' -agreeToLicense yes' ...
        ' -destinationFolder' '"' installto '"'])

%%
cfg.task     = 'run_dcm';
cfg.specfile = '/Volumes/GoogleDrive/My Drive/Research/2201_TVP-integ/03_Output/dcm/sub-72_dcm_spec.mat';
cfg.dcm_path = '/Volumes/GoogleDrive/My Drive/Research/2201_TVP-integ/03_Output/dcm/DCM_sub-72.mat';