
fpath = '/Volumes/GoogleDrive/My Drive/Research/2201_TVP-integ/03_Output/preproc/sub-72_seeg_segmented.mat';
load(fpath); 

D = zeros(size(seeg.data,2));
for d = 1:size(seeg.data,1)
    D = D+corr(squeeze(seeg.data(1,:,:))');
end
D = D/d; 
