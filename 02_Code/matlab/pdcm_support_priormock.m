fpath = '/Volumes/GoogleDrive/My Drive/Research/2201_TVP-integ/03_Output/preproc/sub-72_seeg_segmented.mat';
load(fpath); 

D = zeros(size(seeg.data,2));
for d = 1:size(seeg.data,1)
    D = D+corr(squeeze(seeg.data(1,:,:))');
end
D = D/d; 

%% Transform
x = 0:0.01:1; 
p = [10 0.6];
subplot(121), plot(x, sig(x,p));

tD          = zeros(size(D));
tD(D>0)     = D(D>0); 
tD          = sig(tD,p); 
tD          = exp(tD)-1; 
tD(tD<0.05) = 0; 
subplot(122), imagesc(tD), axis square

function y = sig(x, p), y = 1./(1+exp(-p(1)*(x-p(2)))); end