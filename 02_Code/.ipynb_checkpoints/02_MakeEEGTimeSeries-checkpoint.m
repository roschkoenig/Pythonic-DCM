%% EdgeFC Make Edge and Node Networks

clc
clear all
cd('/media/aswinchari/DATA/EdgeFC');
load('channels.mat');
%% Generate list of patients

filelist = dir('sub*');

%% For each patient, load EDF, make timeseries and save (Day)

for a = 1:length(filelist)
    
    disp(strcat('Current File:',filelist(a).name,' Day'));
    
    % load EDF
    
    cfg             = [];
    filename        = strcat(string(filelist(a).name),'/',string(filelist(a).name(5:6)),'_Day.EDF');
    cfg.dataset     = filename{1};
    cfg.channel    = channels(a).loadchannelsnumber;
    cfg.dftfilter   = 'yes';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 256;
    cfg.continuous  = 'yes';
    if channels(a).Fs == 1024
        cfg.trl         = [1 1228800 0; 1228801 2457600 0; 2457601 3686400 0]; %3686400 = 1 hour at 1024Hz
    elseif channels(a).Fs == 2048
        cfg.trl         = [1 1228800 0; 1228801 2457600 0; 2457601 3686400 0; 3686401 4915200 0; 4915201 6144000 0; 6144001 7372800 0]; %[1 2457600 0; 2457601 4915200 0; 4915201 7372800 0];
    end
    
    eeg = ft_preprocessing(cfg);
    
    % resample
    
    cfg             = [];
    cfg.resamplefs  = 1024;
    cfg.detrend     = 'yes';
    
    eeg = ft_resampledata(cfg,eeg);
    
    % move from 6 to 3 trials if high Fs data
    
    if channels(a).Fs == 2048
        eeg.trial{1,1} = [eeg.trial{1,1} eeg.trial{1,2}];
        eeg.trial{1,2} = [eeg.trial{1,3} eeg.trial{1,4}];
        eeg.trial{1,3} = [eeg.trial{1,5} eeg.trial{1,6}];
        eeg.trial{1,4} = [];
        eeg.trial{1,5} = [];
        eeg.trial{1,6} = [];
    end 

    % save
    
    timeseries(1).timeseries = eeg.trial{1,1};
    timeseries(2).timeseries = eeg.trial{1,2};
    timeseries(3).timeseries = eeg.trial{1,3};
    
    save(strcat(string(filelist(a).name),'/',string(filelist(a).name(5:6)),'_Day.mat'),'timeseries','-v7.3');
    
end

%% For each patient, load EDF, make timeseries and save (Night)

for a = 1:length(filelist)
    
    disp(strcat('Current File:',filelist(a).name,' Night'));
    
    % load EDF
    
    cfg             = [];
    filename        = strcat(string(filelist(a).name),'/',string(filelist(a).name(5:6)),'_Night.EDF');
    cfg.dataset     = filename{1};
    cfg.channel    = channels(a).loadchannelsnumber;
    cfg.dftfilter   = 'yes';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 256;
    cfg.continuous  = 'yes';
    if channels(a).Fs == 1024
        cfg.trl         = [1 1228800 0; 1228801 2457600 0; 2457601 3686400 0]; %3686400 = 1 hour at 1024Hz
    elseif channels(a).Fs == 2048
        cfg.trl         = [1 1228800 0; 1228801 2457600 0; 2457601 3686400 0; 3686401 4915200 0; 4915201 6144000 0; 6144001 7372800 0]; %[1 2457600 0; 2457601 4915200 0; 4915201 7372800 0];
    end
    
    eeg = ft_preprocessing(cfg);
    
    % resample
    
    cfg             = [];
    cfg.resamplefs  = 1024;
    cfg.detrend     = 'yes';
    
    eeg = ft_resampledata(cfg,eeg);
    
    % move from 6 to 3 trials if high Fs data
    
    if channels(a).Fs == 2048
        eeg.trial{1,1} = [eeg.trial{1,1} eeg.trial{1,2}];
        eeg.trial{1,2} = [eeg.trial{1,3} eeg.trial{1,4}];
        eeg.trial{1,3} = [eeg.trial{1,5} eeg.trial{1,6}];
        eeg.trial{1,4} = [];
        eeg.trial{1,5} = [];
        eeg.trial{1,6} = [];
    end 

    % save
    
    timeseries(1).timeseries = eeg.trial{1,1};
    timeseries(2).timeseries = eeg.trial{1,2};
    timeseries(3).timeseries = eeg.trial{1,3};
    
    save(strcat(string(filelist(a).name),'/',string(filelist(a).name(5:6)),'_Night.mat'),'timeseries','-v7.3');
    
end
