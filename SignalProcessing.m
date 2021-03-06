
clear all;
close all;
clc;
%% pre-processing script

addpath('/home/benmc/matlab/');
eeglabpath = '/home/benmc/matlab/eeglab11_0_5_4b';
rmpath(genpath(eeglabpath));
addpath('/home/benmc/matlab/fieldtrip-20130507');
addpath('/home/benmc/matlab/shadedErrorBar');

dat_dir =  '/home/benmc/SuperWM/Data';  % input data


sub_details(1) = {'s04_1'};
sub_details(2) = {'s04_2'};
sub_details(3) = {'s04_3'};
sub_details(4) = {'s06_1'};
sub_details(5) = {'s06_2'};
sub_details(6) = {'s06_3'};
sub_details(7) = {'s07_1'};
sub_details(8) = {'s07_2'};
sub_details(9) = {'s07_3'};
sub_details(10) = {'s08_1'};
sub_details(11) = {'s08_2'};
sub_details(12) = {'s08_3'};
sub_details(13) = {'s09_1'};
sub_details(14) = {'s09_2'};
sub_details(15) = {'s09_3'};
sub_details(16) = {'s10_1'};
sub_details(17) = {'s10_2'};
sub_details(18) = {'s10_3'};
sub_details(19) = {'s11_1'};
sub_details(20) = {'s11_2'};
sub_details(21) = {'s11_3'};
sub_details(22) = {'s12_1'};
sub_details(23) = {'s12_2'};
sub_details(24) = {'s12_3'};
sub_details(25) = {'s13_1'};
sub_details(26) = {'s13_2'};
sub_details(27) = {'s13_3'};
sub_details(28) = {'s14_1'};
sub_details(29) = {'s14_2'};
sub_details(30) = {'s14_3'};
sub_details(31) = {'s15_1'};
sub_details(32) = {'s15_2'};
sub_details(33) = {'s15_3'};
sub_details(34) = {'s16_1'};
sub_details(35) = {'s16_2'};
sub_details(36) = {'s16_3'};
sub_details(37) = {'s17_1'};
sub_details(38) = {'s17_2'};
sub_details(39) = {'s17_3'};
sub_details(40) = {'s18_1'};
sub_details(41) = {'s18_2'};
sub_details(42) = {'s18_3'};
sub_details(43) = {'s19_1'};
sub_details(44) = {'s19_2'};
sub_details(45) = {'s19_3'};
sub_details(46) = {'s20_1'};
sub_details(47) = {'s20_2'};
sub_details(48) = {'s20_3'};
sub_details(49) = {'s21_1'};
sub_details(50) = {'s21_2'};
sub_details(51) = {'s21_3'};
sub_details(52) = {'s22_1'};
sub_details(53) = {'s22_2'};
sub_details(54) = {'s22_3'};
sub_details(55) = {'s23_1'};
sub_details(56) = {'s23_2'};
sub_details(57) = {'s23_3'};
sub_details(58) = {'s24_1'};
sub_details(59) = {'s24_2'};
sub_details(60) = {'s24_3'};
sub_details(61) = {'s25_1'};
sub_details(62) = {'s25_2'};
sub_details(63) = {'s25_3'};
sub_details(64) = {'s26_1'};
sub_details(65) = {'s26_2'};
sub_details(66) = {'s26_3'};
sub_details(67) = {'s27_1'};
sub_details(68) = {'s27_2'};
sub_details(69) = {'s27_3'};


nsubs=length(sub_details);
exclude_subs = [20,21,46,47,48];
dosubs = [1:nsubs];
dosubs = setdiff(dosubs,exclude_subs);
                                                  

for sub = dosubs
    subDir = fullfile(dat_dir,sub_details{sub});    
    load(fullfile(subDir,'condList.mat'))
    
    display(sub)    

    fname_in = ['HP_Epoched_Delay_EEG.set'];
    fname_out = ['ERP_HP.mat'];
    
    if 0%exist(fullfile(subDir,fname_out),'file')
        fprintf('convert trials to fieldtrip already done\n')
        load(fullfile(subDir,fname_out));
    else
        cfg = [];
        cfg.dataset = fullfile(subDir,fname_in);
        cfg.continuous = 'no';
        cfg.baselinewindow = [-0.5 -0.1];
        cfg.demean = 'no'; %perform baseline correction
        cfg.detrend = 'no';
        dat = ft_preprocessing(cfg);
        
        save(fullfile(subDir,fname_out),'dat');
    end
    
    fname_out = ['TFR_HP.mat'];
    
    if 0%exist(fullfile(subDir,fname_out),'file')
        fprintf('convert trials to fieldtrip already done\n')
        
    else
        cfg              = [];
        %cfg.trials       = ismember(condList,ss2);
        cfg.output       = 'pow';
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:1:30;                         % analysis 2 to 30 Hz in steps of 1 Hz
        cfg.t_ftimwin    = 5./cfg.foi;                     % length of time window = 0.5 sec
        cfg.toi          = dat.time{1}(1):0.05:dat.time{1}(end);                 
        cfg.keeptrials   = 'yes';
        tfdat = ft_freqanalysis(cfg, dat);
        
        tfdat.powspctrm = single(tfdat.powspctrm);
        
        display('saving tfr')
        
        save(fullfile(subDir,fname_out),'tfdat','-v7.3');
    end

end

