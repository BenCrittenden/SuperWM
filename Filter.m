
clear all;
close all;
clc;
%% pre-processing script

% Add and remove appropriate paths

% rmpath(genpath(['/home/benmc/matlab/spm8/external/fieldtrip']));
% rmpath(genpath('/opt/analysis/spm8/'))

addpath('/home/benmc/matlab/');
eeglabpath = '/home/benmc/matlab/eeglab/eeglab11_0_5_4b';
addpath(genpath(eeglabpath));

dat_dir =  '/home/benmc/SuperWM/Data';  % input data

% Define subject files
%===Add subjects as they come in, don't forget to change the number X in
%sub_details(X)
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


%======================

%===Select subjects to run
nsubs=length(sub_details);
exclude_subs = [20,21,46,47,48];
dosubs = [1:nsubs];
dosubs = setdiff(dosubs,exclude_subs);


for sub = dosubs
    subDir = fullfile(dat_dir ,sub_details{sub});
    
    % Add EEG lab to path
    if ~exist('ALLEEG','var')
        ALLEEG = eeglab;
    end
    

        fname_in = 'rawEEG.cnt';
        fname_out = 'EEG_filtHP.set';

    
%     catch if this subject has already been run
    if 0%exist(fullfile(subDir,fname_out),'file')
        fprintf('filter already done\n')
        EEG = pop_loadset(fname_out, subDir);
    else
        
        %Load the data, specifity the format of the data
        EEG = pop_loadcnt(fullfile(subDir,fname_in),'dataformat', 'int32');
        
        %Specify the subjects files
        substrg = sub_details{sub};
        EEG.filename = sub_details{sub};
        EEG.setname  = sprintf('S%s',substrg);
        
        %store the EEG datset
        [ALLEEG,EEG] = eeg_store(ALLEEG,EEG,1);
        eeglab redraw
        
        %% Reorder EEG/EOG channels
        %close all
        
        chaneeg = [30 64 28 62 63 59 60 29 61 27 58 21 55 22 56 23 25 57 26 54 20 ...
            51 16 52 17 24 53 19 49 15 50 8  43 10 18 47 13 48 14 39 6  42 ...
            9  11 46 12 35 4  37 5  40 7  45 34 36 38 44 1  2 3];
        chaneog = [65 66];
        chanemg = [31 41];
        chanref = [32];
        chanrmv = [33];
        chanall = [chaneeg chaneog chanemg chanref];
        
        
        EEG.chanlocs = EEG.chanlocs(chanall);
        EEG.data     = EEG.data(chanall,:);

        
        EEG.setname = sprintf('S%s reordered',substrg);
        
        
        [ALLEEG,EEG] = eeg_store(ALLEEG,EEG,1);
        eeglab redraw
        %%
        chaneeg = 01:60;
        chaneog = 61:62;
        chanemg = 63:64;
        chanref = 65;
        
        % Rereference EEG/EOG channels
        EEG.data              = EEG.data-repmat(mean(EEG.data(chanref,:),1),EEG.nbchan,1);
        EEG.chanlocs(chanref) = [];
        EEG.data(chanref,:)   = [];
        EEG.nbchan            = EEG.nbchan-length(chanref);
        EEG.setname = sprintf('S%s rereferenced',substrg);
        
        [ALLEEG,EEG] = eeg_store(ALLEEG,EEG,1);
        eeglab redraw
        
        
        EEG.chanlocs(chanemg(1)).labels = 'EMGL';
        EEG.chanlocs(chanemg(2)).labels = 'EMGR';
        
        [ALLEEG,EEG] = eeg_store(ALLEEG,EEG,1);
        eeglab redraw
        % Resample to 250 Hz
        if EEG.srate ~= 250
            EEG = pop_resample(EEG, 250);
            EEG.setname = sprintf('S%s resampled',substrg);
            [ALLEEG,EEG] = eeg_store(ALLEEG,EEG,1);
        end
        eeglab redraw
        %% Band-pass filter EEG channels
        do_filtering = true;
        if do_filtering
            hpfreq = 0.03; % high-pass frequency (Hz)
            lpfreq = 45; % low-pass frequency (Hz)
            
            data = EEG.data(chaneeg,:);
            data = eegfilt(data,EEG.srate,hpfreq,0);
            data = eegfilt(data,EEG.srate,0,lpfreq);
            EEG.data(chaneeg,:) = data;

            
            data = EEG.data(chaneog,:);
            data = eegfilt(data,EEG.srate,hpfreq,0);
            data = eegfilt(data,EEG.srate,0,lpfreq);
            EEG.data(chaneog,:) = data;
            
         
            data = EEG.data(chanemg,:);
            data = eegfilt(data,EEG.srate,hpfreq,0);
            data = eegfilt(data,EEG.srate,0,lpfreq);
            EEG.data(chanemg,:) = data;
            
            EEG.setname = sprintf('S%s filtered',substrg);
            
            [ALLEEG,EEG] = eeg_store(ALLEEG,EEG,1);
            eeglab redraw
        end
        
        
        
        
        %% save filtered data
        pop_saveset( EEG, 'filename',fname_out,'filepath',subDir);
        clear data
    end
    
 end