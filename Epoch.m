
clear all;
close all;
clc;
%% pre-processing script

addpath('/home/benmc/matlab/');
eeglabpath = '/home/benmc/matlab/eeglab/eeglab11_0_5_4b';
rmpath(genpath(eeglabpath));

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

% Select subjects to run
nsubs=length(sub_details);
exclude_subs = [20,21,46,47,48];
dosubs = [1:nsubs];
dosubs = setdiff(dosubs,exclude_subs);

%creat matrix of event codes
M = zeros(7,2);
M(1,:) = M(1,:) + 10; 
M(2,:) = M(2,:) + 20;  
M(3,:) = M(3,:) + 30; 
M(4,:) = M(4,:) + 40; 
M(5,:) = M(5,:) + 50; 
M(6,:) = M(6,:) + 60; 
M(7,:) = M(7,:) + 70; 
M(:,1) = M(:,1) + 1; % left
M(:,2) = M(:,2) + 2; % right

selecEvents = reshape(M,numel(M),1) + 100; %+1000 because first trigger is for auditory cue
nEvents = length(selecEvents);


for sub = dosubs
    %%
        
    subDir = fullfile(dat_dir,sub_details{sub});
    load(fullfile(subDir,'SuperWM.mat'));
    if ~exist('ALLEEG','var')
        ALLEEG = eeglab;
    end
     
    fname_in = 'EEG_filtHP.set';
    fname_out = ['HP_Epoched_Delay_EEG.set'];
    if 0%exist(fullfile(subDir,fname_out),'file')
        fprintf('filter already done\n')
        EEG = pop_loadset(fname_out, subDir);
    else
        epoch

       % first convert back to eeglab
        EEG = pop_loadset(fname_in, subDir);

        if size(EEG.event,2)~=size(output_data,1)*3
            
            if strcmp('s08_2',sub_details{sub})
                %started the EEG recording after first trial, so remove
                %first trial from behavioural data
                
                output_data = output_data(2:end,:);
                
            elseif strcmp('s22_2',sub_details{sub})
                %started the EEG recording after first trial, so remove
                %first trial from behavioural data
                
                output_data = output_data(2:end,:);
                
            else
            
                for tmp = 1:size(EEG.event,2)
                    EEGcode(tmp) = str2num(EEG.event(1,tmp).type);
                end
                
                
                display('=======================')
                display(' ')
                display(' ')
                display(sub_details{sub})
                display(size(EEG.event,2))
                display(size(output_data,1)*2)
                display(' ')
                display(' ')
                display('=======================')
                error('holy crap!!')
            
            end
        end
        
        
        nTrials = length(output_data);
        ev=1;
        for i=1:nTrials
                % difficult, setsize, location
                code = (output_data(i,2)*10) + output_data(i,10);
                EEG.event(ev).type = code;
                EEG.urevent(ev).type = code;
                ev=ev+1;
                EEG.event(ev).type = code+100; %code for probes trigger, which always followed the stimulus trigger (not needed by me)
                EEG.urevent(ev).type = code+100; %code for probes (not needed by me)
                ev=ev+1;
                EEG.event(ev).type = code+1000; %code for probes trigger, which always followed the stimulus trigger (not needed by me)
                EEG.urevent(ev).type = code+1000; %code for probes (not needed by me)
                ev=ev+1;
        end

        EEG = pop_epoch(EEG, num2cell(selecEvents),[-1 3.5]);
        EEG.setname = fname_out;
        % save epoched data
        pop_saveset(EEG, 'filename',fname_out,'filepath',subDir);
   

    end

    %% find very bad trials
    fname_in = fname_out;
    if exist(fullfile(subDir,'vbadTrials.mat'),'file')
        fprintf('vbadTrials already found\n')
    else
        EEG = pop_loadset(fname_in, subDir);  
         tmpEEG = EEG;
         tmpEEG.data(eeg_chaninds(tmpEEG,{'EMGL'}),:,:)        = zeros(size(tmpEEG.data(eeg_chaninds(tmpEEG,{'EMGL'}),:,:) ));
         tmpEEG.data(eeg_chaninds(tmpEEG,{'EMGR'}),:,:)        = zeros(size(tmpEEG.data(eeg_chaninds(tmpEEG,{'EMGR'}),:,:) ));
        
        if exist([dat_dir '/vbadTrials.mat']); 
            load([dat_dir '/vbadTrials.mat']);
            temp = zeros(tmpEEG.trials,1);
            temp(vbadTrials) = 1;
            tmpEEG.reject.rejmanual = temp';
            tmpEEG.reject.rejmanualE = zeros(tmpEEG.nbchan,tmpEEG.trials);
        end
        
        pop_eegplot(tmpEEG,1,1,0);
        input('press any key to continue')
        vbadTrials = find(ALLEEG.reject.rejmanual);
        save(fullfile(subDir,'vbadTrials.mat'),'vbadTrials')
        clear tmpEEG;
    end
    %%
    condList = [];
    N = 3; %select every third element (because there are three event types)
    ind = [1:N:(length(EEG.epoch).*3)]; %start at second event type
    
    for i=1:length(EEG.epoch)
        condList(i) = EEG.urevent(ind(i)).type;
%         ind=ismember(EEG.epoch(i).eventlatency,[0]);
%         condList(i) = EEG.epoch(i).eventtype(ind);
    end
    
  save(fullfile(subDir,'condList.mat'),'condList')
    
end
    
   