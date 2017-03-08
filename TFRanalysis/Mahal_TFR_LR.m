%Decoding scripts

trash = input('Are you sure?','s');
clear all;
home


%% Mahalinobis Distance

addpath('/home/benmc/SuperWM/Decoding scripts');
addpath('/home/benmc/SuperWM/miscellanious scripts');
dat_dir =  '/home/benmc/SuperWM/Data';  % input data
savepath = '/home/benmc/SuperWM/Distances/TFR/';

if ~exist(savepath)
    mkdir(savepath)
end

pchans = {'PO7';'PO3';'PO4';'PO8';'O2';'O1'};
foi = 8:12;

sub_details(1).in_fname = {'s04_1','s04_2','s04_3'};
sub_details(2).in_fname = {'s06_1','s06_2','s06_3'};
sub_details(3).in_fname = {'s07_1','s07_2','s07_3'};
sub_details(4).in_fname = {'s08_1','s08_2','s08_3'};
sub_details(5).in_fname = {'s09_1','s09_2','s09_3'};
sub_details(6).in_fname = {'s10_1','s10_2','s10_3'};
sub_details(7).in_fname = {'s11_1'};%'s11_2','s11_3'
sub_details(8).in_fname = {'s12_1','s12_2','s12_3'};
sub_details(9).in_fname = {'s13_1','s13_2','s13_3'};
sub_details(10).in_fname = {'s14_1','s14_2','s14_3'};
sub_details(11).in_fname = {'s15_1','s15_2','s15_3'};
sub_details(12).in_fname = {'s16_1','s16_2','s16_3'};
sub_details(13).in_fname = {'s17_1','s17_2','s17_3'};
sub_details(14).in_fname = {'s18_1','s18_2','s18_3'};
sub_details(15).in_fname = {'s19_1','s19_2','s19_3'};
sub_details(16).in_fname = {'s21_1','s21_2','s21_3'};
sub_details(17).in_fname = {'s22_1','s22_2','s22_3'};
sub_details(18).in_fname = {'s23_1','s23_2','s23_3'};
sub_details(19).in_fname = {'s24_1','s24_2','s24_3'};
sub_details(20).in_fname = {'s25_1','s25_2','s25_3'};
sub_details(21).in_fname = {'s26_1','s26_2','s26_3'};
sub_details(22).in_fname = {'s27_1','s27_2','s27_3'};


nsubs = size(sub_details,2);

%% Define conds
%diff(1/10),pos(1/2/3/4),ss(1/2/4),acc(1/2);
%diff: 1=E,10=H
%pos: 1=upperR, 2=lowerR, 3=upperL, 4=lowerR

Conds{1} = 11; % ss1E
Conds{2} = 12; % ss2E

Conds{3} = 21; % ss1E
Conds{4} = 22; % ss2E

Conds{5} = 31; % ss1E
Conds{6} = 32; % ss2E

Conds{7} = 41; % ss1E
Conds{8} = 42; % ss2E

Conds{9} = 51; % ss1E
Conds{10} = 52; % ss2E

Conds{11} = 61; % ss1E
Conds{12} = 62; % ss2E

Conds{13} = 71; % ss1E
Conds{14} = 72; % ss2E

Conds{15} = [11 21 31 41 51 61 71]; % ss1E
Conds{16} = [12 22 32 42 52 62 72]; % ss2E

cond_pairs = [1 2;3 4;5 6;7 8;9 10;11 12;13 14;15 16];


%%
for pair = 1:length(cond_pairs)
    
    for c_sub = 1:nsubs
        
        nsess = size(sub_details(c_sub).in_fname,2);
        clear dat
        
        for c_sess = 1:nsess
            
            display([c_sub c_sess pair])
                         
            sessDir = fullfile(dat_dir,sub_details(c_sub).in_fname{c_sess});
            load(fullfile(sessDir,'condList.mat'))
            load(fullfile(sessDir,'TFRvbadTrials.mat'));
            load(fullfile(sessDir,'TFR_HP.mat'));
            
            ind = 1:length(condList);
            ind = setdiff(ind,vbadTrials);
            condList = condList(ind);
            CondList = condList;
            choi = find(ismember(tfdat.label,pchans));
            
            dat = log10(tfdat.powspctrm(ind,:,:,:));
          
            
%             btoi = tfdat.time>-.51&tfdat.time<-.1; % time window for baseline
%             bdat = nanmean(dat(:,:,:,btoi),4);
%             dat = dat - repmat(bdat,[1,1,1,size(dat,4)]);
            
            foi = tfdat.freq>=8 & tfdat.freq<=12;
            foi = find(foi);
            
            data = squeeze(nanmean(dat(:,choi,foi,:),3));
            data(isnan(data)) = 0; %because Mahal can't deal with NaNs
            
            lab_ind1 = find((ismember(CondList,Conds{cond_pairs(pair,1)}))==1);
            lab_ind2 = find((ismember(CondList,Conds{cond_pairs(pair,2)}))==1);          
            
            dat1 = data(lab_ind1,:,:);
            dat2 = data(lab_ind2,:,:);
            
            data = vertcat(dat1,dat2);
                        
            lab = vertcat(ones(length(lab_ind1),1),(2*ones(length(lab_ind2),1)));
            
            
            [dtw] = mahalFucTrialWise_new(data,lab);
            
            D{c_sub,c_sess} = dtw;
            
            %==================================================================
            %Do the GLM
            %Get RT data for correlat
            load(fullfile(sessDir,'SuperWM.mat')); %load subject data
        
            Accuracy = abs(output_data(ind,6));
            
            Acc = Accuracy([lab_ind1 lab_ind2]);
        
            des_mat = [Acc];
            clear B
            
                
            for tpt = 1:size(dtw,2)
                
                toFit = squeeze(dtw(:,tpt));
                
                B(tpt,:) = glmfit(des_mat,toFit);
                
            end
            
            Betas(c_sub,c_sess,:,:) = B;
            
        end
    
        
    end

    save([savepath 'TFR_LR_' num2str(pair)],'D');
    save([savepath 'TFR_LR_betas_' num2str(pair)],'Betas');
    
end

