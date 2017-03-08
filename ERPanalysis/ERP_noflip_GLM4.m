%%
%The purpose of this GLM is to see where the CDA correlates with Accuracy.
%Done seperately for each set-size.

trash = input('Are you sure?','s');

clear all
close all

%addpath(genpath('/home/benmc/matlab/fieldtrip-20130507'));
addpath('/Users/benc/MATLAB/');

addpath('/Users/benc/MATLAB/fieldtrip-20130507');
addpath('/Users/benc/MATLAB/shadedErrorBar');

dat_dir =  '/Users/benc/Experiments/SuperWM/Data';  % input data


sub_details(1).in_fname = {'s04_1','s04_2','s04_3'};
sub_details(2).in_fname = {'s06_1','s06_2','s06_3'};
sub_details(3).in_fname = {'s07_1','s07_2','s07_3'};
sub_details(4).in_fname = {'s08_1','s08_2','s08_3'};
sub_details(5).in_fname = {'s09_1','s09_2','s09_3'};
sub_details(6).in_fname = {'s10_1','s10_2','s10_3'};
sub_details(7).in_fname = {'s13_1','s13_2','s13_3'};
sub_details(8).in_fname = {'s14_1','s14_2','s14_3'};
sub_details(9).in_fname = {'s15_1','s15_2','s15_3'};
sub_details(10).in_fname = {'s16_1','s16_2','s16_3'};
sub_details(11).in_fname = {'s17_1','s17_2','s17_3'};
sub_details(12).in_fname = {'s18_1','s18_2','s18_3'};
sub_details(13).in_fname = {'s19_1','s19_2','s19_3'};
sub_details(14).in_fname = {'s21_1','s21_2','s21_3'};
sub_details(15).in_fname = {'s22_1','s22_2','s22_3'};
sub_details(16).in_fname = {'s23_1','s23_2','s23_3'};
sub_details(17).in_fname = {'s24_1','s24_2','s24_3'};
sub_details(18).in_fname = {'s25_1','s25_2','s25_3'};
sub_details(19).in_fname = {'s26_1','s26_2','s26_3'};
sub_details(20).in_fname = {'s27_1','s27_2','s27_3'};


Conds = [11 21 31 41 51 61 71 12 22 32 42 52 62 72];

nsubs=length(sub_details);
All_WM = 1:nsubs;
High_WM = [1:7 13 18]; %Used in EEG
Low_WM = [8:12 14:17 19 20]; %Used in EEG
dosubs = All_WM;


AllSubDat = [];
for sub = dosubs
    AllSessDat = [];
    n_trials_in_sess = [];
    for sess = 1:length(sub_details(sub).in_fname)
        
        clear data
        
        fprintf(['Doing ' sub_details(sub).in_fname{sess} '\n'])
        sessDir = fullfile(dat_dir,sub_details(sub).in_fname{sess});
        load(fullfile(sessDir,'condList.mat'))
        load(fullfile(sessDir,'ERPvbadTrials.mat'));
        load(fullfile(sessDir,'ERP_HP.mat'));
        ind = 1:length(condList);
        ind = setdiff(ind,vbadTrials);
        n_trials_in_sess(sess) = length(ind);
        condList = condList(ind);
        dat.trial = dat.trial(ind);
        
        data = nan(size(dat.trial,2),size(dat.trial{1},1),size(dat.trial{1},2));
        for ct = 1:size(dat.trial,2)
            data(ct,:,:) = dat.trial{ct};
        end
        
        
        %baseline
        btoi = find(dat.time{1}<0 & dat.time{1}>-0.005);
        bdat = mean(data(:,:,btoi),3);
        data = data - repmat(bdat,[1,1,size(data,3)]); % remove baseline
        
        
        %Run GLM
        load(fullfile(sessDir,'SuperWM.mat')); %load subject data
        
        Acc = abs(output_data(ind,6));
        SS = output_data(ind,2);
        
        left_chans = {'O1','PO7','PO3'};
        right_chans =  {'O2','PO8','PO4'};
        
        left_ind = find(ismember(dat.label,left_chans));
        right_ind = find(ismember(dat.label,right_chans));
        
        for tpt = 1:size(data,3)
            
            for ss=1:7
                
                switch ss
                    case 1
                        sstrl = [11 12];
                    case 2
                        sstrl = [21 22];
                    case 3
                        sstrl = [31 32];
                    case 4
                        sstrl = [41 42];
                    case 5
                        sstrl = [51 52];
                    case 6
                        sstrl = [61 62];
                    case 7
                        sstrl = [71 72];
                end
                
                sstrl = find(ismember(condList,sstrl));
                
                for trl = sstrl
                    
                    if ismember(condList(trl),[11 21 31 41 51 61 71])
                        cl(trl,:,:) = data(trl,right_ind,tpt);
                        il(trl,:,:) = data(trl,left_ind,tpt);
                    else
                        cl(trl,:,:) = data(trl,left_ind,tpt);
                        il(trl,:,:) = data(trl,right_ind,tpt);
                    end
                    
                end
                
                ci = cl - il;
                
                cl = squeeze(mean(cl,3));
                il = squeeze(mean(il,3));
                ci = squeeze(mean(ci,3));
                
                ci = ci(sstrl);
                
                des_mat = horzcat(ci); %,il,ci);
                
                B_ac(tpt,ss,:) = glmfit(Acc(sstrl),des_mat);
                
                clear cl il ci
                
            end
            
        end
        
        AllSessBeta_ac(sess,:,:,:) = B_ac;
        
        clear B_ss B_ac
        
    end
    
    AllSubBeta_ac(sub,:,:,:) = squeeze(nanmean(AllSessBeta_ac,1));
    
end



%%


trange = [-100 3000];
toi = find((dat.time{1}*1000)>=trange(1)&(dat.time{1}*1000)<=trange(2));

dosubs =High_WM;




%%

smf = 50;
ss=7;
times = dat.time{1}(toi);

AllSubBeta = AllSubBeta_ac;

Hbeta_ci = gsmooth(squeeze(nanmean(AllSubBeta(High_WM,toi,ss,2),1)),smf);
Lbeta_ci = gsmooth(squeeze(nanmean(AllSubBeta(Low_WM,toi,ss,2),1)),smf);



col = [0 0 1];

figure,

hold on

plot(times, Hbeta_ci,'-', 'LineWidth',3,...
    'Color',[1 0 0]);

plot(times, Lbeta_ci,'-', 'LineWidth',3,...
    'Color',[0 0 1]);



% axes = [-0.5 3 -0.02 .01];
% axes = [-0.5 3 -0.8 1];
axes = [-inf inf -inf inf];
axis(axes)

hold off


%%

test_dat1 = squeeze(AllSubBeta(High_WM,toi,2));
test_dat2 = squeeze(AllSubBeta(Low_WM,toi,2));

%smooth data
for s = 1:size(test_dat1,1)
    
    test_dat1(s,:) = gsmooth(test_dat1(s,:),smf);
    
end

for s = 1:size(test_dat2,1)
    
    test_dat2(s,:) = gsmooth(test_dat2(s,:),smf);
    
end

clear H pval tp

for tp = 1:size(test_dat1,2)
    
    x = test_dat1(:,tp);
    y = test_dat2(:,tp);
    
    z(tp) = mean(x) - mean(y);
    
    [H(tp) pval(tp)] = ttest2(x,y,'tail','both');
    
end

figure,imagesc(H)


