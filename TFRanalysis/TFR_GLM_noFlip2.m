trash = input('Are you sure?','s');

clear all

addpath('/home/benmc/matlab/');
eeglabpath = '/home/benmc/matlab/eeglab11_0_5_4b';
rmpath(genpath(eeglabpath));
addpath('/home/benmc/matlab/fieldtrip-20130507');
addpath('/home/benmc/matlab/shadedErrorBar');
addpath('/home/benmc/AttMem/Miscellanious');
addpath('/home/benmc/matlab/paruly/paruly');

dat_dir =  '/home/benmc/SuperWM/Data';  % input data

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

Conds = [11 21 31 41 51 61 71 ; 12 22 32 42 52 62 72];

    
nsubs=length(sub_details);
All_WM = 1:nsubs;
High_WM = [1 2 3 4 5 6 7 9 15 20]; %3 7
Low_WM = [10:14 16:19 21 22]; %12
dosubs = All_WM;


%%          
          
AllSubDat = [];
for sub = dosubs
    AllSessDat = [];
    n_trials_in_sess = [];
    
    clear data
    
    for sess = 1:length(sub_details(sub).in_fname)
        fprintf(['Doing ' sub_details(sub).in_fname{sess} '\n'])
        sessDir = fullfile(dat_dir,sub_details(sub).in_fname{sess});
        load(fullfile(sessDir,'condList.mat'))

        load(fullfile(sessDir,'TFRvbadTrials.mat'));
        load(fullfile(sessDir,'TFR_HP.mat'));
        ind = 1:length(condList);
        ind = setdiff(ind,vbadTrials);
        n_trials_in_sess(sess) = length(ind);
        condList = condList(ind);
        tfdat.powspctrm = single(tfdat.powspctrm(ind,:,:,:));
        data = log10(tfdat.powspctrm);

        % baseline
%         btoi = tfdat.time>-.51&tfdat.time<0; % time window for baseline
%         bdat = nanmean(data(:,:,:,btoi),4);
%         data = data - repmat(bdat,[1,1,1,size(data,4)]); % remove baseline
            
        
        %Run GLM
        load(fullfile(sessDir,'SuperWM.mat')); %load subject data
        
        Acc = abs(output_data(ind,6));
        SS = output_data(ind,2);
        
        chois = [53 54 59 56 57 58];
        
        freqs = [4:7];
        foi = find(tfdat.freq>=freqs(1)&tfdat.freq<=freqs(2));
        
        des_mat = [Acc];
        
        for chn = 1:length(chois);
            c = chois(chn);
            
            for f = 1:length(freqs)
                cf = freqs(f);
                
                for tpt = 1:size(data,4)
                    
                    for side = [1 2]
                        
                        for setsi = 1:7
                    
                            CondsTemp = Conds(:,setsi);
                            
                            trls = find(ismember(condList,CondsTemp(side,:)));
                            
                            toFit = squeeze(data(trls,c,cf,tpt));
                            des_matTemp = des_mat(trls,:);
                            
                            if ~isnan(toFit(1))
                                B(chn,f,tpt,side,setsi,:) = glmfit(des_matTemp,toFit);
                            else
                                B(chn,f,tpt,side,setsi,:) = [NaN NaN];
                            end
                        
                        end
                    
                    end
                    
                end
            end
        end
        
        AllSessBeta(sess,:,:,:,:,:,:) = B;
        
    end
    
    AllSubBeta(sub,:,:,:,:,:,:) = nanmean(AllSessBeta,1);


end




%%

left_chans = {'O1','PO7','PO3'};
right_chans =  {'O2','PO8','PO4'};

left_ind = find(ismember(tfdat.label,left_chans));
right_ind = find(ismember(tfdat.label,right_chans));

[~,Lind] = ismember(left_ind,chois);
[~,Rind] = ismember(right_ind,chois);

frange = [1 30];
foi = find(tfdat.freq>=frange(1)&tfdat.freq<=frange(2));

trange = [-1000 3500];
toi = find((tfdat.time*1000)>=trange(1)&(tfdat.time*1000)<=trange(2));

dosubs = Low_WM;



%%

figure,
% hold on

for ss = 1:7
    
    regr = 2;
    smf = 5;
    times = tfdat.time(toi);
    
    
    betaLc = smooth(squeeze(nanmean(nanmean(nanmean(AllSubBeta(dosubs,Lind,:,toi,1,ss,regr),1),2),3)),smf);
    betaLi = smooth(squeeze(nanmean(nanmean(nanmean(AllSubBeta(dosubs,Rind,:,toi,1,ss,regr),1),2),3)),smf);
    
    betaRc = smooth(squeeze(nanmean(nanmean(nanmean(AllSubBeta(dosubs,Rind,:,toi,2,ss,regr),1),2),3)),smf);
    betaRi = smooth(squeeze(nanmean(nanmean(nanmean(AllSubBeta(dosubs,Lind,:,toi,2,ss,regr),1),2),3)),smf);
    
    
    
    beta_c = (betaLc + betaRc)./2;
    beta_i = (betaLi + betaRi)./2;
    
    subplot(2,4,ss)

    col = [0 0 1];

    hold on
    plot(times, beta_c,         'LineWidth',3,...
        'Color',col);
    plot(times, beta_i,'--',    'LineWidth',3,...
        'Color',col);
    
    plot(times, beta_c-beta_i,  ':', 'LineWidth',1,...
        'Color',col);
    
    % axes = [-0.5 3 -0.0005 0.0005];
    % axes = [-0.5 3 -0.05 0.03];
    axes = [-.5 3 -0.01 0.01];
    axis(axes)
    
    hold off   
    
    
end




% hold off
