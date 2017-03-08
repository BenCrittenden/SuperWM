%%

trash = input('Are you sure?','s');

clear all
close all

%addpath(genpath('/home/benmc/matlab/fieldtrip-20130507'));
addpath('/Users/benc/MATLAB/');

addpath('/Users/benc/MATLAB/shadedErrorBar');
addpath('/Users/benc/Experiments/SuperWM/Miscellaneous');

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


Conds = [11 21 31 41 51 61 71 ; 12 22 32 42 52 62 72];
    
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
        
        chois = [left_ind ;right_ind];
        
        des_mat = [Acc SS.^2 SS.^3]; 
        
        for chn = 1:length(chois);
            c = chois(chn);
            for tpt = 1:size(data,3)
                
                for side = 1:2
                trls = find(ismember(condList,Conds(side,:)));
                        
                toFit = squeeze(data(trls,c,tpt));
                des_matTemp = des_mat(trls,:);
                
                B(chn,tpt,side,:) = glmfit(des_matTemp,toFit); 
                
                end
            end            
        end
                
        AllSessBeta(sess,:,:,:,:) = B;
                
        
        
    end
    
    AllSubBeta(sub,:,:,:,:) = squeeze(nanmean(AllSessBeta,1));
    
end








%%

left_chans = {'O1','PO7','PO3'};
right_chans =  {'O2','PO8','PO4'};

left_ind = find(ismember(dat.label,left_chans));
right_ind = find(ismember(dat.label,right_chans));

[~,Lind] = ismember(left_ind,chois);
[~,Rind] = ismember(right_ind,chois);

trange = [-1000 3500];
toi = find((dat.time{1}*1000)>=trange(1)&(dat.time{1}*1000)<=trange(2));

dosubs =Low_WM;




%%

regr = 2;
smf = 25;
times = dat.time{1}(toi);


betaLc = gsmooth(squeeze(nanmean(nanmean(AllSubBeta(dosubs,Lind,toi,2,regr),1),2)),smf);
betaLi = gsmooth(squeeze(nanmean(nanmean(AllSubBeta(dosubs,Rind,toi,2,regr),1),2)),smf);

betaRc = gsmooth(squeeze(nanmean(nanmean(AllSubBeta(dosubs,Rind,toi,1,regr),1),2)),smf);
betaRi = gsmooth(squeeze(nanmean(nanmean(AllSubBeta(dosubs,Lind,toi,1,regr),1),2)),smf);



beta_c = (betaLc + betaRc)./2;
beta_i = (betaLi + betaRi)./2;

beta_R = (betaRc + betaRi)./2;
beta_

col = [0 0 1];

figure,

hold on

plot(times, beta_c,         'LineWidth',3,...
                            'Color',col);
plot(times, beta_i,'--',    'LineWidth',3,...
                            'Color',col);
                     
plot(times, beta_c-beta_i,  ':', 'LineWidth',1,...
                                  'Color',col);
                   
% axes = [-0.5 3 -0.02 .01];
% axes = [-0.5 3 -0.8 1];
axes = [-inf inf -inf inf];
axis(axes)

hold off


%%  

rgsr = 3;
smf = 25;

trange = [-100 3100];
toi = find((dat.time{1}*1000)>=trange(1)&(dat.time{1}*1000)<=trange(2));

times = dat.time{1};

allBeta = gsmooth(squeeze(nanmean(nanmean(nanmean(AllSubBeta(All_WM,:,:,:,rgsr),1),2),4)),smf);
HiBeta = gsmooth(squeeze(nanmean(nanmean(nanmean(AllSubBeta(High_WM,:,:,:,rgsr),1),2),4)),smf);
LoBeta = gsmooth(squeeze(nanmean(nanmean(nanmean(AllSubBeta(Low_WM,:,:,:,rgsr),1),2),4)),smf);



figure,
hold on

axes = [-inf inf -inf inf];
axes = [-0.5 3.5 -0.5 0.25];
% axes = [-0.5 3.5 -0.02 0.005];
% axes = [-0.5 3.5 -2.5 5];

%add stimulus lines
plot([0 0],[axes(3) axes(4)],...
                            'Color',[.4 .4 .4],...
                            'LineWidth',2);
plot([-.4 times(end)],[0 0],...
                            'Color',[.4 .4 .4],...
                            'LineWidth',2);

plot(times(toi), HiBeta(toi),'LineWidth',3,...
                         'Color',[1 0 0]);
plot(times(toi), LoBeta(toi),'LineWidth',3,...
                         'Color',[0 0 1]);
                   
% axes = [-0.5 3.5 -0.025 0.01];
axes = [-0.5 3 -.75 .5];
% axes = [-0.5 3 -3 9];
axis(axes)

% totest1 = squeeze(nanmean(AllSubBeta(High_WM,:,:,rgsr),2));
% totest2 = squeeze(nanmean(AllSubBeta(Low_WM(1:10),:,:,rgsr),2));
% 
% totest = totest1-totest2;
% 
% [datobs,datrnd] = cluster_test_helper(totest',nPerm);
% [h1,p,clusterinfo] = cluster_test(datobs,datrnd,Tail,clusThresh,SigThresh);
% plot_cluster = single(h1);
% plot_cluster(plot_cluster==0) = NaN;
% plot_cluster = plot_cluster*1;
% plot(times, plot_cluster,'*','Color',[0 1 0],'lineWidth',5);
% 
% plot([times(1) times(end)],[0 0],'k','LineWidth',0.5);

            
% set(leg1,'Location','NorthWest');
set(gca, 'XTick', 0:1:3)
set(gca, 'YTick', [0 .25])
whitebg('w')

Xlabel = xlabel('time, ms');
Ylabel = ylabel('beta');

set([gca Xlabel Ylabel],...
    'FontName',         'Helvetica',...
    'FontSize',         24,...
    'FontWeight',       'bold');

set(gca,...
    'Box',                  'off',...
    'TickDir',              'out',...
    'TickLength',           [0.005 0.005],...
    'XMinorTick',           'off',...
    'YMinorTick',           'off',...
    'YGrid',                'off',...
    'XGrid',                'off',...
    'XColor',               [0 0 0],...
    'YColor',               [0 0 0],...
    'LineWidth',             1);

hold off


%% Get the single sub data for stats and error bars

for s = 1:size(AllSubBeta,1)
    
    Sub_beta(s,:) = gsmooth(squeeze(mean(mean(AllSubBeta(s,:,:,:,rgsr),2),4)),smf);
    
    
end

Sb_error = errorBarCalc(Sub_beta,'se');
Sb_error_hi = errorBarCalc(Sub_beta(High_WM,:),'se');
Sb_error_lo = errorBarCalc(Sub_beta(Low_WM,:),'se');

hold on

shadedErrorBar(times(toi),HiBeta(toi),Sb_error_hi(toi),'r',1); 
shadedErrorBar(times(toi),LoBeta(toi),Sb_error_lo(toi),'b',1); 

hold off


for t = 1:size(Sb_error,2)
    
   at = Sub_beta(:,t);
   ht = Sub_beta(High_WM,t);
   lt = Sub_beta(Low_WM,t);
   
%    [H(t) p(t)] = ttest2(ht,lt);
   
   zt = zeros(length(ht),1);
   [H(t) p(t)] = ttest(ht,zt);
    
    
end

hold on

H(H==0) = NaN;
sig = H.*.4;
plot(times(toi),sig(toi),'LineWidth',3,...
                    'Color',[1 .7 .7],...
                    'LineStyle','-');
    
 hold off
 
 
 %% Lateralisation
 
regr = 3;
smf = 25;
times = dat.time{1}(toi);


betaLc = squeeze(nanmean(AllSubBeta(:,Lind,toi,2,regr),2));
betaLi = squeeze(nanmean(AllSubBeta(:,Rind,toi,2,regr),2));

betaRc = squeeze(nanmean(AllSubBeta(:,Rind,toi,1,regr),2));
betaRi = squeeze(nanmean(AllSubBeta(:,Lind,toi,1,regr),2));

Lat_L = betaLc-betaLi;
Lat_R = betaRc-betaRi;


for i = 1:size(AllSubBeta,1)
    
%     Lc(i,:) = gsmooth(betaLc(i,:),smf);
%     Li(i,:) = gsmooth(betaLi(i,:),smf);
%     Rc(i,:) = gsmooth(betaRc(i,:),smf);
%     Ri(i,:) = gsmooth(betaRi(i,:),smf);
    
    Ll(i,:) = gsmooth(Lat_L(i,:),smf);
    Lr(i,:) = gsmooth(Lat_R(i,:),smf);
    
    Lat(i,:) = gsmooth(mean([Lat_L(i,:); Lat_R(i,:)]),smf);

end

% 
% Ll_error = errorBarCalc(Ll,'se');
% Lr_error = errorBarCalc(Lr,'se');
% 
% Ll_error_hi = errorBarCalc(Ll(High_WM,:),'se');
% Lr_error_hi = errorBarCalc(Lr(High_WM,:),'se');
% 
% Ll_error_lo = errorBarCalc(Ll(Low_WM,:),'se');
% Lr_error_lo = errorBarCalc(Lr(Low_WM,:),'se');
% 
% Ll_mean = mean(Ll,1);
% Lr_mean = mean(Lr,1);
% 
% Ll_mean_hi = mean(Ll(High_WM,:),1);
% Lr_mean_hi = mean(Lr(High_WM,:),1);
% 
% Ll_mean_lo = mean(Ll(Low_WM,:),1);
% Lr_mean_lo = mean(Lr(Low_WM,:),1);
% 
% figure
% 
% hold on
% 
% shadedErrorBar(times,Lr_mean_hi,Lr_error_hi,'r',1); 
% shadedErrorBar(times,Lr_mean_lo,Lr_error_lo,'b',1); 
% 
% hold off

%

Lat_error = errorBarCalc(Lat,'se');
Lat_mean = mean(Lat,1);

Lat_error_hi = errorBarCalc(Lat,'se');
Lat_error_lo = errorBarCalc(Lat,'se');

Lat_mean_hi = mean(Lat(High_WM,:),1);
Lat_mean_lo = mean(Lat(Low_WM,:),1);


figure,

hold on

axes = [-inf inf -inf inf];
axes = [-0.5 3.5 -0.005 0.005];


%add stimulus lines
plot([0 0],[axes(3) axes(4)],...
                            'Color',[.4 .4 .4],...
                            'LineWidth',2);
plot([-.4 times(end)],[0 0],...
                            'Color',[.4 .4 .4],...
                            'LineWidth',2);

plot(times, Lat_mean_hi,'LineWidth',3,...
                         'Color',[1 0 0]);
                     
plot(times, Lat_mean_lo,'LineWidth',3,...
                         'Color',[0 0 1]);
                     
% plot(times, Lat_mean,'LineWidth',3,...
%                          'Color',[.5 .5 .5]);                     
                     
                   
axes = [-0.5 3.5 -.75 .5];



            
% set(leg1,'Location','NorthWest');
set(gca, 'XTick', 0:1:3)
set(gca, 'YTick', [0 .0025])
whitebg('w')

Xlabel = xlabel('time, ms');
Ylabel = ylabel('beta');

set([gca Xlabel Ylabel],...
    'FontName',         'Helvetica',...
    'FontSize',         24,...
    'FontWeight',       'bold');

set(gca,...
    'Box',                  'off',...
    'TickDir',              'out',...
    'TickLength',           [0.005 0.005],...
    'XMinorTick',           'off',...
    'YMinorTick',           'off',...
    'YGrid',                'off',...
    'XGrid',                'off',...
    'XColor',               [0 0 0],...
    'YColor',               [0 0 0],...
    'LineWidth',             1);


shadedErrorBar(times,Lat_mean_hi,Lat_error_hi,'r',1); 
shadedErrorBar(times,Lat_mean_lo,Lat_error_lo,'b',1); 
% shadedErrorBar(times,Lat_mean,Lat_error,'k',1); 


hold off

clear H p

for t = 1:size(Lat,2)
    
   ht = Lat(High_WM,t);
   lt = Lat(Low_WM,t);
   bt = Lat(:,t);
   
   [H(t) p(t)] = ttest2(ht,lt);
   
%    zt = zeros(length(ht),1);
%    [H(t) p(t)] = ttest(ht,zt);
    
    
end

hold on

H(H==0) = NaN;
sig = H.*.005;
plot(times,sig,'LineWidth',3,...
                    'Color',[.5 .5 .5],...
                    'LineStyle','-');
    
hold off

