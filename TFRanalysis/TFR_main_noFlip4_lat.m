trash = input('Are you sure?','s');

clear all

addpath('/home/benmc/matlab/');
eeglabpath = '/home/benmc/matlab/eeglab11_0_5_4b';
rmpath(genpath(eeglabpath));
addpath('/Users/benc/MATLAB/');
addpath('/Users/benc/MATLAB/shadedErrorBar');
addpath('/Users/benc/Experiments/SuperWM/Miscellaneous')

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
        %         data = tfdat.powspctrm;
        
        % baseline
%         btoi = tfdat.time>-.51&tfdat.time<0; % time window for baseline
%         bdat = nanmean(data(:,:,:,btoi),4);
%         data = data - repmat(bdat,[1,1,1,size(data,4)]); % remove baseline
        
        for con = 1:numel(Conds)
            selectTrials = find(condList==Conds(con));
            AllSessDat(sess,con,:,:,:) = squeeze(nanmean(data(selectTrials,:,:,:),1));
        end
        
    end
    
    AllSubDat(sub,:,:,:,:) = nanmean(AllSessDat,1);
    
    
end

%AllSubDat is a N subs x 6 conditions x 60 channels x 29 frequencies x 80 timepoints matrix

%%

% save('TFR_dat.mat','AllSubDat','-v7.3')

%%

% load('/home/benmc/SuperWM/Pipeline/TFRanalysis/TFR_High_dat.mat')
% load('/home/benmc/SuperWM/Data/s04_1/TFR.mat');

%%

left_chans = {'O1','PO7','PO3'};
right_chans =  {'O2','PO8','PO4'};


left_ind = find(ismember(tfdat.label,left_chans));
right_ind = find(ismember(tfdat.label,right_chans));

frange = [8 12]; %8 12
foi = find(tfdat.freq>=frange(1)&tfdat.freq<=frange(2));

trange = [-1000 3500];
toi = find((tfdat.time*1000)>=trange(1)&(tfdat.time*1000)<=trange(2));

%%
    
    trange = [-200 3000];
    toi = find((tfdat.time*1000)>=trange(1)&(tfdat.time*1000)<=trange(2));
    
    ss1 = 1:7;
    ss2 = ss1+7;
%     ss1 = 3;
%     ss2 = ss1+7;
    %%
    
    conlc = squeeze(nanmean(nanmean(AllSubDat(:,ss2,right_ind,foi,toi),2),3));
    conli = squeeze(nanmean(nanmean(AllSubDat(:,ss2,left_ind,foi,toi),2),3));
    
    conrc = squeeze(nanmean(nanmean(AllSubDat(:,ss1,left_ind,foi,toi),2),3));
    conri = squeeze(nanmean(nanmean(AllSubDat(:,ss1,right_ind,foi,toi),2),3));
    
    
    %% For data plotting
    
    lat = ((conrc - conri) + (conlc - conli))./2;
    
    latc = (conrc + conlc)./2;
    lati = (conri + conli)./2;
    
    all = (conrc + conlc + conri + conli)./4;
    
    %% Average across frequencies
    
    lat = squeeze(nanmean(lat,2));
    latc = squeeze(nanmean(latc,2));
    lati = squeeze(nanmean(lati,2));
    all = squeeze(nanmean(all,2));
    
    %% Plotting set-size effect
    
    clear error_mat
    times = tfdat.time(toi);
    smf = 4;
    
%     trange = [-500 3000];
%     toi = find((tfdat.time*1000)>=trange(1)&(tfdat.time*1000)<=trange(2));
    
%     lat = all;
    
    H_error = gsmooth(errorBarCalc(lat(High_WM,:),'se'),smf);
    L_error = gsmooth(errorBarCalc(lat(Low_WM,:),'se'),smf);
    A_error = gsmooth(errorBarCalc(lat,'se'),smf);
    
    clear line1 line2 line3
    
    lineH = gsmooth(squeeze(nanmean(lat(High_WM,:),1)),smf);
    lineL = gsmooth(squeeze(nanmean(lat(Low_WM,:),1)),smf);
    lineA = gsmooth(squeeze(nanmean(lat(:,:),1)),smf);
    
    
    
    figure,
    hold on
    
    %add stimulus lines
    axes = [-.25 3 -0.01 0.01];
    % axes = [-.25 3 -10 5];
    
    plot([0 0],[axes(3) axes(4)],'LineWidth',2,...
        'Color',[.5 .5 .5]);
    plot([axes(1) axes(2)],[0 0],'LineWidth',2,...
        'Color',[.5 .5 .5]);
    
    h(1) = plot(times, lineH,'LineWidth',3,...
        'Color',[1 0 0]);
    h(2) = plot(times, lineL, 'LineWidth',3,...
        'Color',[0 0 1]);
    h(3) = plot(times, lineA,'LineWidth',3,...
                             'Color',[.3 .3 .3]);
    
    shadedErrorBar(times,lineH,H_error,'r',1);
    shadedErrorBar(times,lineL,L_error,'b',1);
    shadedErrorBar(times,lineA,A_error,'k',1);
    
    
%     axes = [-1 3.5 -1 1.5];
%     axes = [-1 3.5 -inf inf];
    % axes = [-1 3.5 -25 10];
    axes = [-1 3.5 -0.04 0.02];
    axis(axes)
    
    
    % set(leg1,'Location','NorthWest');
    set(gca, 'XTick', 0:1:3)
    set(gca, 'YTick', [0 0.01])
    % set(gca, 'YTick', [0 10])
    whitebg('w')
    
    Xlabel = xlabel('time, ms');
    Ylabel = ylabel('power');
    
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
    
    %%
    
    clear H p
    
    for t = 1:size(lat,2)
        
        for s = 1:size(lat,1)
            slat(s,:) = gsmooth(lat(s,:),smf);
        end
        
        Ht = slat(High_WM,t);
        Lt = slat(Low_WM,t);
        At = lat(:,t);
        
        [H(t) p(t)] = ttest2(Ht,Lt);
        
        %     zt = zeros(size(At,1),1);
        %     [H(t) p(t)] = ttest(Ht,Lt);
        
        
    end
    
    hold on
    H(H==0) = NaN;
    H = H .* -0.6;
    % H = H .* 5;
    plot(times,H, 'LineWidth',3,...
        'Color',[.3 .3 .3]);
    
    %% find the minimum desync
    
    trange = [50 750];
    toi2 = find((times(toi)*1000)>=trange(1)&(times(toi)*1000)<=trange(2));
    
    Hi = slat(High_WM,toi2);
    Lo = slat(Low_WM,toi2);
    
    Hi_min(ss,:) = min(Hi');
    Lo_min(ss,:) = min(Lo');
    
    Hi_max(ss,:) = max(Hi');
    Lo_max(ss,:) = max(Lo');
    

    

%%
Hi = mean((Hi_max-Hi_min)');
Lo = mean((Lo_max-Lo_min)');

% Hi = mean(Hi_max');
% Lo = mean(Lo_max');

figure,
hold on

bar([Hi' Lo'],1.5);
legend({'High WM','Low WM'});
ylabel('peak alpha desync')
xlabel('Set Size')

ax = gca;
ax.YTick = [0 0.5];
ax.XTick = [1:7];

% fA = vertcat([ones(size(Hi_max,2),1) Hi_max'],[2*ones(size(Lo_max,2),1) Lo_max']);
fA = vertcat([ones(size(Hi_max,2),1) (Hi_max-Hi_min)'],[2*ones(size(Lo_max,2),1) (Lo_max-Lo_min)']);

figure,
hold on
plot(Lo_min,'b')
plot(Hi_min,'r')
hold off


%% find power at trial onset

trange = [-100 0];
toi3 = find((tfdat.time*1000)>=trange(1)&(tfdat.time*1000)<=trange(2));

allb = squeeze(nanmean(nanmean(nanmean(nanmean(AllSubDat(:,:,:,foi,toi3),2),3),4),5));

Hi = mean(allb(High_WM));
Lo = mean(allb(Low_WM));

figure,
bar([allb(High_WM)' NaN allb(Low_WM)']);

figure,
hold on
bar([mean(allb(High_WM));NaN],'r')
bar([NaN ; mean(allb(Low_WM))],'b');
axis([-inf inf 0 1])
legend({'High WM','Low WM'})
hold off

% load the parameter estimates from the VP_A model
load('/Users/benc/Experiments/SuperWM/old/code/code/model_params_gr.mat');

figure,
scatter(allb,params(:,1));

figure
hold on
scatter(allb(High_WM),params(High_WM,4),'r');
scatter(allb(Low_WM),params(Low_WM,4),'b');
axis([-inf inf -inf inf])
% axis([-0.5 1.5 -100 500])
hold off

[rho pval] = corr(allb,params(:,1));
[rho_hi pval_hi] = corr(allb(High_WM),params(High_WM,1));
[rho_lo pval_lo] = corr(allb(Low_WM),params(Low_WM,1));



[H p] = ttest2(allb(High_WM),allb(Low_WM));

%%
load('/Users/benc/Experiments/SuperWM/old/code/code/model_params_gr.mat');

%params has the model fit parameters. 1 is the interesting one

figure,

for sp = 1:7
    subplot(3,3,sp)
    hold on
    axis([0.5 3 -1.5 .5])
    scatter(saveparamHi(1,:),Hi_min(sp,:),'r');
    scatter(saveparamLo(1,:),Lo_min(sp,:),'b');
    
    [Hc_Hi(sp) pc_Hi(sp)] = corr(saveparamHi(1,:)',Hi_min(sp,:)');
    
    coeffs = polyfit(saveparamHi(1,:)', Hi_min(sp,:)', 1);
    fittedX = linspace(1,3,10);
    fittedY = polyval(coeffs, fittedX);
    plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
 
    [Hc_Lo(sp) pc_Lo(sp)] = corr(saveparamLo(1,:)',Lo_min(sp,:)');
    [Hc(sp) pc(sp)] = corr([saveparamHi(1,:) saveparamLo(1,:)]',...
                            [Hi_min(sp,:) Lo_min(sp,:)]');
                        
    coeffs = polyfit(saveparamLo(1,:)', Lo_min(sp,:)', 1);
    fittedX = linspace(1,3,10);
    fittedY = polyval(coeffs, fittedX);
    plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
    
    
    hold off
end

figure,
hold on
scatter(saveparamHi(1,:),mean(Hi_min,1),'r');
scatter(saveparamLo(1,:),mean(Lo_min,1),'b');
axis([0 3 -1 1])
hold off

[Hcm_Hi pcm_Hi] = corr(saveparamHi(1,:)',mean(Hi_min,1)');
[Hcm_Lo pcm_Lo] = corr(saveparamLo(1,:)',mean(Lo_min,1)');
[Hcm pcm] = corr([saveparamHi(1,:) saveparamLo(1,:)]',...
                            mean([Hi_min(:,:) Lo_min(:,:)])');


%% Plot power spectrum

figure,

for csub = 1:size(lat1,1)
    
    subplot(5,3,csub)
    
    imagesc(squeeze(lat(csub,:,:)));
    caxis([0 2]);
    axis xy
    
end




%%

load('/home/benmc/SuperWM/LoWMscores.mat')

trange = [0 3000];
toi = find((tfdat.time*1000)>=trange(1)&(tfdat.time*1000)<=trange(2));



%%

figure,
[a b] = sort(LoWMscores(:,1));
for p = 1:11
    
    c = b(p);
    subplot(3,5,p)
    
    plot(squeeze(nanmean(lat1(c,:,toi),2))')
    axis([0 length(toi) -0.5 2])
    
end




