%%

trash = input('Are you sure?','s');

clearvars
close all

%addpath(genpath('/home/benmc/matlab/fieldtrip-20130507'));
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
        btoi = find(dat.time{1}<0 & dat.time{1}>-0.1); %100ms represents 1 alpha cycle
        bdat = mean(data(:,:,btoi),3);
        data = data - repmat(bdat,[1,1,size(data,3)]); % remove baseline
        
        nCond = size(Conds,2);
        for c=1:nCond
            selectTrials = find(condList==Conds(c));
            AllSessDat(sess,c,:,:) = squeeze(nanmean(data(selectTrials,:,:),1));
        end
    end
    
    AllSubDat(sub,:,:,:) = squeeze(nanmean((AllSessDat),1));
    
end



%%     

left_chans = {'O1','PO7','PO3'};
right_chans =  {'O2','PO8','PO4'};

left_ind = find(ismember(dat.label,left_chans));
right_ind = find(ismember(dat.label,right_chans));

trange = [-100 3100];
toi = find((dat.time{1}*1000)>=trange(1)&(dat.time{1}*1000)<=trange(2));

ss = 1;

%% 

Hconrc = squeeze(nanmean(AllSubDat(High_WM,ss,right_ind,toi),3));
Hconli = squeeze(nanmean(AllSubDat(High_WM,ss,left_ind,toi),3));

Hconri = squeeze(nanmean(AllSubDat(High_WM,ss+7,right_ind,toi),3));
Hconlc = squeeze(nanmean(AllSubDat(High_WM,ss+7,left_ind,toi),3));

Lconrc = squeeze(nanmean(AllSubDat(Low_WM,ss,right_ind,toi),3));
Lconli = squeeze(nanmean(AllSubDat(Low_WM,ss,left_ind,toi),3));

Lconri = squeeze(nanmean(AllSubDat(Low_WM,ss+7,right_ind,toi),3));
Lconlc = squeeze(nanmean(AllSubDat(Low_WM,ss+7,left_ind,toi),3));

Aconrc = squeeze(nanmean(AllSubDat(:,ss,right_ind,toi),3));
Aconli = squeeze(nanmean(AllSubDat(:,ss,left_ind,toi),3));

Aconri = squeeze(nanmean(AllSubDat(:,ss+7,right_ind,toi),3));
Aconlc = squeeze(nanmean(AllSubDat(:,ss+7,left_ind,toi),3));

%%

% Hlatr = ((Hconrc - Hconri)./(Hconrc + Hconri)).*100;
% Hlatl = ((Hconlc - Hconli)./(Hconlc + Hconli)).*100;
% Hlat = (Hlatr+Hlatl)./2;
% 
% Llatr = ((Lconrc - Lconri)./(Lconrc + Lconri)).*100;
% Llatl = ((Lconlc - Lconli)./(Lconlc + Lconli)).*100;
% Llat = (Llatr+Llatl)./2;

Hlat = ((Hconrc - Hconri) + (Hconlc - Hconli))./2;
Llat = ((Lconrc - Lconri) + (Lconlc - Lconli))./2;

% Hlat = ((Hconrc + Hconlc)./2) - ((Hconri + Hconli)./2);
% Llat = ((Lconrc + Lconlc)./2) - ((Lconri + Lconli)./2);
% 
Hc = (Hconrc + Hconlc)./2;
Lc = (Lconrc + Lconlc)./2;

Hi = (Hconri + Hconli)./2;
Li = (Lconri + Lconli)./2;

% Hc = Hconlc;
% Lc = Lconlc;

% Hi = Hconli;
% Li = Hconri;
% 
Hlat = Hi;
Llat = Li;

Alat = ((Aconrc + Aconlc)./2) - ((Aconri + Aconli)./2);

%% Plotting set-size effect

clear error_mat
times = dat.time{1};
times = times(toi);
smf = 20;

Hlat_error = errorBarCalc(Hlat,'se');
Llat_error = errorBarCalc(Llat,'se');
Alat_error = errorBarCalc(Alat,'se');

lineH = gsmooth(squeeze(nanmean(Hlat,1)),smf);
lineL = gsmooth(squeeze(nanmean(Llat,1)),smf);
lineA = gsmooth(squeeze(nanmean(Alat,1)),smf);

H_error = gsmooth(Hlat_error,smf);
L_error = gsmooth(Llat_error,smf);
A_error = gsmooth(Alat_error,smf);

figure,
hold on

% axes = [-inf inf -inf inf];                     
axes = [-0.5 3 -0.5 3.5];

%add stimulus lines
plot([0 0],[axes(3) axes(4)],...
                            'Color',[.4 .4 .4],...
                            'LineWidth',2);
plot([-.4 times(end)],[0 0],...
                            'Color',[.4 .4 .4],...
                            'LineWidth',2);

% h(1) = plot(times, lineH,'LineWidth',3,...
%                          'Color',[1 0 0]);
%                      
% h(2) = plot(times, lineL,'LineWidth',3,...
%                          'Color',[0 0 1]);
                     
h(3) = plot(times, lineA,'LineWidth',3,...
                         'Color',[0 0 0]);
                                          
% shadedErrorBar(times,lineH,H_error,'r',1);
% shadedErrorBar(times,lineL,L_error,'b',1);
shadedErrorBar(times,lineA,A_error,'k',1);                    



% axis(axes)

            
% set(leg1,'Location','NorthWest');
set(gca, 'XTick', 0:1:3)
set(gca, 'YTick', [0 3])
whitebg('w')

Xlabel = xlabel('time, ms');
Ylabel = ylabel('voltage');

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

for tp = 1:size(Hlat,2)
    
    
%     Htest = Hlat(:,tp);
%     Ltest = Llat(:,tp);
    
    for r=1:size(Hlat,1)
        x(r,:) = gsmooth(Hlat(r,:),smf);
    end
    
    for r=1:size(Llat,1)
        y(r,:) = gsmooth(Llat(r,:),smf);
    end
    
    for r=1:size(Alat,1)
        z(r,:) = gsmooth(Alat(r,:),smf);
    end

    Htest = x(:,tp);
    Ltest = y(:,tp);
    Atest = z(:,tp);
    
    Ztest = zeros(size(z,1),1);
    
    
%     [H(tp) pval(tp)] = ttest2(Htest,Ltest,'tail','right');
    [H(tp) pval(tp)] = ttest(Atest,Ztest,'tail','right');
    
end

H(H==0)=NaN;

hold on

plot(times,-0.5.*H,...
                   'Color',[.4 .4 .4],...
                   'LineWidth',4);

hold off
               
% figure,imagesc(H);



%% Find the CDA peak

trange = [0 1000];
toi = find((times*1000)>=trange(1)&(times*1000)<=trange(2));
Max_H = max((x(:,toi))');
Max_L = max((y(:,toi))');

figure,bar([Max_H NaN Max_L]);

[H p] = ttest2(Max_H,Max_L,'tail','right')



