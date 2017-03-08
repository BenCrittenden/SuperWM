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

        % baseline
        btoi = tfdat.time>-.51&tfdat.time<0; % time window for baseline
        bdat = nanmean(data(:,:,:,btoi),4);
        data = data - repmat(bdat,[1,1,1,size(data,4)]); % remove baseline
        
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

frange = [8 12];
foi = find(tfdat.freq>=frange(1)&tfdat.freq<=frange(2));

trange = [-1000 3500];
toi = find((tfdat.time*1000)>=trange(1)&(tfdat.time*1000)<=trange(2));

dosubs = Low_WM;

%%

con1lc = squeeze(nanmean(AllSubDat(dosubs,1,right_ind,foi,toi),3));
con1li = squeeze(nanmean(AllSubDat(dosubs,1,left_ind,foi,toi),3));

con2lc = squeeze(nanmean(AllSubDat(dosubs,2,right_ind,foi,toi),3));
con2li = squeeze(nanmean(AllSubDat(dosubs,2,left_ind,foi,toi),3));

con3lc = squeeze(nanmean(AllSubDat(dosubs,3,right_ind,foi,toi),3));
con3li = squeeze(nanmean(AllSubDat(dosubs,3,left_ind,foi,toi),3));

con4lc = squeeze(nanmean(AllSubDat(dosubs,4,right_ind,foi,toi),3));
con4li = squeeze(nanmean(AllSubDat(dosubs,4,left_ind,foi,toi),3));

con5lc = squeeze(nanmean(AllSubDat(dosubs,5,right_ind,foi,toi),3));
con5li = squeeze(nanmean(AllSubDat(dosubs,5,left_ind,foi,toi),3));

con6lc = squeeze(nanmean(AllSubDat(dosubs,6,right_ind,foi,toi),3));
con6li = squeeze(nanmean(AllSubDat(dosubs,6,left_ind,foi,toi),3));

con7lc = squeeze(nanmean(AllSubDat(dosubs,7,right_ind,foi,toi),3));
con7li = squeeze(nanmean(AllSubDat(dosubs,7,left_ind,foi,toi),3));

con1rc = squeeze(nanmean(AllSubDat(dosubs,8,left_ind,foi,toi),3));
con1ri = squeeze(nanmean(AllSubDat(dosubs,8,right_ind,foi,toi),3));

con2rc = squeeze(nanmean(AllSubDat(dosubs,9,left_ind,foi,toi),3));
con2ri = squeeze(nanmean(AllSubDat(dosubs,9,right_ind,foi,toi),3));

con3rc = squeeze(nanmean(AllSubDat(dosubs,10,left_ind,foi,toi),3));
con3ri = squeeze(nanmean(AllSubDat(dosubs,10,right_ind,foi,toi),3));

con4rc = squeeze(nanmean(AllSubDat(dosubs,11,left_ind,foi,toi),3));
con4ri = squeeze(nanmean(AllSubDat(dosubs,11,right_ind,foi,toi),3));

con5rc = squeeze(nanmean(AllSubDat(dosubs,12,left_ind,foi,toi),3));
con5ri = squeeze(nanmean(AllSubDat(dosubs,12,right_ind,foi,toi),3));

con6rc = squeeze(nanmean(AllSubDat(dosubs,13,left_ind,foi,toi),3));
con6ri = squeeze(nanmean(AllSubDat(dosubs,13,right_ind,foi,toi),3));

con7rc = squeeze(nanmean(AllSubDat(dosubs,14,left_ind,foi,toi),3));
con7ri = squeeze(nanmean(AllSubDat(dosubs,14,right_ind,foi,toi),3));




%% For data plotting

[color1E,color2E,color4E,color1H,color2H,color4H,color1,color2,color4,colorE,colorH] = AttMemColors;
% 
% lat1r = ((con1rc - con1ri)./(con1rc + con1ri)).*100;
% lat1l = ((con1lc - con1li)./(con1lc + con1li)).*100;
% lat1 = (lat1r+lat1l)./2;
% 
% lat2r = ((con2rc - con2ri)./(con2rc + con2ri)).*100;
% lat2l = ((con2lc - con2li)./(con2lc + con2li)).*100;
% lat2 = (lat2r+lat2l)./2;
% 
% lat3r = ((con3rc - con3ri)./(con3rc + con3ri)).*100;
% lat3l = ((con3lc - con3li)./(con3lc + con3li)).*100;
% lat3 = (lat3r+lat3l)./2;
% 
% lat4r = ((con4rc - con4ri)./(con4rc + con4ri)).*100;
% lat4l = ((con4lc - con4li)./(con4lc + con4li)).*100;
% lat4 = (lat4r+lat4l)./2;
% 
% lat5r = ((con5rc - con5ri)./(con5rc + con5ri)).*100;
% lat5l = ((con5lc - con5li)./(con5lc + con5li)).*100;
% lat5 = (lat5r+lat5l)./2;
% 
% lat6r = ((con6rc - con6ri)./(con6rc + con6ri)).*100;
% lat6l = ((con6lc - con6li)./(con6lc + con6li)).*100;
% lat6 = (lat6r+lat6l)./2;
% 
% lat7r = ((con7rc - con7ri)./(con7rc + con7ri)).*100;
% lat7l = ((con7lc - con7li)./(con7lc + con7li)).*100;
% lat7 = (lat7r+lat7l)./2;

% lat1 = ((con1rc - con1ri) + (con1lc - con1li))./2;
% lat2 = ((con2rc - con2ri) + (con2lc - con2li))./2;
% lat3 = ((con3rc - con3ri) + (con3lc - con3li))./2;
% lat4 = ((con4rc - con4ri) + (con4lc - con4li))./2;
% lat5 = ((con5rc - con5ri) + (con5lc - con5li))./2;
% lat6 = ((con6rc - con6ri) + (con6lc - con6li))./2;
% lat7 = ((con7rc - con7ri) + (con7lc - con7li))./2;

lat1 = ((con1rc) + (con1lc))./2;
lat2 = ((con2rc) + (con2lc))./2;
lat3 = ((con3rc) + (con3lc))./2;
lat4 = ((con4rc) + (con4lc))./2;
lat5 = ((con5rc) + (con5lc))./2;
lat6 = ((con6rc) + (con6lc))./2;
lat7 = ((con7rc) + (con7lc))./2;

% lat1 = ((con1ri) + (con1li))./2;
% lat2 = ((con2ri) + (con2li))./2;
% lat3 = ((con3ri) + (con3li))./2;
% lat4 = ((con4ri) + (con4li))./2;
% lat5 = ((con5ri) + (con5li))./2;
% lat6 = ((con6ri) + (con6li))./2;
% lat7 = ((con7ri) + (con7li))./2;


% lat1 = con1lc;
% lat2 = con2lc;
% lat3 = con3lc;
% lat4 = con4lc;
% lat5 = con5lc;
% lat6 = con6lc;
% lat7 = con7lc;

% lat1 = con1li;
% lat2 = con2li;
% lat3 = con3li;
% lat4 = con4li;High
% lat5 = con5li;
% lat6 = con6li;
% lat7 = con7li;


%% Plotting set-size effect

clear error_mat
times = tfdat.time;
smf = 5;

clear line1 line2 line3 line4

line1 = smooth(squeeze(nanmean(nanmean(lat1,2),1)),smf);
line2 = smooth(squeeze(nanmean(nanmean(lat2,2),1)),smf);
line3 = smooth(squeeze(nanmean(nanmean(lat3,2),1)),smf);
line4 = smooth(squeeze(nanmean(nanmean(lat4,2),1)),smf);
line5 = smooth(squeeze(nanmean(nanmean(lat5,2),1)),smf);
line6 = smooth(squeeze(nanmean(nanmean(lat6,2),1)),smf);
line7 = smooth(squeeze(nanmean(nanmean(lat7,2),1)),smf);


% figure,
hold on

h(1) = plot(times, line1,'LineWidth',3,...
                         'Color',[1 0 0]);
h(2) = plot(times, line2,'b', 'LineWidth',3,...
                         'Color',[1 .15 .15]);
h(3) = plot(times, line3,'b', 'LineWidth',3,...
                         'Color',[1 .3 .3]);
h(4) = plot(times, line4,'b', 'LineWidth',3,...
                         'Color',[1 .45 .45]);
h(5) = plot(times, line5,'b', 'LineWidth',3,...
                         'Color',[1 .6 .6]);
h(6) = plot(times, line6,'b', 'LineWidth',3,...
                         'Color',[1 .75 .75]);
h(6) = plot(times, line7,'b', 'LineWidth',3,...
                         'Color',[1 .9 .9]); 
                     
% h(1) = plot(times, line1,'LineWidth',3,...
%                          'Color',[0 0 1]);
% h(2) = plot(times, line2,'b', 'LineWidth',3,...
%                          'Color',[.15 .15 1]);
% h(3) = plot(times, line3,'b', 'LineWidth',3,...
%                          'Color',[.3 .3 1]);
% h(4) = plot(times, line4,'b', 'LineWidth',3,...
%                          'Color',[.45 .45 1]);
% h(5) = plot(times, line5,'b', 'LineWidth',3,...
%                          'Color',[.6 .6 1]);
% h(6) = plot(times, line6,'b', 'LineWidth',3,...
%                          'Color',[.75 .75 1]);
% h(6) = plot(times, line7,'b', 'LineWidth',3,...
%                          'Color',[.9 .9 1]);
                     
                     
% h(1) = plot(times, line1,'LineWidth',3,...
%                          'Color',[0 0 0]);
% h(2) = plot(times, line2,'b', 'LineWidth',3,...
%                          'Color',[.15 .15 .15]);
% h(3) = plot(times, line3,'b', 'LineWidth',3,...
%                          'Color',[.3 .3 .3]);
% h(4) = plot(times, line4,'b', 'LineWidth',3,...
%                          'Color',[.45 .45 .45]);
% h(5) = plot(times, line5,'b', 'LineWidth',3,...
%                          'Color',[.6 .6 .6]);
% h(6) = plot(times, line6,'b', 'LineWidth',3,...
%                          'Color',[.75 .75 .75]);
% h(6) = plot(times, line7,'b', 'LineWidth',3,...
%                          'Color',[.9 .9 .9]);                     


% axes = [-.5 3 -.05 .1];                     
axes = [-.5 3 -.25 0.2];
% axes = [-.5 3 -.025 0.02];
% axes = [-.5 3 -.1 0.5];
% axes = [-.5 3 0 1]
axis(axes)

%add stimulus lines
plot([0 0],[axes(3) axes(4)],'k:','LineWidth',0.5);
plot([500 500],[axes(3) axes(4)],'k:','LineWidth',0.5);
plot([1000 1000],[axes(3) axes(4)],'k:','LineWidth',0.5);
plot([1500 1500],[axes(3) axes(4)],'k:','LineWidth',0.5);
plot([2000 2000],[axes(3) axes(4)],'k:','LineWidth',0.5);
plot([2500 2500],[axes(3) axes(4)],'k:','LineWidth',0.5);

plot([times(1) times(end)],[0 0],'k','LineWidth',0.5);

            
% set(leg1,'Location','NorthWest');
set(gca, 'XTick', 0:500:3000)
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



%% Plot power spectrum

figure,

for csub = 1:size(lat1,1)
    
    subplot(5,3,csub)
    
    imagesc(squeeze(lat1(csub,:,:)));
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



