%% pre-processing script

% Add and remove appropriate paths

% rmpath(genpath(['/home/benmc/matlab/spm8/external/fieldtrip']));
% rmpath(genpath('/opt/analysis/spm8/'))

clear
addpath('/home/benmc/matlab/');
dat_dir =  '/home/benmc/SuperWM/Data';  % input data

% Define subject files

sub_dir = 's06';
sub_details = {'st_2','st_3',};
cd(fullfile(dat_dir,sub_dir))

temp1 = load(fullfile(dat_dir,sub_dir,'raw_data',sub_details{1},'SuperWM.mat'));
temp2 = load(fullfile(dat_dir,sub_dir,'raw_data',sub_details{2},'SuperWM.mat'));
% temp3 = load(fullfile(dat_dir,sub_dir,'raw_data',sub_details{3},'SuperWM.mat'));
% temp4 = load(fullfile(dat_dir,sub_dir,'raw_data',sub_details{4},'SuperWM.mat'));


%Any lines need removing?
temp2.output_data = temp2.output_data(1:(end-1),:);

output_data = vertcat(temp1.output_data,temp2.output_data);%,temp3.output_data);
% SuperWM = vertcat(temp1.output_data,temp2.output_data,temp3.output_data,temp4.output_data);

if 0%size(SuperWM,1)~=(224*3)
    error('not 3 full blocks, check')
else
    save(fullfile(dat_dir,sub_dir,'SuperWM2.mat'),'output_data');
end


