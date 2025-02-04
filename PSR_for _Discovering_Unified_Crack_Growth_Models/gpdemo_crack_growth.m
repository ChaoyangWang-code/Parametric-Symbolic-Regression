clc;
clear all

%% if HPC, use:
%folderPath = 'workdir';     %example:    Windows:'C:\Users\Desktop\test250106\'   linux:'/gs/home/gp_test/test1';
%addpath(genpath(folderPath));
% if isempty(gcp('nocreate'))
%     parpool(48);
% end

gp = rungp(@gpdemo_crack_growth_config);
paretoid(gp);

quit

