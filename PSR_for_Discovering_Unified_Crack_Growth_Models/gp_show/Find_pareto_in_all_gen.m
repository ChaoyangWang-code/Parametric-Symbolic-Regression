

clc
clear all
Path='workdir';     %Select the result path,example:    Windows:'C:\Users\Desktop\test250106\'   linux:'/gs/home/gp_test/test1/';
parfor gen_count=1:35
    
    str=num2str(gen_count);
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);
       
    foldername_gen=[Path, str];
    %show_results_gen_i(figuresize,gp,foldername_gen,str)
    pareto_all_gen{gen_count,1}=Find_pareto_pop_in_gen_i(data_geni{gen_count}.gp);
end


if ~exist(Path, 'dir')
    mkdir(Path);  
end


filename = fullfile(Path, 'pareto_all_gen.mat');

save(filename, 'pareto_all_gen');