
clc
clear all

color = [
    1, 0, 0;        
    0, 1, 0;        
    0, 0, 1;        
    1, 0.843, 0;    
    0.580, 0, 0.827; 
    1, 0.647, 0;    
    0.529, 0.808, 0.922; 
    0.133, 0.545, 0.133; 
    1, 0.078, 0.576; 
    0, 0.808, 0.820  
];

%figuresize=[100, 100, 1200, 1000]; %[left bottom width height]
for gen_count=60
    figuresize=[100, 100, 1000, 800]; %[left bottom width height]
    %color = [1 0.38 0;1 0 0; 0 0 1;0.63 0.13 0.94;0.5 1 0;1 0.5 0.31; 0.53 0.15 0.34; 1 0 1;1 1 0; 0.45 0.29 0.07;0 0 1;0.03 0.18 0.33];  
    colors=[color;color;color;color];
    Path='workdir';      %Select the result path,example:    Windows:'C:\Users\Desktop\test250106\'  
    show_figure='on';  
    str=num2str(gen_count);
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);
    
    
    foldername_gen=[Path, str];
    %show_results_gen_i(figuresize,gp,foldername_gen,str)
    best_pop=MULTI_R_show_results_of_all_pareto_in_gen_i(colors,show_figure,figuresize,data_geni{gen_count}.gp,foldername_gen,str);
    clc
end



function [ best_pop ] = MULTI_R_show_results_of_all_pareto_in_gen_i(colors,show_figure,figuresize,gp,foldername_gen,str)

tour_ind = 1:gp.runcontrol.pop_size;
tour_ind=tour_ind';
tour_fitness = gp.fitness.values;
tour_comp  = gp.fitness.complexity;
%perform fast pareto sort on tournament members
if gp.fitness.minimisation
    mo = [tour_fitness tour_comp];
else
    mo = [ (-1 *tour_fitness) tour_comp];
end
%remove any 'infs' from consideration and recall if necessary
infs= isinf(mo(:,1));
mo(infs,:)=[];
if isempty(mo)
    error('All inf');
end
tour_ind(infs,:)=[];
rank = ndfsort_rank1(mo);
rank1_tour_ind_temp = tour_ind(rank == 1);
num_rank1 = numel(rank1_tour_ind_temp);
pareto_fit=tour_fitness(rank1_tour_ind_temp,1);
pareto_comp=tour_comp(rank1_tour_ind_temp,1);
%%
p_all=[pareto_fit,pareto_comp];
pu=unique(p_all,'rows','stable');             
pu = sortrows(pu, 1);          
pu_num=size(pu,1);
rank1_tour_ind=[];
for ui=1:pu_num
    row_to_find = pu(ui,:);
    matches = ismember(p_all, row_to_find, 'rows');
    [row, ~] = find(matches);
    rank1_tour_ind(ui,1)=rank1_tour_ind_temp(row(1),1);
end



%% 
returns=gp.fitness.returnvalues;

pareto_i=1;

str_pareto=num2str(pareto_i);
pvals=rank1_tour_ind(pareto_i);
pop_now=returns{pvals};
parameter_now=pop_now{1};
numGenes=pop_now{2};
num_const=pop_now{3};    
F_eq=pop_now{5};
for i=1:numGenes+2
    syms(['f',num2str(i)]);
end
for i=1:num_const
    syms(['c',num2str(i)]);
end

F_eq_show_temp=F_eq;
syms delta_K Kmax k1 k2
eq=['log10(f1)'];

for i=1:numGenes
    x1_str=['(','delta_K/','k1',')'];
    x2_str=['(','Kmax/','k2',')'];
    F_eq_show_temp{i} = strrep(F_eq_show_temp{i},'x1',x1_str);
    F_eq_show_temp{i} = strrep(F_eq_show_temp{i},'x2',x2_str);
    %eq_part{i}=['f',num2str(i+1),'*',F_eq_show_temp{i}];
    eq_part{i}=['f',num2str(i+1),'*','log10(',F_eq_show_temp{i},')',];
    eq=[eq,'+',eq_part{i}];
end
eq=[eq,'+','f',num2str(numGenes+2),'*','log10(delta_K)'];
eq10=['10^(',eq,')'];

pat2 = 'c(\d+)';
Const_pair_now=parameter_now(1,4:3+num_const);   
%         eq = regexprep(eq,pat2,'Const_pair_now($1)');
%         eq10 = regexprep(eq10,pat2,'Const_pair_now($1)');
%         F_eq_show=char(simplify(eval(eq)));
%         F_eq_show10=char(simplify(eval(eq10)));
F_eq_show=eq;
F_eq_show10=eq10;

materials_N=size(gp.multidata,2);                 
data_N=size(gp.all_data_without_collection,2);    
data_all=gp.all_data_without_collection;
data_all_collection=gp.multidata;
foldername_pareto_i_fig=[foldername_gen,'\','pareto_',num2str(pareto_i),'_fig'];
if ~exist(foldername_pareto_i_fig, 'dir')
    mkdir(foldername_pareto_i_fig);
end        
foldername_pareto_i_png=[foldername_gen,'\','pareto_',num2str(pareto_i),'_png'];
if ~exist(foldername_pareto_i_png, 'dir')
    mkdir(foldername_pareto_i_png);
end  

%% 
theta_temp=parameter_now(:,4+num_const:4+num_const+numGenes+1); 
theta = [theta_temp(:, 1), theta_temp(:, end), theta_temp(:, 2:end-1)];
% 
num_coefficients = size(theta, 2);
% 
coeff_labels = arrayfun(@(x) ['f' num2str(x)], 1:num_coefficients, 'UniformOutput', false);
h2 = figure('Visible', show_figure);
BOX = boxplot(theta, 'Labels', coeff_labels);
set(BOX, 'LineWidth', 2); 

set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
%xlabel('Parameters', 'FontName', 'Times New Roman', 'FontSize', 16);
ylabel('Value', 'FontName', 'Times New Roman', 'FontSize', 16);
%title(['Distribution of Coefficients for all labeled data'], 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
% 
xticks(1:num_coefficients);
xticklabels(coeff_labels);
% 
%set(gca, 'Box', 'on', 'LineWidth', 2);
grid on;
ylim([-12, 10]);  


set(gca, 'XTick', []);

set(gca, 'GridLineStyle', '-'); 
set(gca, 'GridColor', [0.6 0.6 0.6]); 



%% 
load('NASGRO_fit_nondimensional.mat');
NASGRO_fit=NASGRO_fit_nondimensional;
% 
num_coefficients = size(NASGRO_fit, 2);
% 
coeff_labels = arrayfun(@(x) ['f' num2str(x)], 1:num_coefficients, 'UniformOutput', false);
h2 = figure('Visible', show_figure);
%BOX = boxplot(NASGRO_fit);
BOX = boxplot(NASGRO_fit, 'Labels', coeff_labels);
set(BOX, 'LineWidth', 2); 

set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
%xlabel('Parameters', 'FontName', 'Times New Roman', 'FontSize', 16);
ylabel('Value', 'FontName', 'Times New Roman', 'FontSize', 16);
%title(['Distribution of Coefficients for all labeled data'], 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');

%set(gca, 'Box', 'on', 'LineWidth', 2);
grid on;
ylim([-12, 10]);  

set(gca, 'XTick', []);

set(gca, 'GridLineStyle', '-'); 
set(gca, 'GridColor', [0.6 0.6 0.6]); 



end

