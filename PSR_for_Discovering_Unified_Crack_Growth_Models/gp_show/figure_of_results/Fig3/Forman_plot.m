
clc
clear all

colors = [
    0 0.4470 0.7410;    
    0.8500 0.3250 0.0980; 
    0.9290 0.6940 0.1250; 
    0.4940 0.1840 0.5560; 
    0.4660 0.6740 0.1880; 
    0.3010 0.7450 0.9330; 
    0.6350 0.0780 0.1840; 
    0.25 0.25 0.25;       
    1 0 1;                
    0 1 1;                
    0 1 0;                
    0.5 0 0;              
    0 0 0.5;              
];

%figuresize=[100, 100, 1200, 1000]; %[left bottom width height]
for gen_count=60
    figuresize=[100, 100, 1000, 800]; %[left bottom width height]
    %color = [1 0.38 0;1 0 0; 0 0 1;0.63 0.13 0.94;0.5 1 0;1 0.5 0.31; 0.53 0.15 0.34; 1 0 1;1 1 0; 0.45 0.29 0.07;0 0 1;0.03 0.18 0.33];  
    %colors=[color;color;color;color];
    Path='workdir';      %Select the result path,example:    Windows:'C:\Users\Desktop\test250106\'  
    
    show_figure='on';  
    str=num2str(gen_count);
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);
    
    
    foldername_gen=[Path, str];
    %show_results_gen_i(figuresize,gp,foldername_gen,str)
    MULTI_R_show_results_of_all_pareto_in_gen_i(colors,show_figure,figuresize,data_geni{gen_count}.gp,foldername_gen,str);
    clc
end

function [  ] = MULTI_R_show_results_of_all_pareto_in_gen_i(colors,show_figure,figuresize,gp,foldername_gen,str)

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
F_eq_show=eq;
F_eq_show10=eq10;

materials_N=size(gp.multidata,2);                  
data_N=size(gp.all_data_without_collection,2);     
data_all=gp.all_data_without_collection;
data_all_collection=gp.multidata;


%% 
N_PRE=100;

GP_results=cell(2,data_N);
Experiment_results=cell(2,data_N);

for i=1:data_N
    
    material_index=data_all{2,i};
    parameter_K=parameter_now(material_index,2:3);
    delta_kth=parameter_K(1);
    kc=parameter_K(2);
    Const_pair_now=parameter_now(material_index,4:3+num_const);
    theta=parameter_now(material_index,4+num_const:4+num_const+numGenes+1);           
    %
    D_i=data_all{1,i};

    delta_K_pre=linspace(min(data_all{1,i}(:,1)), max(data_all{1,i}(:,1)),N_PRE);
    Kmax_pre=linspace(min(data_all{1,i}(:,4)), max(data_all{1,i}(:,4)),N_PRE);
    delta_K_pre=delta_K_pre';
    Kmax_pre=Kmax_pre';
    xtrain=[delta_K_pre./delta_kth,Kmax_pre./kc];
    F_eq=pop_now{5};
    pat1 = 'x(\d+)';
    pat2 = 'c(\d+)';
    F_eq = regexprep(F_eq,pat2,'Const_pair_now($1)');
    F_eq = regexprep(F_eq,pat1,'xtrain(:,$1)');
    geneOutputs = ones(N_PRE,numGenes+2);
    for j = 1:numGenes
        ind = j + 1;
        gene_temp=eval_exe(F_eq{j},xtrain,Const_pair_now);
        geneOutputs(:,ind)=lg(gene_temp);
    end
    geneOutputs(:,numGenes+2)=lg(delta_K_pre);               
    ypredtrain = geneOutputs * theta';

    % h1 = figure('Visible', show_figure); 
    % plot(lg(D_i(:,1)), lg(D_i(:,2)),'+','Color', colors(1, :),'MarkerSize', 10);
    % hold on; grid on
    % plot(lg(delta_K_pre),ypredtrain,'-','Color', colors(3, :),'LineWidth',1.5)  
    % xlabel('$\Delta K\ (MPa\sqrt{m})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14,'FontWeight', 'bold');
    % ylabel('$da/dN\ (mm/cycle)$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14,'FontWeight', 'bold');
    % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    GP_results{1,i}=[delta_K_pre,10.^ypredtrain];
    GP_results{2,i}=material_index;
    GP_results{3,i}=data_all{1,i}(1,3);
    Experiment_results{1,i}=[D_i(:,1),D_i(:,2)];
    Experiment_results{2,i}=material_index;
    Experiment_results{3,i}=data_all{1,i}(1,3);
end

%% 

numSlices = materials_N;
% 
Y_slices = linspace(1, numSlices, numSlices); 


planeAlpha = 0.1; 
planeColor = [0.8 0.8 0.8]; 
point_size=50;

delta_K_pre_all=[];
ypre_all=[];
for i=1:data_N
    delta_K_pre_all=[delta_K_pre_all;GP_results{1,i}(:,1)];
    ypre_all=[ypre_all;GP_results{1,i}(:,2)];
end


R_ALL_temp=[-0.2,0.3,0.5,0.1,0.2,0.7,0.4,0.6,0.8,-1.0,0.01];
R_ALL = unique(R_ALL_temp);

x_min = min(delta_K_pre_all)*0.9;
x_max = max(delta_K_pre_all)*1.1;
z_min = min(ypre_all)*0.9;
z_max = max(ypre_all)*1.1;


markerShapes = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', '+', '*','x'};
if numSlices > length(markerShapes)
    warning('Repeat markerShapes');
end

% 
%figure(1);
figure('Position', [100, 100, 600, 600]);
hold on;
% 
set(gca, 'XScale', 'log', 'ZScale', 'log');
% 
%[az, el] = view
view(10, 30); 

% 
for i = 1:numSlices
    Y_val = Y_slices(i);
    patch([x_min x_max x_max x_min], ...      
          [Y_val Y_val Y_val Y_val], ...      
          [z_min z_min z_max z_max], ...      
          planeColor, ...                      
          'FaceAlpha', planeAlpha, ...        
          'EdgeColor', 'k', ...               
          'LineWidth', 1);                     
end


for i = 1:data_N
    % 
    R=Experiment_results{3,i};
    R_index = find(R_ALL == R);
    marker = markerShapes{R_index};   

    material_index=Experiment_results{2,i}; 
    Experiment_results_sparse=Experiment_results{1,i}(1:4:end,:);
    N_data=length(Experiment_results_sparse);
    X_all=Experiment_results_sparse(:,1);
    Y_all=ones(N_data,1)*Y_slices(material_index);
    Z_all=Experiment_results_sparse(:,2);
    % 
    scatter3(X_all, Y_all, Z_all, point_size, colors(material_index,:),marker,'LineWidth', 2)
    %scatter3(X_all, Y_all, Z_all, point_size, colors(material_index,:), marker, 'filled');

    N_data=length(GP_results{1,i});
    X_curve=GP_results{1,i}(:,1);
    Y_curve=ones(N_data,1)*Y_slices(material_index);
    Z_curve=GP_results{1,i}(:,2);
    plot3(X_curve, Y_curve, Z_curve, 'r--', 'LineWidth', 2);

end


set(gca, 'YTick', [], 'YColor', 'none');


xlim([x_min, x_max]);
ylim([min(Y_slices)-0.05, max(Y_slices)+0.05]); 
zlim([z_min, z_max]);

% 
axis tight;

hold off;


figure
hold on;
% 
set(gca, 'XScale', 'log', 'YScale', 'log');
for i = 1:data_N
    R=Experiment_results{3,i};
    R_index = find(R_ALL == R);
    marker = markerShapes{R_index};   
    material_index=Experiment_results{2,i}; 
    Experiment_results_sparse=Experiment_results{1,i}(1:4:end,:);
    X_all=Experiment_results_sparse(:,1);
    Z_all=Experiment_results_sparse(:,2);
    % 
    scatter(X_all, Z_all, point_size, colors(material_index,:),marker,'LineWidth', 2)
end





end









