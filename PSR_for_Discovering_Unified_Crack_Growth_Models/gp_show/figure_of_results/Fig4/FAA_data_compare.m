clc
clear all
load('dadN_NASGRO.mat')
load('FAA_NASGRO_experiment_data.mat')
Num_R=size(dadN_NASGRO,2);
%Num_R=66;

Path='workdir';      %Select the result path,example:    Windows:'C:\Users\Desktop\test250106\'  
gen_count=60;
pareto_i=1;

%Path='C:\Users\wzy59\Desktop\10_1\';
show_figure='on';  


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


str=num2str(gen_count);
data_path=[Path, str,'\gen',str, '.mat'];   
data_geni=load(data_path);
foldername_gen=[Path, str];
gp=data_geni.gp;


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
returns=gp.fitness.returnvalues;



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
foldername_pareto_i_fig=[foldername_gen,'\','pareto_',num2str(pareto_i),'_fig'];
if ~exist(foldername_pareto_i_fig, 'dir')
    mkdir(foldername_pareto_i_fig);
end        
foldername_pareto_i_png=[foldername_gen,'\','pareto_',num2str(pareto_i),'_png'];
if ~exist(foldername_pareto_i_png, 'dir')
    mkdir(foldername_pareto_i_png);
end 


%% 
GP_ypredtrain_ALLDATA=[];
NASGRO_ypredtrain_ALLDATA=[];
yexperiment_ALLDATA=[];

GP_results=cell(1,13);
NASGRO_results=cell(1,13);
Experiment_results=cell(1,13);
material_last_index=0;
for i=1:Num_R   
    material_index=FAA_NASGRO_experiment_data_all_R{2,i};

    index=dadN_NASGRO{1,i}(:,3);
    data_sort=sortrows(FAA_NASGRO_experiment_data_all_R{1,i},1);
    compare_data2_temp=data_sort(:,2);
    compare_data2=compare_data2_temp(index);                    
    material_index=FAA_NASGRO_experiment_data_all_R{2,i};
    parameter_K=parameter_now(material_index,2:3);
    delta_kth=parameter_K(1);
    kc=parameter_K(2);
    Const_pair_now=parameter_now(material_index,4:3+num_const);
    theta=parameter_now(material_index,4+num_const:4+num_const+numGenes+1);           
    deltaK_temp=data_sort(:,1);                 
    deltaK=deltaK_temp(index);
    Kmax_temp=data_sort(:,4);
    Kmax=Kmax_temp(index);
    xtrain=[deltaK./delta_kth,Kmax./kc];
    num_of_data_i=size(deltaK,1);
    pat1 = 'x(\d+)';
    pat2 = 'c(\d+)';
    F_eq = regexprep(F_eq,pat2,'Const_pair_now($1)');
    F_eq = regexprep(F_eq,pat1,'xtrain(:,$1)');
    geneOutputs = ones(num_of_data_i,numGenes+2);
    for j = 1:numGenes
        ind = j + 1;
        gene_temp=eval_exe(F_eq{j},xtrain,Const_pair_now);
        geneOutputs(:,ind)=lg(gene_temp);
    end
    
    geneOutputs(:,numGenes+2)=lg(deltaK);               
    ypredtrain = geneOutputs * theta';

    compare_data1=10.^(ypredtrain);             
    compare_data3=dadN_NASGRO{1,i}(:,2);        

    if material_index==material_last_index
        GP_material_i=[GP_material_i;compare_data1]; 
        NASGRO_material_i=[NASGRO_material_i;compare_data3]; 
        Experiment_material_i=[Experiment_material_i;compare_data2]; 
        if i==Num_R
            GP_results{1,material_index}=GP_material_i;
            NASGRO_results{1,material_index}=NASGRO_material_i;
            Experiment_results{1,material_index}=Experiment_material_i;
        end
    else
        if material_last_index~=0
            
            GP_results{1,material_last_index}=GP_material_i;
            NASGRO_results{1,material_last_index}=NASGRO_material_i;
            Experiment_results{1,material_last_index}=Experiment_material_i;
        end
        GP_material_i=compare_data1;
        NASGRO_material_i=compare_data3;
        Experiment_material_i=compare_data2;
    end


    GP_ypredtrain_ALLDATA=[GP_ypredtrain_ALLDATA;compare_data1];
    NASGRO_ypredtrain_ALLDATA=[NASGRO_ypredtrain_ALLDATA;compare_data3];
    yexperiment_ALLDATA=[yexperiment_ALLDATA;compare_data2];

    material_last_index=material_index;   


end

GP_ypredtrain_ALLDATA=log10(GP_ypredtrain_ALLDATA);
NASGRO_ypredtrain_ALLDATA=log10(NASGRO_ypredtrain_ALLDATA);
Experiment_ALLDATA=log10(yexperiment_ALLDATA);

%% R2
GP_R2=1-sum((GP_ypredtrain_ALLDATA-Experiment_ALLDATA).^2)/sum((GP_ypredtrain_ALLDATA-mean(Experiment_ALLDATA)).^2)
NASGRO_R2=1-sum((NASGRO_ypredtrain_ALLDATA-Experiment_ALLDATA).^2)/sum((NASGRO_ypredtrain_ALLDATA-mean(Experiment_ALLDATA)).^2)


%% 


numSlices = material_index;

Y_slices = linspace(1, numSlices, numSlices); 

planeAlpha = 0.0; 
planeXZ45=0.3;
% 
planeColor = [0.8 0.8 0.8]; 
plane45Color = [0 0.4470 0.7410]; 

point_size = 36;

diagonal_slope = 1; 
diagonal_intercept = 0; 

log_min=min(NASGRO_ypredtrain_ALLDATA);
log_max=max(NASGRO_ypredtrain_ALLDATA);

x_min = 10^log_min*0.9;
x_max = 10^log_max*1.1;
z_min = 10^log_min*0.9;
z_max = 10^log_max*1.1;


%figure(1);
figure('Position', [100, 100, 600, 600]);
hold on;

% 
%[az, el] = view
view(10, 30); 


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


for i = 1:numSlices
    % 
    N_data=length(GP_results{1,i});
    X_all=GP_results{1,i};
    Y_all=ones(N_data,1)*Y_slices(i);
    Z_all=Experiment_results{1,i};
    scatter3(X_all, Y_all, Z_all, point_size, colors(i,:), 'filled');
    %scatter3(X_all, Y_all, Z_all, point_size, colors(i,:));
end


[X_plane, Y_plane] = meshgrid(linspace(x_min, x_max, 2), linspace(min(Y_slices)-0.05, max(Y_slices)+0.05, 2));
Z_plane = X_plane; 
surf(X_plane, Y_plane, Z_plane, 'FaceColor', plane45Color, ...    
                             'FaceAlpha', planeXZ45, ...       
                             'EdgeColor', 'none');        


set(gca, 'XScale', 'log', 'ZScale', 'log');


set(gca, 'YTick', [], 'YColor', 'none');

% 
%xlabel('GP-Predictions');
%zlabel('Experiments');

ax = gca;  
ax.XAxis.FontName = 'Times New Roman';
ax.XAxis.FontSize = 12;
ax.YAxis.FontName = 'Times New Roman';
ax.YAxis.FontSize = 12;

% 
axis tight;
hold off;


% 
figure('Position', [100, 100, 600, 600]);
hold on;
% 
%[az, el] = view
view(10, 30); 

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

% colors = hsv(numSlices); 

% 
for i = 1:numSlices
    % 
    N_data=length(NASGRO_results{1,i});
    X_all=NASGRO_results{1,i};
    Y_all=ones(N_data,1)*Y_slices(i);
    Z_all=Experiment_results{1,i};
    scatter3(X_all, Y_all, Z_all, point_size, colors(i,:), 'filled');
    %scatter3(X_all, Y_all, Z_all, point_size, colors(i,:));
end


[X_plane, Y_plane] = meshgrid(linspace(x_min, x_max, 2), linspace(min(Y_slices)-0.05, max(Y_slices)+0.05, 2));
Z_plane = X_plane; 
surf(X_plane, Y_plane, Z_plane, 'FaceColor', plane45Color, ...    
                             'FaceAlpha', planeXZ45, ...       
                             'EdgeColor', 'none');        


set(gca, 'XScale', 'log', 'ZScale', 'log');


set(gca, 'YTick', [], 'YColor', 'none');

% 
%xlabel('NASGRO-Predictions');
%zlabel('Experiments');

ax = gca;  
ax.XAxis.FontName = 'Times New Roman';
ax.XAxis.FontSize = 12;
ax.YAxis.FontName = 'Times New Roman';
ax.YAxis.FontSize = 12;

% 
axis tight;
hold off;




