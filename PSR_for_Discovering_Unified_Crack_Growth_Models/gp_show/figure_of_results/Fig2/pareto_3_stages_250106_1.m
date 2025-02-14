
clc
clear all

color = [
    191, 192, 255;       
    22, 122, 184;                
    249, 134, 18;   
    194, 47, 47; 
];

%figuresize=[100, 100, 1200, 1000]; %[left bottom width height]

%color = [1 0.38 0;1 0 0; 0 0 1;0.63 0.13 0.94;0.5 1 0;1 0.5 0.31; 0.53 0.15 0.34; 1 0 1;1 1 0; 0.45 0.29 0.07;0 0 1;0.03 0.18 0.33];  
colors=color./255;
Path='workdir';     %Select the result path,example:    Windows:'C:\Users\Desktop\test250106\'   linux:'/gs/home/gp_test/test1/';

%show_results_gen_i(figuresize,gp,foldername_gen,str)
MULTI_R_show_results_of_all_pareto_in_gen_i(Path,colors);
clc


function [ ] = MULTI_R_show_results_of_all_pareto_in_gen_i(Path,colors)

%% 
choose_stage=1;
markerSize=100;

%
transparency_1=0.5;   
transparency_2=0.5;   

transparency_pareto_1=1.0;   
transparency_pareto_2=1.0;   

if choose_stage==1

    gen_count=1;                   
    str=num2str(gen_count);
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);    
    gp=data_geni{gen_count}.gp;
    for i=1:gp.runcontrol.pop_size
        A(i,1)=gp.fitness.complexity(i,1);
        A(i,2)=gp.fitness.values_alldata{i,1}(1,1);
        A(i,3)=gp.fitness.values_alldata{i,1}(1,2);
        A(i,4)=gp.fitness.values_alldata{i,1}(1,3);
    end

    rows_with_inf = any(isinf(A), 2);

    A_cleaned = A(~rows_with_inf, :);
    B = unique(A_cleaned, 'rows', 'stable');
    B1_eff=[];
    B1_ineff=[];
    [numRows, ~] = size(B);
    B1_eff=B; 

    tour_ind = 1:numRows;
    tour_ind=tour_ind';
    mo = [B1_eff(:,1) B1_eff(:,2)];            
    rank = ndfsort_rank1(mo);
    rank1_tour_ind_temp = tour_ind(rank == 1);
    num_rank1 = numel(rank1_tour_ind_temp);
    pareto_comp=mo(rank1_tour_ind_temp,1);
    pareto_fit=mo(rank1_tour_ind_temp,2);
    p_all=[pareto_comp,pareto_fit];
    pu1=unique(p_all,'rows','stable');            

    clearvars -except Path colors markerSize transparency_1 transparency_2 B1_eff B1_ineff pu1 transparency_pareto_1 transparency_pareto_2
    %clear A gp rows_with_inf A_cleaned B numRows mo rank 

    gen_count=25;                  
    str=num2str(gen_count);
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);    
    gp=data_geni{gen_count}.gp;
    for i=1:gp.runcontrol.pop_size
        A(i,1)=gp.fitness.complexity(i,1);
        A(i,2)=gp.fitness.values_alldata{i,1}(1,1);
        A(i,3)=gp.fitness.values_alldata{i,1}(1,2);
        A(i,4)=gp.fitness.values_alldata{i,1}(1,3);
    end

    rows_with_inf = any(isinf(A), 2);

    A_cleaned = A(~rows_with_inf, :);
    B = unique(A_cleaned, 'rows', 'stable');
    B2_eff=[];
    B2_ineff=[];
    [numRows, ~] = size(B);
    B2_eff=B; 

    tour_ind = 1:numRows;
    tour_ind=tour_ind';
    mo = [B2_eff(:,1) B2_eff(:,2)];          
    rank = ndfsort_rank1(mo);
    rank1_tour_ind_temp = tour_ind(rank == 1);
    num_rank1 = numel(rank1_tour_ind_temp);
    pareto_comp=mo(rank1_tour_ind_temp,1);
    pareto_fit=mo(rank1_tour_ind_temp,2);
    p_all=[pareto_comp,pareto_fit];
    pu2=unique(p_all,'rows','stable');             
    figure
    hold on;
    % if ~isempty(B1_ineff)
    %     scatter(B1_ineff(:,1), B1_ineff(:,2), markerSize, [0.8, 0.8, 0.8], 'filled', 'o');
    % end
    scatter(B1_eff(:,1), B1_eff(:,2), markerSize, colors(1, :), 'filled', 'o','MarkerFaceAlpha', transparency_1, 'MarkerEdgeAlpha', transparency_2);   
    % if ~isempty(B2_ineff)
    %     scatter(B2_ineff(:,1), B2_ineff(:,2), markerSize, [0.8, 0.8, 0.8], 'filled', 'o');
    % end
    scatter(B2_eff(:,1), B2_eff(:,2), markerSize, colors(2, :), 'filled', 'o','MarkerFaceAlpha', transparency_1, 'MarkerEdgeAlpha', transparency_2);
    scatter(pu1(:,1), pu1(:,2), markerSize, colors(3, :),'filled', 'o','MarkerFaceAlpha', transparency_pareto_1, 'MarkerEdgeAlpha', transparency_pareto_2); 
    scatter(pu2(:,1), pu2(:,2), markerSize, colors(4, :),'filled', 'o','MarkerFaceAlpha', transparency_pareto_1, 'MarkerEdgeAlpha', transparency_pareto_2); 
    xlim([-10 550]);
elseif choose_stage==2
    MAX_of_Y=20;
    gen_count=1;                   
    str=num2str(gen_count);
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);    
    gp=data_geni{gen_count}.gp;
    for i=1:gp.runcontrol.pop_size
        A(i,1)=gp.fitness.complexity(i,1);
        A(i,2)=gp.fitness.values_alldata{i,1}(1,1);
        A(i,3)=gp.fitness.values_alldata{i,1}(1,2);
        A(i,4)=gp.fitness.values_alldata{i,1}(1,3);
    end

    rows_with_inf = any(isinf(A), 2);

    A_cleaned = A(~rows_with_inf, :);
    B = unique(A_cleaned, 'rows', 'stable');
    B1_eff=[];
    B1_ineff=[];
    [numRows, ~] = size(B);
    for i=1:numRows
        if (B(i,2)<gp.fitness.loss_MSE_limit)
            B1_eff=[B1_eff;B(i,:)];
        else
            B1_ineff=[B1_ineff;B(i,:)];
        end
    end

    B1_eff(B1_eff(:, 3) > MAX_of_Y, 3) = MAX_of_Y;
    B1_ineff(B1_ineff(:, 3) > MAX_of_Y, 3) = MAX_of_Y;    

    tour_ind = 1:numRows;
    tour_ind=tour_ind';
    mo = [B1_eff(:,1) B1_eff(:,3)];               
    rank = ndfsort_rank1(mo);
    rank1_tour_ind_temp = tour_ind(rank == 1);
    num_rank1 = numel(rank1_tour_ind_temp);
    pareto_comp=mo(rank1_tour_ind_temp,1);
    pareto_fit=mo(rank1_tour_ind_temp,2);
    p_all=[pareto_comp,pareto_fit];
    pu1=unique(p_all,'rows','stable');             

    clearvars -except Path colors markerSize transparency_1 transparency_2 B1_eff B1_ineff pu1 MAX_of_Y transparency_pareto_1 transparency_pareto_2

    gen_count=50;                  
    str=num2str(gen_count);
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);    
    gp=data_geni{gen_count}.gp;
    for i=1:gp.runcontrol.pop_size
        A(i,1)=gp.fitness.complexity(i,1);
        A(i,2)=gp.fitness.values_alldata{i,1}(1,1);
        A(i,3)=gp.fitness.values_alldata{i,1}(1,2);
        A(i,4)=gp.fitness.values_alldata{i,1}(1,3);
    end

    rows_with_inf = any(isinf(A), 2);

    A_cleaned = A(~rows_with_inf, :);
    B = unique(A_cleaned, 'rows', 'stable');
    B2_eff=[];
    B2_ineff=[];
    [numRows, ~] = size(B);
    for i=1:numRows
        if (B(i,2)<gp.fitness.loss_MSE_limit)
            B2_eff=[B2_eff;B(i,:)];
        else
            B2_ineff=[B2_ineff;B(i,:)];
        end
    end

    B2_eff(B2_eff(:, 3) > MAX_of_Y, 3) = MAX_of_Y;
    B2_ineff(B2_ineff(:, 3) > MAX_of_Y, 3) = MAX_of_Y;

    tour_ind = 1:numRows;
    tour_ind=tour_ind';
    mo = [B2_eff(:,1) B2_eff(:,3)];              
    rank = ndfsort_rank1(mo);
    rank1_tour_ind_temp = tour_ind(rank == 1);
    num_rank1 = numel(rank1_tour_ind_temp);
    pareto_comp=mo(rank1_tour_ind_temp,1);
    pareto_fit=mo(rank1_tour_ind_temp,2);
    p_all=[pareto_comp,pareto_fit];
    pu2=unique(p_all,'rows','stable');             

    figure
    hold on;
    % if ~isempty(B1_ineff)
    %     scatter(B1_ineff(:,1), B1_ineff(:,3), markerSize, [0.8, 0.8, 0.8], 'filled', 'o');
    % end
    scatter(B1_eff(:,1), B1_eff(:,3), markerSize, colors(1, :), 'filled', 'o','MarkerFaceAlpha', transparency_1, 'MarkerEdgeAlpha', transparency_2);   
    % if ~isempty(B2_ineff)
    %     scatter(B2_ineff(:,1), B2_ineff(:,3), markerSize, [0.8, 0.8, 0.8], 'filled', 'o');
    % end
    scatter(B2_eff(:,1), B2_eff(:,3), markerSize, colors(2, :), 'filled', 'o','MarkerFaceAlpha', transparency_1, 'MarkerEdgeAlpha', transparency_2);
    scatter(pu1(:,1), pu1(:,2), markerSize, colors(3, :),'filled', 'o','MarkerFaceAlpha', transparency_pareto_1, 'MarkerEdgeAlpha', transparency_pareto_2); 
    scatter(pu2(:,1), pu2(:,2), markerSize, colors(4, :),'filled', 'o','MarkerFaceAlpha', transparency_pareto_1, 'MarkerEdgeAlpha', transparency_pareto_2); 
    xlim([-10 550]);
    ylim([-1 22]);
else

    gen_count=1;                  
    str=num2str(gen_count); 
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);    
    gp=data_geni{gen_count}.gp;
    for i=1:gp.runcontrol.pop_size
        A(i,1)=gp.fitness.complexity(i,1);
        A(i,2)=gp.fitness.values_alldata{i,1}(1,1);
        A(i,3)=gp.fitness.values_alldata{i,1}(1,2);
        A(i,4)=gp.fitness.values_alldata{i,1}(1,3);
    end

    rows_with_inf = any(isinf(A), 2);

    A_cleaned = A(~rows_with_inf, :);
    B = unique(A_cleaned, 'rows', 'stable');
    B1_eff=[];
    B1_ineff=[];
    [numRows, ~] = size(B);
    for i=1:numRows
        if (B(i,2)<gp.fitness.loss_MSE_limit)&&(B(i,3)<gp.fitness.loss_VAR_limit)
            B1_eff=[B1_eff;B(i,:)];
        else
            B1_ineff=[B1_ineff;B(i,:)];
        end
    end

    tour_ind = 1:numRows;
    tour_ind=tour_ind';
    mo = [B1_eff(:,1) B1_eff(:,4)];              
    rank = ndfsort_rank1(mo);
    rank1_tour_ind_temp = tour_ind(rank == 1);
    num_rank1 = numel(rank1_tour_ind_temp);
    pareto_comp=mo(rank1_tour_ind_temp,1);
    pareto_fit=mo(rank1_tour_ind_temp,2);
    p_all=[pareto_comp,pareto_fit];
    pu1=unique(p_all,'rows','stable');           

    clearvars -except Path colors markerSize transparency_1 transparency_2 B1_eff B1_ineff pu1 transparency_pareto_1 transparency_pareto_2

    gen_count=60;                    
    str=num2str(gen_count); 
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);    
    gp=data_geni{gen_count}.gp;
    for i=1:gp.runcontrol.pop_size
        A(i,1)=gp.fitness.complexity(i,1);
        A(i,2)=gp.fitness.values_alldata{i,1}(1,1);
        A(i,3)=gp.fitness.values_alldata{i,1}(1,2);
        A(i,4)=gp.fitness.values_alldata{i,1}(1,3);
    end
 
    rows_with_inf = any(isinf(A), 2);

    A_cleaned = A(~rows_with_inf, :);
    B = unique(A_cleaned, 'rows', 'stable');
    B2_eff=[];
    B2_ineff=[];
    [numRows, ~] = size(B);
    for i=1:numRows
        if (B(i,2)<gp.fitness.loss_MSE_limit)&&(B(i,3)<gp.fitness.loss_VAR_limit)
            B2_eff=[B2_eff;B(i,:)];
        else
            B2_ineff=[B2_ineff;B(i,:)];
        end
    end

    tour_ind = 1:numRows;
    tour_ind=tour_ind';
    mo = [B2_eff(:,1) B2_eff(:,4)];                
    rank = ndfsort_rank1(mo);
    rank1_tour_ind_temp = tour_ind(rank == 1);
    num_rank1 = numel(rank1_tour_ind_temp);
    pareto_comp=mo(rank1_tour_ind_temp,1);
    pareto_fit=mo(rank1_tour_ind_temp,2);
    p_all=[pareto_comp,pareto_fit];
    pu2=unique(p_all,'rows','stable');          

    figure
    hold on;
    % if ~isempty(B1_ineff)
    %     scatter(B1_ineff(:,1), B1_ineff(:,4), markerSize, [0.8, 0.8, 0.8], 'filled', 'o');
    % end
    scatter(B1_eff(:,1), B1_eff(:,4), markerSize, colors(1, :), 'filled', 'o','MarkerFaceAlpha', transparency_1, 'MarkerEdgeAlpha', transparency_2);   
    % if ~isempty(B2_ineff)
    %     scatter(B2_ineff(:,1), B2_ineff(:,4), markerSize, [0.8, 0.8, 0.8], 'filled', 'o');
    % end
    scatter(B2_eff(:,1), B2_eff(:,4), markerSize, colors(2, :), 'filled', 'o','MarkerFaceAlpha', transparency_1, 'MarkerEdgeAlpha', transparency_2);
    scatter(pu1(:,1), pu1(:,2), markerSize, colors(3, :),'filled', 'o','MarkerFaceAlpha', transparency_pareto_1, 'MarkerEdgeAlpha', transparency_pareto_2); 
    scatter(pu2(:,1), pu2(:,2), markerSize, colors(4, :),'filled', 'o','MarkerFaceAlpha', transparency_pareto_1, 'MarkerEdgeAlpha', transparency_pareto_2); 
    xlim([15 240]);
    ylim([2.8 7.2]);
end

    ax = gca; 
    ax.FontName = 'Arial';  
    ax.FontSize = 18;       






end














