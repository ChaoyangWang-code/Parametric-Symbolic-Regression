
clc
clear all

%% 
gen_count=25;
%pop_choose=2;

color = [
    191, 192, 255;       
    22, 122, 184;                
    249, 134, 18;   
    194, 47, 47; 
];
colors=color./255;


Path='workdir';      %Select the result path,example:    Windows:'C:\Users\Desktop\test250106\'  

%load('single_data1_all_pareto_g25.mat')
load('single_data2_all_pareto_g25.mat')


str=num2str(gen_count);
data_path=[Path, str,'\gen',str, '.mat'];   
pareto_path=[Path,'pareto_all_gen.mat'];  
load(data_path);
load(pareto_path);

comp=cell2mat(pareto_all_gen{gen_count,1}(:,5));
fit=cell2mat(pareto_all_gen{gen_count,1}(:,4));
figure
plot(comp,fit,'bo')

N_pareto=size(theta_train_all_pareto,1);

numSlices = 3;


Y_slices = linspace(1, numSlices, numSlices); 
% 
planeAlpha = 0.0; 
planeColor = [0.8 0.8 0.8]; 
point_size=50;

x_max = 0;
z_min = 10000000;
z_max = 0;
for pop_choose=1:numSlices
    theta_train=theta_train_all_pareto{pop_choose,1};
    theta_validation = theta_validation_all_pareto{pop_choose,1}; 
    theta_combine=[theta_train;theta_validation];
    % 
    if any(isinf(theta_combine(:))) || any(isnan(theta_combine(:)))
        a = 0;
    else
        gene_N=size(theta_combine,2);
        % 
        maxValue = max(theta_combine(:));     
        % 
        minValue = min(theta_combine(:));
        if gene_N>x_max
           x_max=gene_N;
        end
        if maxValue>z_max
           z_max=maxValue;
        end
        if minValue<z_min
           z_min=minValue;
        end
        if maxValue>100 || minValue<-100
            warning('z warning')

        end
    end
end

x_min = 1*0.9;
x_max = x_max;
z_min = z_min*1.1;
z_max = z_max*1.1;



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


for pop_choose = 1:numSlices
    theta_train=theta_train_all_pareto{pop_choose,1};
    theta_validation=theta_validation_all_pareto{pop_choose,1};
    num_coefficients = size(theta_train, 2);
    X_all=1:num_coefficients;
    Y_all=ones(1,num_coefficients)*Y_slices(pop_choose);
    scatter3(X_all, Y_all, theta_train, point_size, colors(2,:),'o','filled','LineWidth', 2)
    for j=1:3
        scatter3(X_all, Y_all, theta_validation(j,:), point_size, colors(4,:),'o','LineWidth', 2)
    end
end



theta_train_all_pareto=theta_train_all_pareto(end-numSlices+1:end, :);
theta_validation_all_pareto=theta_validation_all_pareto(end-numSlices+1:end, :);

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


for pop_choose = 1:numSlices
    theta_train=theta_train_all_pareto{pop_choose,1};
    theta_validation=theta_validation_all_pareto{pop_choose,1};
    num_coefficients = size(theta_train, 2);
    X_all=1:num_coefficients;
    Y_all=ones(1,num_coefficients)*Y_slices(pop_choose);
    scatter3(X_all, Y_all, theta_train, point_size, colors(2,:),'o','filled','LineWidth', 2)
    for j=1:3
        scatter3(X_all, Y_all, theta_validation(j,:), point_size, colors(4,:),'o','LineWidth', 2)
    end
end








