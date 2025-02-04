clc
clear all


%% 
gen_count=25;
pop_choose=1;
Path='workdir';      %Select the result path,example:    Windows:'C:\Users\Desktop\test250106\'  


str=num2str(gen_count);
data_path=[Path, str,'\gen',str, '.mat'];   
pareto_path=[Path,'pareto_all_gen.mat'];  
load(data_path);
load(pareto_path);
%load('D:\gp_test\gaosuan\test240626\test240626_2_22_2\25\gen25.mat')
syms x1 x2 c1 c2

pareto_ID=cell2mat(pareto_all_gen{gen_count,1}(:,1));

N_pareto=size(pareto_ID,1);

genes=pareto_all_gen{gen_count,1}{pop_choose,2}{1,2};
cons=pareto_all_gen{gen_count,1}{pop_choose,2}{1,3};



%% 
load('FAA_data_240605.mat')

multidata_validation=testdata;
gp.multidata=multidata_validation;


pop_ID=pareto_ID(pop_choose,1);
evalstrs_VALIDATION=tree2evalstr(gp.pop{pop_ID},gp);
[fitness,tempgp] = feval(fitfun,evalstrs_VALIDATION,gp);       
returnvals_validation{pop_choose,1} = tempgp.fitness.returnvalues;
%returnvals_validation{i,2}=fitness;


%% 
validation_results=returnvals_validation{pop_choose,1}{1,1}(1:3,:);
loss_MSE=max(validation_results(:,1));
gene=genes(pop_choose,1);
con=cons(pop_choose,1);
theta=validation_results(:,(4+con):(4+con+gene+1));
var_temp=var(theta);
loss_VAR=max(var_temp);
loss_num=gene+2;
fitness_validation(pop_choose,1)=loss_MSE;
fitness_validation(pop_choose,2)=loss_VAR;
fitness_validation(pop_choose,3)=loss_num;

%%

theta_all_temp=returnvals_validation{pop_choose,1}{1,1}(:,(4+con):(4+con+gene+1));
theta_all = [theta_all_temp(:, 1), theta_all_temp(:, end), theta_all_temp(:, 2:end-1)];    

theta_train_temp=pareto_all_gen{gen_count,1}{1,2}{1,1}((4+con):(4+con+gene+1));
theta_train = [theta_train_temp(:, 1), theta_train_temp(:, end), theta_train_temp(:, 2:end-1)];


theta_validation=theta_all(1:3,:);
num_coefficients = size(theta_train, 2);
% 
coeff_labels = arrayfun(@(x) ['f' num2str(x)], 1:num_coefficients, 'UniformOutput', false);
%h2 = figure('Visible', 'on');

for i=1:1
    h1=plot(1:num_coefficients, theta_train(i, :), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'b');
    hold on; grid on;
end


for i = 1:3
    h2=plot(1:num_coefficients, theta_validation(i, :), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
end


%set(gca, 'LineWidth', 2); 
%xlabel('Parameters', 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
%ylabel('Value', 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
%title(['Distribution of Coefficients'], 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
% 
% xticks(1:num_coefficients);
% xticklabels(coeff_labels);
% 
%set(gca, 'Box', 'on', 'LineWidth', 2);
% grid on;
% 
set(gca, 'XTick', {});
set(gca, 'TickLength', [0.02, 0.025]);  
set(gca, 'TickDir', 'in');             
set(gca, 'XColor', [0 0 0]);             
set(gca, 'YColor', [0 0 0]);            
set(gca, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

xlim([0.5, 8.5]);  
ylim([-18, 10]);  
set(gca, 'GridLineStyle', '-'); 
set(gca, 'GridColor', [0.6 0.6 0.6]); 
%legend([h1, h2], {'train set', 'test set'}, 'Location', 'southeast');


























