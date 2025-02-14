

clc
clear all;
close all;

load('single_data_compare_gen25_testdata.mat')

for data_i=1:13
    validation_results=single_data_compare{data_i,1}{2,2};
    rows = size(validation_results, 1);
    %model_index=1;
    for i=1:rows
        if all(isfinite(validation_results(i, :))) && all(isreal(validation_results(i, :)))
            model_index=i;
            break
        end
    end
    A(data_i,1)=single_data_compare{data_i,1}{2,1}(model_index,1);    
    A(data_i,2)=single_data_compare{data_i,1}{2,2}(model_index,1);    
    A(data_i,3)=single_data_compare{data_i,1}{2,2}(model_index,2);    
    A(data_i,4)=single_data_compare{data_i,1}{2,2}(model_index,3);    
    A(data_i,5)=model_index;    
end
a=[0.0261	0.0261	0.6893 3 1];
A=[A;a];


A(:, 3) = min(A(:, 3), 10);      %Set the truncation value of L_VAR


train_MSE = A(:, 1);
test_MSE = A(:, 2);
L_VAR = A(:, 3);
L_PNUM = A(:, 4);


x = 1:size(A, 1);
width = 0.18; 


figure;


yyaxis left;
h1 = bar(x - width, train_MSE, width, 'FaceColor', [68, 114, 196]./255, 'EdgeColor', 'none', 'DisplayName', '$\mathcal{L}_{\mathrm{MSE}} (\mathrm{Train \, data})$');
hold on;
h2 = bar(x, test_MSE, width, 'FaceColor', [157, 195, 230]./255, 'EdgeColor', 'none', 'DisplayName', '$\mathcal{L}_{\mathrm{MSE}} (\mathrm{Test \, data})$');
ylabel(['$\mathcal{L}_{\mathrm{MSE}}$'], ...
       'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylim([0, max([train_MSE; test_MSE])*1.8]);


yyaxis right;
h3 = bar(x + width, L_VAR, width, 'FaceColor', [248, 203, 173]./255, 'EdgeColor', 'none', 'DisplayName', '$\mathcal{L}_{\mathrm{VAR}} (\mathrm{Test \, data})$');
h4 = bar(x + 2*width, L_PNUM, width, 'FaceColor', [211, 122, 50]./255, 'EdgeColor', 'none', 'DisplayName', '$\mathcal{L}_{\mathrm{PNUM}}$');
ylabel(['$\mathcal{L}_{\mathrm{VAR}} \, \& \, \mathcal{L}_{\mathrm{PNUM}}$'], ...
       'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylim([0, max([L_VAR; L_PNUM])*1.2]);


legend([h1, h2, h3, h4], {'$\mathcal{L}_{\mathrm{MSE}} (\mathrm{Train \, data})$', '$\mathcal{L}_{\mathrm{MSE}} (\mathrm{Test \, data})$', '$\mathcal{L}_{\mathrm{VAR}} (\mathrm{Test \, data})$', '$\mathcal{L}_{\mathrm{PNUM}}$'}, ...
       'Interpreter', 'latex', 'Location', 'northwest', 'NumColumns', 2, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');


xlabel('$\mathrm{Model_{out}}$ for Training Data', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
grid on;


set(gca, 'XTick', x);
xticklabels = get(gca, 'XTickLabel'); 
xticklabels{end} = 'All'; 
set(gca, 'XTickLabel', xticklabels); 






