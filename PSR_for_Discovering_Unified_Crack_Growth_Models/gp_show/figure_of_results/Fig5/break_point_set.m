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

A_sorted = sortrows(A, 3);


log10(A_sorted(:,3))




























