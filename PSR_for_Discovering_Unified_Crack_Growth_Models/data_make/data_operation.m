clear testdata

data_num=size(testdata_temp,2);
material_num=testdata_temp{2,end};
% for i=2:data_num 
%     if testdata_temp{2,i}>material_num
%         material_num=testdata_temp{2,i};
%     end
% end


for i=1:material_num                      %总共多少种材料
    count_begin=1;
    clear material_i_data
    for j=1:data_num                 %总共用多少条数据
        if testdata_temp{2,j}==i          %所有数据中查找属于材料i的数据
            if count_begin==1
                material_i_data=testdata_temp{1,j};
                count_begin=0;
            else
                material_i_data=[material_i_data;testdata_temp{1,j}];
            end
        end  
    end
    testdata{1,i}=material_i_data;
    testdata{2,i}=i;
end

%% plot
color = [
    0.00, 0.45, 0.74;  % 蓝色
    0.85, 0.33, 0.10;  % 红色
    0.93, 0.69, 0.13;  % 黄色
    0.49, 0.18, 0.56;  % 紫色
    0.47, 0.67, 0.19;  % 绿色
    0.30, 0.75, 0.93;  % 天蓝色
    0.64, 0.08, 0.18;  % 暗红
    0.00, 0.50, 0.00;  % 暗绿
    1.00, 0.07, 0.65;  % 桃红
    0.30, 0.20, 0.20;  % 棕色
    0.60, 0.40, 0.70;  % 浅紫
    1.00, 0.84, 0.00;  % 金色
    0.00, 0.00, 0.50;  % 海军蓝
];
colors=[color;color];
multimaterial_num=size(testdata,2);
%figure
for i=1:multimaterial_num
    data_all_material_i=testdata{1,i};
    delta_K=data_all_material_i(:,1);
    dadN=data_all_material_i(:,2);
    figure(1)
    loglog(delta_K,dadN,'o','Color', colors(i, :))
    hold on;grid on;
    title(['material-', num2str(i)])
    % xlabel('Coefficients', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    % ylabel('Value', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    % zlabel('Value', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    % title('Distribution of Coefficients', 'FontName', 'Times New Roman', 'FontWeight', 'normal');    
    figure(2)   
    scatter3(log10(data_all_material_i(:,1)),data_all_material_i(:,3),log10(data_all_material_i(:,2)),'o','MarkerEdgeColor', colors(i, :))
    hold on;
    %title(['material-', num2str(i)])
    %xlabel('$\Delta K\ (MPa\sqrt{m})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14,'FontWeight', 'bold');
    xlabel('$\mathrm{\Delta K}\ (\mathrm{ksi\sqrt{in}})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('R', 'FontName', 'Times New Roman','FontSize', 14, 'FontWeight', 'bold');
    zlabel('$\mathrm{da/dN}\ (\mathrm{in/cycle})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14,'FontWeight', 'bold');    
    %title('Distribution of Coefficients', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    %legend({'data-1', 'data-2', 'data-3', 'data-4', 'data-5', 'data-6', 'data-7', 'data-8', 'data-9', 'data-10', 'data-11', 'data-12', 'data-13'})
end
%view(15,30)


for i=1:multimaterial_num
    figure
    data_all_material_i=testdata{1,i};
    delta_K=data_all_material_i(:,1);
    dadN=data_all_material_i(:,2);
    loglog(delta_K,dadN,'o','Color', colors(i, :))
    hold on;grid on;
end











