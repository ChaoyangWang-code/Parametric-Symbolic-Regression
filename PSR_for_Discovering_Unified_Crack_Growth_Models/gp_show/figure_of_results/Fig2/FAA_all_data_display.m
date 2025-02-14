
clc 
clear all

load('FAA_data_240605.mat')
data_num=size(testdata_temp,2);
material_num=testdata_temp{2,end};
% for i=2:data_num 
%     if testdata_temp{2,i}>material_num
%         material_num=testdata_temp{2,i};
%     end
% end

for i=1:material_num                      
    count_begin=1;
    clear material_i_data
    for j=1:data_num                 
        if testdata_temp{2,j}==i          
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
% color = [
%     0.00, 0.45, 0.74;  
%     0.85, 0.33, 0.10;  
%     0.93, 0.69, 0.13;  
%     0.49, 0.18, 0.56;  
%     0.47, 0.67, 0.19;  
%     0.30, 0.75, 0.93;  
%     0.64, 0.08, 0.18;  
%     0.00, 0.50, 0.00;  
%     1.00, 0.07, 0.65;  
%     0.30, 0.20, 0.20;  
%     0.60, 0.40, 0.70;  
%     1.00, 0.84, 0.00;  
%     0.00, 0.00, 0.50;  
% ];

color = [
    0   114 178;   % Navy Blue
    0   158 115;   % Teal
    86  180 233;   % Sky Blue
    230 159   0;   % Orange
    213 94    0;   % Vermillion
    204 121 167;   % Reddish Purple
    240 228  66;   % Yellow
    153 153 153;   % Gray
    100 140  40;   % Olive Green
    170 170 230;   % Lilac
    251 128 114;   % Coral
    179 222 105;   % Light Green
    166  86  40;   % Brown
];
color = color / 255;
colors=[color;color];

solidMarkers  = {'o','+','*','x','s','d','^','v','>','<','p','h'}; 
hollowMarkers = {'o','s','d','^','v','>','<','p','h'}; 


markerSize = 5;
lineWidth = 1.0;

solidAttr = cell(12,1);
for i = 1:12
    shape = solidMarkers{i};
    
    if ismember(shape, {'o','s','d','^','v','>','<','p','h'})
        
        faceColor = 'auto'; 
    else
       
        faceColor = 'none';
    end
    
    solidAttr{i} = struct('shape',shape, ...
                          'faceColor',faceColor, ...
                          'edgeColor','auto', ...
                          'markerSize',markerSize, ...
                          'lineWidth',lineWidth);
end

hollowAttr = cell(9,1);
for i = 1:9
    shape = hollowMarkers{i};
    hollowAttr{i} = struct('shape',shape, ...
                           'faceColor','none', ...
                           'edgeColor','auto', ...
                           'markerSize',markerSize, ...
                           'lineWidth',lineWidth);
end


%attrMarkers = [solidAttr; hollowAttr];

attrMarkers = [hollowAttr; solidAttr];

R_ALL=[];
N_DATA=size(testdata_temp,2);

for i=1:N_DATA
    R=testdata_temp{1,i}(1,3);
    R_ALL=[R_ALL,R];
end
unique_R = unique(R_ALL);

material_last_index=0;
%figure
for i=1:N_DATA
    material_index=testdata_temp{2,i};
    if material_index==material_last_index
        R_index=R_index+1; 
    else
        R_index=1;
    end

    data_i=testdata_temp{1,i};
    % R=data_i(1,3);
    % R_index=find(unique_R==R);
    delta_K=data_i(:,1);
    dadN=data_i(:,2);

    m = attrMarkers{R_index};
    if strcmp(m.faceColor,'auto')
        m.faceColor = colors(R_index,:);
    end    
    %figure(material_index)
    figure(1)
    % loglog(delta_K,dadN,'LineStyle', 'none','Marker',attrMarkers{R_index}.shape,'Color', colors(color_index, :))
    loglog(delta_K, dadN, ...
        'LineStyle', 'none', ...
        'Marker', m.shape, ...
        'MarkerEdgeColor', colors(R_index,:), ...
        'MarkerFaceColor', m.faceColor, ...
        'MarkerSize', m.markerSize, ...
        'LineWidth', m.lineWidth);    
    hold on;grid on;
    %title(['material-', num2str(R_index)])
    % xlabel('Coefficients', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    % ylabel('Value', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    % zlabel('Value', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    % title('Distribution of Coefficients', 'FontName', 'Times New Roman', 'FontWeight', 'normal');    
    % figure(2)   
    % scatter3(log10(data_i(:,1)),data_i(:,3),log10(data_i(:,2)),'o','MarkerEdgeColor', colors(color_index, :))
    % hold on;
    % %title(['material-', num2str(i)])
    % %xlabel('$\Delta K\ (MPa\sqrt{m})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14,'FontWeight', 'bold');
    % xlabel('$\mathrm{\Delta K}\ (\mathrm{ksi\sqrt{in}})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    % ylabel('R', 'FontName', 'Times New Roman','FontSize', 14, 'FontWeight', 'bold');
    % zlabel('$\mathrm{da/dN}\ (\mathrm{in/cycle})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14,'FontWeight', 'bold');    
    % title('Distribution of Coefficients', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    % legend({'data-1', 'data-2', 'data-3', 'data-4', 'data-5', 'data-6', 'data-7', 'data-8', 'data-9', 'data-10', 'data-11', 'data-12', 'data-13'})
    material_last_index=material_index;
    %     
    % set(gca, 'Color', 'none');
    % set(gcf, 'Color', 'none');
end
%view(15,30)

%xlabel('$\Delta K\ (ksi\sqrt{in})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold');
%ylabel('$\mathrm{da/dN}\ (\mathrm{in/cycle})$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold');
set(gca, 'FontSize', 14);  









