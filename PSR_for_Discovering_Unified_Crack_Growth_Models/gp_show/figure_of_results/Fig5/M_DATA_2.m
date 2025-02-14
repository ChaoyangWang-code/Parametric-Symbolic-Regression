clc
clear all
% 
model1 = [0.0249, 0.2993, 0.0277, 0.2287, 3]; % M_multi: TrainMSE, TrainVAR, TestMSE, TestVAR, PNUM
model2 = [0.0294, 0.0205, 13.0976, 8];         % M_single: TrainMSE, TestMSE, TestVAR, PNUM

color=[
68, 114, 196;
157, 195, 230;
248, 203, 173;
211, 122, 50;
207  93  92
];

color=color./255;

% 
figure('Position', [100, 100, 560, 420]);
width = 0.15; 
offset = 0.2; 
x = [1, 2.5];   


yyaxis left;

bar(x(1)-offset*2, model1(1), width, 'FaceColor', color(1,:),'EdgeColor', 'none'); hold on;
bar(x(1)-offset, model1(3), width, 'FaceColor', color(2,:),'EdgeColor', 'none');

bar(x(2)-offset*1.5, model2(1), width, 'FaceColor', color(1,:),'EdgeColor', 'none');
bar(x(2)-offset*0.5, model2(2), width, 'FaceColor', color(2,:),'EdgeColor', 'none');
ylim([0, 0.05]);


yyaxis right;
% 
h=bar(x(1), model1(2), width, 'FaceColor', color(5,:),'EdgeColor', 'none');
%hatchfill(h, 'cross', 'FaceColor', color(3,:)); 
bar(x(1)+offset, model1(4), width, 'FaceColor', color(3,:),'EdgeColor', 'none');
bar(x(1)+offset*2, model1(5), width, 'FaceColor', color(4,:),'EdgeColor', 'none');
% 
bar(x(2)+offset*0.5, model2(3), width, 'FaceColor', color(3,:),'EdgeColor', 'none');
bar(x(2)+offset*1.5, model2(4), width, 'FaceColor', color(4,:),'EdgeColor', 'none');
ylim([0, 15]);
% 
xlim([0.2 3.3]);
% xticks(x);
set(gca, 'XTick', {});
set(gca, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

grid on;