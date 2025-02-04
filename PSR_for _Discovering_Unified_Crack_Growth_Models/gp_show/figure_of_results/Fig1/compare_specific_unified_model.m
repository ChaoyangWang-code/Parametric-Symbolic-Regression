clc 
clear all
close all

data_all=[
3.00	2.56	2.1
8.8	    2.98	3.6
0.4	    2.76	4.2
];

colors = [
    0 0.4470 0.7410;      
    0.8500 0.3250 0.0980; 
    %0.9290 0.6940 0.1250; 
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


noise_level=0.2;

for i=1:3
    C=data_all(i,1)*1e-10;
    m=data_all(i,2);
    th=data_all(i,3);

    delta_K_temp=log10(th+0.7):0.02:log10(80);
    delta_K_line=log10(th+0.2):0.02:log10(120);
    delta_K=10.^(delta_K_temp);
    delta_K_line=10.^(delta_K_line);
    rn=size(delta_K,2);
    dadN=C.*( delta_K-th ).^m;
    dadN_line=C.*( delta_K_line-th ).^m;
    noise = noise_level * randn(1, rn); 
    dadN_temp=log10(dadN)+noise;
    dadN=10.^(dadN_temp);
    figure(1)
    if i<3
        plot(log10(delta_K),log10(dadN),'o','MarkerEdgeColor', colors(1, :),'MarkerFaceColor', colors(1, :) ,'HandleVisibility', 'off')
    else
        plot(log10(delta_K),log10(dadN),'o','MarkerEdgeColor', colors(1, :),'MarkerFaceColor', colors(1, :),'DisplayName', 'Combine of datasets')
    end
    hold on; 
    figure(2)
    plot(log10(delta_K),log10(dadN),'o','MarkerEdgeColor', colors(i, :),'MarkerFaceColor', colors(i, :))
    hold on;
    if i<3
        plot(log10(delta_K_line),log10(dadN_line),'r--','LineWidth',2.5,'HandleVisibility', 'off')
    else
        plot(log10(delta_K_line),log10(dadN_line),'r--','LineWidth',2.5,'DisplayName', 'A unified model')        
    end
end

figure(1)
C=0.49
m=3.12
th=2.03
delta_K_temp=log10(th+0.2):0.02:log10(120);
delta_K=10.^(delta_K_temp);
dadN=C*1e-10.*( delta_K-th ).^m;
plot(log10(delta_K),log10(dadN),'r-','LineWidth',2.5,'DisplayName', 'A specific model')
xlim([log10(1.5) log10(160)]); 
x_ticks = log10([1 10 100]);
x_labels = {'1', '10', '100'};               
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);
ylim([log10(1e-13) log10(1e-2)]);
y_ticks = log10([1e-14 1e-10 1e-6 1e-2 ]);
y_labels = {'10^{-14}', '10^{-10}', '10^{-6}', '10^{-2}'};                
set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels);
% lgd1 =legend('Location', 'southeast');
% % 
% lgd1.FontSize = 16; 
% lgd1.FontName = 'Times New Roman'; 
% 
ax = gca; 
ax.FontSize = 14; 
ax.FontName = 'Times New Roman'; 
% 
ax.Box = 'off'; 
ax.XAxis.LineWidth = 1.5; 
ax.YAxis.LineWidth = 1.5; 



figure(2)
xlim([log10(1.5) log10(160)]); 
x_ticks = log10([1 10 100]);
x_labels = {'1', '10', '100'};               
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);
ylim([log10(1e-13) log10(1e-2)]);
y_ticks = log10([1e-14 1e-10 1e-6 1e-2 ]);
y_labels = {'10^{-14}', '10^{-10}', '10^{-6}', '10^{-2}'};                
set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels);
% lgd2 =legend('Location', 'southeast');
% % 
% lgd2.FontSize = 16; 
% lgd2.FontName = 'Times New Roman'; 
% 
ax = gca; 
ax.FontSize = 14; 
ax.FontName = 'Times New Roman'; 
% 
ax.Box = 'off'; 
ax.XAxis.LineWidth = 1.5; 
ax.YAxis.LineWidth = 1.5; 





