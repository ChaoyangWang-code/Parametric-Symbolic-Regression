clc
clear all

load('FAA_NASGRO_experiment_data.mat')
load('NASGRO_FIT.mat')
load('WALKER_FIT.mat')
load('FAA_data_SPLINE_FIT_240418.mat')


%%
 load('D:\gp_test\gaosuan\test240626\test240626_2_2\60\gen60.mat')
% F=returnvals{221,1}{1,5}{1,1};
%load('D:\gp_test\test240626\test240626_2_17\25\gen25.mat')

pop_ID=221;

syms x1 x2 c1 c2

F=returnvals{pop_ID,1}{1,5}{1,1};
FF=simplify(eval(F))


theta=returnvals{pop_ID,1}{1,1}(:,2:6);
k1=theta(:,1);
k2=theta(:,2);
f1=theta(:,3)-theta(:,4).*log10(k1);
f2=theta(:,4);
f3=theta(:,5)+theta(:,4);
newf=[f1,f2,f3];

R=-1:0.01:1;

for i=1:13
    kapa=(1-R).*k2(i)/k1(i);
    KAPA=kapa+1./kapa;
    U=KAPA.^(f2(i)/f3(i));
    SS=1-U.*(1-R);
    figure(1)
    plot(R,SS,'-','LineWidth',1.5)
    hold on; 
    %grid on;
    %xlabel('R', 'FontName', 'Times New Roman', 'FontWeight', 'normal','FontSize', 18);
    %ylabel('So/Smax', 'FontName', 'Times New Roman', 'FontWeight', 'normal','FontSize', 18);
    %title('GP','FontName', 'Times New Roman', 'FontWeight', 'normal','FontSize', 18);  
    %set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
end

axis([-1,1,0,1])
% 
ax = gca; 
ax.XAxisLocation = 'origin'; 
ax.YAxisLocation = 'origin'; 

ax.XTick = -1:0.2:1; 
ax.YTick = 0:0.1:1; 
box off;
set(gcf,'color','w');
% 
%xlabel('R', 'FontName', 'Times New Roman', 'FontWeight', 'normal','FontSize', 16);
%ylabel('So/Smax', 'FontName', 'Times New Roman', 'FontWeight', 'normal','FontSize', 16);
%text(-0.95, 0.95, '\alpha  = 1 (Plane stress)', 'FontSize', 9, 'FontWeight', 'bold');

%% f_NASGRO
alpha = 2;
n = [0.3]'; 
A0 = (0.825-0.34*alpha+0.05*alpha^2)*(cos(pi/2*n)).^(1/alpha);
A1 = (0.415-0.071*alpha)*n;
A3 = 2*A0+A1-1;
A2 = 1-A0-A1-A3;
%
R1 = linspace(-1,0,100);
R2 = linspace(0,1,100);
S0_Smax2 = A0*ones(1,100)+A1*R2 + A2*R2.^2 +A3*R2.^3;
S0_Smax1 = A0*ones(1,100)+A1*R1;
%
S0_Smax = [S0_Smax1(:,1:end-1),S0_Smax2];
R = [R1(:,1:end-1),R2];

figure(2)
for i = 1:length(n)
    plot(R,S0_Smax(i,:),'r',"LineWidth",1.5)
    hold on
end



alpha = 1;
n = [0.0 0.2 0.4 0.6]'; 
A0 = (0.825-0.34*alpha+0.05*alpha^2)*(cos(pi/2*n)).^(1/alpha);
A1 = (0.415-0.071*alpha)*n;
A3 = 2*A0+A1-1;
A2 = 1-A0-A1-A3;
%
R1 = linspace(-1,0,100);
R2 = linspace(0,1,100);
S0_Smax2 = A0*ones(1,100)+A1*R2 + A2*R2.^2 +A3*R2.^3;
S0_Smax1 = A0*ones(1,100)+A1*R1;
%
S0_Smax = [S0_Smax1(:,1:end-1),S0_Smax2];
R = [R1(:,1:end-1),R2];

for i = 1:length(n)
    plot(R,S0_Smax(i,:),'k--',"LineWidth",1.0)
    hold on
end




axis([-1,1,0,1])
% 
ax = gca; 
ax.XAxisLocation = 'origin'; 
ax.YAxisLocation = 'origin'; 

% 
ax.XTick = -1:0.2:1; 
ax.YTick = 0:0.1:1; 
box off;
set(gcf,'color','w');
% 
%xlabel('R', 'FontName', 'Times New Roman', 'FontWeight', 'normal','FontSize', 16);
%ylabel('So/Smax', 'FontName', 'Times New Roman', 'FontWeight', 'normal','FontSize', 16);
%text(0.2, 0.3, '\alpha  = 1 (Plane stress)', 'FontName', 'Times New Roman', 'FontWeight', 'normal','FontSize', 16);





