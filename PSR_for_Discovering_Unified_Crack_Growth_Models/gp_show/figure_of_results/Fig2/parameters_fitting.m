clear; clc; close all;

% 
theta1 = linspace(-3,3,100);
theta2 = linspace(-3,3,100);
[Theta1, Theta2] = meshgrid(theta1, theta2);

L = peaks(Theta1, Theta2);


control_colors = [0.2 0.7 0.6;   
                  0.3 0.8 0.7;   
                  0.4 0.85 0.8;  
                  0.5 0.6 0.8;   
                  0.6 0.5 0.7];  
num_colors = 256;
cmap = interp1(linspace(0,1,size(control_colors,1)), control_colors, linspace(0,1,num_colors));

fig = figure('Color','w');
ax = axes(fig);

colormap(ax, cmap);

s = surf(ax, Theta1, Theta2, L, 'EdgeColor','none', 'FaceAlpha', 0.95);
shading interp;
hold on;

%camlight left; lighting phong;
xlabel('\theta_1','FontSize',12);
ylabel('\theta_2','FontSize',12);
zlabel('L(\theta)','FontSize',12);
% title('3D Surface with 2D Contour Projection','FontSize',14);

floorZ = min(L(:)) - 3;

contour_levels = 20; 
[C,hContour] = contour3(Theta1, Theta2, L, contour_levels, 'k');

hContour.ZData = floorZ * ones(size(hContour.XData));


floorSurface = surf(Theta1, Theta2, floorZ*ones(size(L)), L, ...
    'EdgeColor','none', 'FaceColor','interp', 'FaceAlpha',0.8);



[Fx, Fy] = gradient(peaks(Theta1, Theta2), theta1(2)-theta1(1), theta2(2)-theta2(1));


initial_points = [-2 -2; 2 2; -2 2; 2 -2];
num_steps = 10;
step_size = 0.2;

path_colors = [1 0.5 0;    
               0.9 0.3 0;  
               1 0.4 0.1;  
               1 0.6 0.2]; 

for i = 1:size(initial_points,1)
    t_current = initial_points(i,:);
    traj = t_current; 
    for step = 1:num_steps
        [gx, gy] = getGradient(t_current(1), t_current(2), theta1, theta2, Fx, Fy);
        t_new = t_current - step_size * [gx gy];
        traj = [traj; t_new];
        t_current = t_new;
    end
    traj_loss = peaks(traj(:,1), traj(:,2));
    
end

view(ax, [45 50]);
axis off; 

% % 
function [gx, gy] = getGradient(t1, t2, theta1, theta2, Fx, Fy)
    [~, i1] = min(abs(theta1 - t1));
    [~, i2] = min(abs(theta2 - t2));
    gx = Fx(i2, i1);
    gy = Fy(i2, i1);
end






