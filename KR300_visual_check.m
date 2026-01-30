function KR300_DH_frames_visual()
% ============================================================
% KR300 R2700 – DH CLASSIC
% Visualizacion de FRAMES DH ubicados EN LAS JUNTAS
% (no barras físicas, no dinámica, solo DH correcto)
% ============================================================

clc; close all;

%% Joint values (rad)
q = [0 -pi/2 -pi/2 0 0 0];

%% DH parameters (Craig) – FRAMES AT JOINTS
% i :  ai [m]   di [m]    alpha_i
a     = [0.330   1.150   0.115   0     0     0.240];
d     = [0.645   0       0       1.220 0     0];
alpha = [pi/2    0      pi/2    pi/2 -pi/2  0];

%% Forward kinematics – DH
T = eye(4);

O = zeros(3,7);      % origins of frames {0}..{6}
R = cell(7,1);       % rotations of frames

O(:,1) = [0;0;0];
R{1}   = eye(3);

for i = 1:6
    T = T * dh_T(a(i), alpha(i), d(i), q(i));
    O(:,i+1) = T(1:3,4);
    R{i+1}   = T(1:3,1:3);
end

%% Plot frames only
figure('Color','w'); hold on; grid on; axis equal;
view(90,0);   % lateral X–Z view

xlabel('X [m]');
ylabel('Z [m]');
title('KR300 R2700 – DH frames at joints');

scale = 0.15;

for i = 1:7
    draw_frame(O(:,i), R{i}, scale, i-1);
end

end

%% ============================================================
%% DH homogeneous transform (Craig)
%% ============================================================
function T = dh_T(a,alpha,d,theta)
T = [ cos(theta)  -sin(theta)*cos(alpha)   sin(theta)*sin(alpha)  a*cos(theta);
      sin(theta)   cos(theta)*cos(alpha)  -cos(theta)*sin(alpha)  a*sin(theta);
      0            sin(alpha)               cos(alpha)             d;
      0            0                        0                      1 ];
end

%% ============================================================
%% Draw coordinate frame
%% ============================================================
function draw_frame(O, R, s, idx)

x = O + s*R(:,1);
y = O + s*R(:,2);
z = O + s*R(:,3);

plot([O(1) x(1)], [O(3) x(3)], 'r', 'LineWidth', 2); % X
plot([O(1) y(1)], [O(3) y(3)], 'g', 'LineWidth', 2); % Y
plot([O(1) z(1)], [O(3) z(3)], 'b', 'LineWidth', 2); % Z

text(O(1), O(3), sprintf('  {%d}',idx), 'FontSize', 10);

end

