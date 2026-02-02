function KR300_frames_KUKA_style()
% ============================================================
% KR300 R2700 – KUKA-style kinematics
% Frames físicos en las juntas
% Transformaciones homogéneas puras (NO DH)
% ============================================================

clc; close all;

%% =============================
%% Joint values (rad)
%% =============================
q = [0 -pi/2 pi/2 0 0 0];   % cambiar para testear

%% =============================
%% Geometry (from your drawing)
%% All expressed in LOCAL frames
%% =============================

%% =============================
%% Forward kinematics
%% =============================

% Base → A1
T01 = transl(0.330, 0, 0.645)*rotx(pi/2);

% A1 → A2
T12 =  transl(0, 1.150, 0) ;

% A2 → A3
T23 = transl(-0.115, 0, 0) * rotx(-pi/2)* rotz(pi/2);

% A3 → A4
T34 = transl(0, 0, 1.220)* rotx(pi/2) ;

% A4 → A5
T45 = rotx(-pi/2);

% A5 → A6 (flange)
T56 = transl(0, 0, 0.240);

T_offset = {T01,T12,T23,T34,T45,T56};

%% =============================
%% Forward kinematics
%% =============================
T = eye(4);

O = zeros(3,7);
R = cell(7,1);

O(:,1) = [0;0;0];
R{1} = eye(3);

for i = 1:6
    T = T * T_offset{i} * rotz(q(i));
    O(:,i+1) = T(1:3,4);
    R{i+1}   = T(1:3,1:3);
end

%% =============================
%% Plot
%% =============================
figure('Color','w'); hold on; grid on; axis equal;
view(3);

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('KR300 R2700 – KUKA-style frames at joints');

scale = 0.15;

for i = 1:7
    draw_frame(O(:,i), R{i}, scale, i-1);
end

end

%% ============================================================
%% Helpers
%% ============================================================

function T = transl(x,y,z)
T = eye(4);
T(1:3,4) = [x;y;z];
end

function T = rotx(a)
T = [1 0 0 0;
     0 cos(a) -sin(a) 0;
     0 sin(a)  cos(a) 0;
     0 0 0 1];
end

function T = rotz(a)
T = [cos(a) -sin(a) 0 0;
     sin(a)  cos(a) 0 0;
     0       0      1 0;
     0       0      0 1];
end

function draw_frame(O,R,s,idx)
x = s*R(:,1);
y = s*R(:,2);
z = s*R(:,3);

plot3([O(1) O(1)+x(1)], [O(2) O(2)+x(2)], [O(3) O(3)+x(3)], 'r', 'LineWidth',2);
plot3([O(1) O(1)+y(1)], [O(2) O(2)+y(2)], [O(3) O(3)+y(3)], 'g', 'LineWidth',2);
plot3([O(1) O(1)+z(1)], [O(2) O(2)+z(2)], [O(3) O(3)+z(3)], 'b', 'LineWidth',2);

text(O(1),O(2),O(3),sprintf('{%d}',idx),'FontSize',10);
end

