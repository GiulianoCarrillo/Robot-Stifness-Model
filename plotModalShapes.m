function plotModalShapes(Angcord, Phi_all, a, d, alpha, p, scale)
% =========================================================================
% plotModalShapes
%
%   q_def = q0 + scale * Phi(:,k)
%
% INPUTS:
%   Angcord : [N x 6] posture matrix (rad)
%   Phi_all : cell{N} Matrix modes [6 x 6]
%   a,d,alpha : DH parameters
%   p      : Index of the posture to visualize
%
% =========================================================================


q0  = Angcord(p,:)';     % Nominal posture
Phi = Phi_all{p};        % Modes in the indicated posture
n   = length(q0);

for k = 1:n

    % Deformated posture from mode k
    q_def = q0 + scale * Phi(:,k);

    figure('Color','w'); hold on;
    axis equal; grid on;
    view(3);

    % Nominal posture of the robot
    plot_robot(q0, a, d, alpha, 'k-o', 2);

    % Robot after deformation.
    plot_robot(q_def, a, d, alpha, 'r-o', 2);

    title(sprintf('Mode %d – Posture %d (scale %.2f)', k, p, scale), ...
          'FontSize', 12);

    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');

    legend('Nomimnal Robot','Robot after deformation', 'Location','best');

end

end
function plot_robot(q, a, d, alpha, style, lw)
% Draw the robot with lines between articulations.

n = length(q);
T = eye(4);

P = zeros(3,n+1);
P(:,1) = [0;0;0];

for i = 1:n
    T = T * dh_num(a(i), alpha(i), d(i), q(i));
    P(:,i+1) = T(1:3,4);
end

plot3(P(1,:), P(2,:), P(3,:), style, ...
      'LineWidth', lw, 'MarkerSize', 6);

end
function T = dh_num(a,alpha,d,theta)

T = [ ...
 cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
 sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
 0           sin(alpha)             cos(alpha)            d;
 0           0                      0                     1];

end

