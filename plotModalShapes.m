function plotModalShapes(Angcord, Phi_all, a, d, alpha, p, scale)
% =========================================================================
% plotModalShapes
%
% Visualiza los modos de vibración de un robot serial (stick model)
% usando perturbaciones articulares:
%
%   q_def = q0 + scale * Phi(:,k)
%
% INPUTS:
%   Angcord : [N x 6] matriz de posturas (rad)
%   Phi_all : cell{N} con matrices de modos [6 x 6]
%   a,d,alpha : parámetros DH
%   p      : índice de postura a visualizar
%   scale  : factor de exageración modal (rad)
%
% =========================================================================

close all;

q0  = Angcord(p,:)';     % postura nominal
Phi = Phi_all{p};        % modos en esa postura
n   = length(q0);

for k = 1:n

    % Postura deformada según el modo k
    q_def = q0 + scale * Phi(:,k);

    figure('Color','w'); hold on;
    axis equal; grid on;
    view(3);

    % Robot nominal
    plot_robot(q0, a, d, alpha, 'k-o', 2);

    % Robot deformado
    plot_robot(q_def, a, d, alpha, 'r-o', 2);

    title(sprintf('Modo %d – postura %d (escala %.2f)', k, p, scale), ...
          'FontSize', 12);

    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');

    legend('Robot nominal','Modo deformado', 'Location','best');

end

end
function plot_robot(q, a, d, alpha, style, lw)
% Dibuja el robot como barras entre articulaciones

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

