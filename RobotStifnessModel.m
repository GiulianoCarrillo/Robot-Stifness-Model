function RobotStifnessModel()

clear; clc;

%% =============================================================================
%% Trajectory (q1..q6)
%% =============================================================================
Angcord1 = [ ...
   -3.97   -9.32    36.81   -91.37   -92.63   332.47;
   -4.30  -21.87    62.29   -92.14   -92.51   319.54;
   -4.69  -30.27    79.89   -92.81   -92.39   310.33;
   -5.16  -36.97    94.43   -93.51   -92.23   302.47;
   -5.73  -42.59   107.24   -94.27   -92.02   295.27;
   -6.44  -47.36   118.92   -95.16   -91.72   288.36;
   -7.36  -51.30   129.77   -96.23   -91.27   281.45;
   -8.57  -54.23   139.98   -97.55   -90.56   274.22;
  -10.26  -55.70   149.59   -99.24   -89.37   266.16;
  -12.77  -54.53   158.52  -101.43   -87.17   256.30 ];

xpos = [ 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8];

Angcord = deg2rad(Angcord1);

n = 6;

%% =============================================================================
%% ------------------ Parameters -----------------------------------------------
%% =============================================================================

RPM = 2000
z = 4

Toolw = 1500*(2*pi/60)*4


a     = [0.33 1.15 0.115 0 0 0];
d     = [0.645 0 0 1.22 0 0.24];
alpha = [pi/2 0 -pi/2 pi/2 -pi/2 0];

lc    = [0.35 0.60 0.55 0.15 0.12 0.20];
lc_y  = [0.00 0.00 0.00 0.00 0.00 0.00];
lc_z  = [0.00 0.00 0.00 0.00 0.00 0.00];

m     = [315 315 165 20 35 46];

Ixx = [ 12.0  18.0  10.5  0.40  0.80  1.10 ];
Iyy = [  9.0  14.0   8.0  0.30  0.60  0.90 ];
Izz = [  6.0   9.5   5.5  0.20  0.40  0.70 ];

Ixy = [ 0 0 0 0 0 0 ];
Ixz = [ 0 0 0 0 0 0 ];
Iyz = [ 0 0 0 0 0 0 ];

Ks = diag([20e6 20e6 20e6 50e4 50e4 20e4]);
Cs= diag ([20e2 20e2 20e2 50e2 50e2 20e2]);
g  = 9.81;

%% =========================================================================
%% Modal analysis
%% =========================================================================
N = size(Angcord,1);
wn = zeros(N,n);
wd   = zeros(N,n);   % [rad/s] Dampened frecuencies
zeta = zeros(N,n);   % [-]     Dampened factor

eps = 1e-5;

for p = 1:N

    q0 = Angcord(p,:)';

    % --- Inertia matrix (numeric, exact) ---
    D = inertia_numeric(q0,a,d,alpha,lc,lc_y,lc_z,m,Ixx,Iyy,Izz,Ixy,Ixz,Iyz);

    % --- Gravity stiffness (numeric linearization) ---
    Kg = zeros(n);

    G0 = gravity_numeric(q0,a,d,alpha,lc,lc_y,lc_z,m,g);

    for i=1:n %  Numeric derivate
        dq = zeros(n,1); % A column vector with zeros.
        dq(i)=eps; %we apply the perturbation just to the joint k
        Gp = gravity_numeric(q0+dq,a,d,alpha,lc,lc_y,lc_z,m,g); %valuation of the gravity function with a perturbation just in the joiint k.
        Gm = gravity_numeric(q0-dq,a,d,alpha,lc,lc_y,lc_z,m,g);
        Kg(:,i) = (Gp-Gm)/(2*eps); % Derivate
    end

    Kt = Ks + Kg;



    % --- Modal problem ---
    [~,W2] = eig(Kt,D);
    wn(p,:) = sqrt(real(diag(W2)))';

    % ============================
% Estado-espacio amortiguado
% ============================
Z = zeros(n);
I = eye(n);

A = [ Z         I ;
     -D\Kt    -D\Cs ];

lam_all = eig(A);

% 1) Candidatos complejos (solo una mitad del par conjugado)
lam_c = lam_all(imag(lam_all) > 0);

% 2) Si faltan modos (polos reales), completamos con reales
if numel(lam_c) < n
    lam_r = lam_all(abs(imag(lam_all)) < 1e-12);   % reales
    % ordenarlos por |real| (los más “lentos” primero o como prefieras)
    [~,idxr] = sort(abs(real(lam_r)), 'ascend');
    lam_r = lam_r(idxr);
    lam_sel = [lam_c; lam_r(1:min(n-numel(lam_c), numel(lam_r)))];
else
    % si sobran complejos, elegir los n con mayor frecuencia (|imag|)
    [~,idxc] = sort(abs(imag(lam_c)), 'descend');
    lam_sel = lam_c(idxc(1:n));
end

% asegurar tamaño n
lam_sel = lam_sel(:).';
if numel(lam_sel) < n
    lam_sel(end+1:n) = 0;
end

wd(p,:)   = abs(imag(lam_sel));                          % [rad/s]
zeta(p,:) = -real(lam_sel) ./ sqrt(real(lam_sel).^2 + imag(lam_sel).^2);  % [-]
              % [-]

end

%% =========================================================================
%% Plot wn (no amortiguadas)
%% =========================================================================
figure; hold on;
colors = lines(6);
for i=1:6
    plot(xpos(:), wn(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
plot(xpos(:), Toolw*ones(size(xpos(:))), 'k--', 'LineWidth', 2);

xlabel('X (m)');
ylabel('\omega_n  [rad/s]');
title('Natural frequencies (undamped)');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','Tool line');
grid on;

%% =========================================================================
%% Plot wd (amortiguadas)
%% =========================================================================
figure; hold on;
colors = lines(6);
for i=1:6
    plot(xpos(:), wd(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
plot(xpos(:), Toolw*ones(size(xpos(:))), 'k--', 'LineWidth', 2);

xlabel('X (m)');
ylabel('\omega_d  [rad/s]');
title('Damped natural frequencies');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','Tool line');
grid on;

%% =========================================================================
%% Damped factor
%% =========================================================================

figure; hold on;
colors = lines(6);
for i=1:6
    plot(xpos(:), zeta(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
xlabel('X (m)');
ylabel('\zeta [-]');
title('Modal damping ratio along trajectory');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6');
grid on;

end

%% =========================================================================
%% NUMERIC INERTIA MATRIX
%% =========================================================================
function D = inertia_numeric(q,a,d,alpha,lc,lc_y,lc_z,m,Ixx,Iyy,Izz,Ixy,Ixz,Iyz)

n = length(q);
T = eye(4); %identity matrix

O = cell(n+1,1); %cell or empty arrange for origin of reference systems
z = cell(n+1,1); %cell or empty arrange for z-axes
R = cell(n,1); %cell or empty arrange for rotation matrix

O{1} = [0;0;0]; % Origin of base reference frame2im
z{1} = [0;0;1]; % Z axes of the reference frame.

for i=1:n
    T = T * dh_num(a(i),alpha(i),d(i),q(i));
    R{i}   = T(1:3,1:3);
    O{i+1} = T(1:3,4);
    z{i+1} = R{i}*[0;0;1];
end

D = zeros(n);

for i=1:n
    rcm = [lc(i); lc_y(i); lc_z(i)];
    Ocm = O{i+1} + R{i}*rcm;

    Jv = zeros(3,n);
    Jw = zeros(3,n);

    for j=1:i
        Jv(:,j) = cross(z{j},Ocm-O{j});
        Jw(:,j) = z{j};
    end

    Ii = [ Ixx(i) -Ixy(i) -Ixz(i);
          -Ixy(i)  Iyy(i) -Iyz(i);
          -Ixz(i) -Iyz(i)  Izz(i) ];

    D = D + m(i)*(Jv'*Jv) + Jw'*R{i}*Ii*R{i}'*Jw;
end
end

%% =========================================================================
%% GRAVITY MATRIX
%% =========================================================================
function G = gravity_numeric(q,a,d,alpha,lc,lc_y,lc_z,m,g)

n = length(q);
G = zeros(n,1);

epsG = 1e-6; %We define the perturbation used to apply the mathematical formula of the derivate

for k = 1:n
    dqk = zeros(n,1); % A column vector with zeros.
    dqk(k) = epsG; %we apply the perturbation just to the joint k

    Pp = 0; % Potential energy with positive perturbation
    T  = eye(4); % Transformation Matrix storage.
    for i = 1:n
        T = T * dh_num(a(i),alpha(i),d(i),q(i) + dqk(i)); % Transformation matrix with perturbation in the joint k.
        R = T(1:3,1:3);
        O = T(1:3,4);
        rcm = [lc(i); lc_y(i); lc_z(i)]; %Center of mass coordinates for each link.
        Ocm = O + R*rcm;
        Pp = Pp + m(i)*g*Ocm(3); %Potential Energi with Influence on the joint k. (positive perturbation)
    end

    Pm = 0; %Potential energy with negative perturbation
    T  = eye(4); % Transformation Matrix storage.
    for i = 1:n
        T = T * dh_num(a(i),alpha(i),d(i),q(i) - dqk(i)); % Transformation matrix with perturbation in the joint k.
        R = T(1:3,1:3);
        O = T(1:3,4);
        rcm = [lc(i); lc_y(i); lc_z(i)]; %Center of mass coordinates for each link.
        Ocm = O + R*rcm;
        Pm = Pm + m(i)*g*Ocm(3); %Potential Energi with Influence on the joint k.  (negative perturbation)
    end

    G(k) = (Pp - Pm) / (2*epsG); % Derivate of the potential Energi For each joint. (Column vector)
end

end

%% =========================================================================
%% DH TRANSFORM
%% =========================================================================
function T = dh_num(a,alpha,d,theta)
T = [ ...
 cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
 sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
 0           sin(alpha)             cos(alpha)            d;
 0           0                      0                     1];
end

