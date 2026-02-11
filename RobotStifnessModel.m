function RobotStifnessModel()



clear; clc;
pkg load control

%% =============================================================================
%% ------------------ Parameters -----------------------------------------------
%% =============================================================================

RPM = 2000;                % [rev/min]
z_teeth = 4;               % [-] number of teeth
Toolw = RPM*(2*pi/60)*z_teeth;    % [rad/s] tooth passing frequency

a     = [0.33 1.15 0.115 0 0 0];             % [m]
d     = [0.645 0 0 1.22 0 0.24];             % [m]
alpha = [pi/2 0 pi/2 pi/2 -pi/2 0];         % [rad]



lc    = [0.35 0.60 0.55 0.15 0.12 0.20];     % [m]
lc_y  = [0.00 0.00 0.00 0.00 0.00 0.00];     % [m]
lc_z  = [0.00 0.00 0.00 0.00 0.00 0.00];     % [m]

m     = [315 315 165 20 35 46];              % [kg]

Ixx = [ 12.0  18.0  10.5  0.40  0.80  1.10 ];  % [kg*m^2]
Iyy = [  9.0  14.0   8.0  0.30  0.60  0.90 ];  % [kg*m^2]
Izz = [  6.0   9.5   5.5  0.20  0.40  0.70 ];  % [kg*m^2]

Ixy = [ 0 0 0 0 0 0 ];  % [kg*m^2]
Ixz = [ 0 0 0 0 0 0 ];  % [kg*m^2]
Iyz = [ 0 0 0 0 0 0 ];  % [kg*m^2]

Ks = diag([1409800 400760 935280 360000 370000 380000]);    % [N*m/rad]
Cs = diag([2e1  2e2  2e1  5e1  5e1  2e1 ]);     % [N*m*s/rad]
g  = 9.81;                                       % [m/s^2]

%% =============================================================================
%% Trajectory (q1..q6)
%% =============================================================================
xpos = [ 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 ].';

p = [2.6,0,0.2;
2.4,0,0.2;
2.2,0,0.2;
2.0,0,0.2;
1.8,0,0.2;
1.6,0,0.2;
1.4,0,0.2;
1.2,0,0.2;
1.0,0,0.2];   %[m]

Rot = [0,pi/2,0];                               % [rad]

xpos = p(:,1);

Angcord = GetAngles(a, d, alpha, p, Rot);
Angcord(:,2)=-Angcord(:,2)
Angcord(:,3)=Angcord(:,3)+pi

n = 6;

%% =========================================================================
%% Modal analysis
%% =========================================================================
N = size(Angcord,1);

wn   = zeros(N,n);   % [rad/s]
wd   = zeros(N,n);   % [rad/s]
zeta = zeros(N,n);   % [-]

DomModeX = zeros(N,1);
DomModeY = zeros(N,1);

wn_domX = zeros(N,1);
wn_domY = zeros(N,1);

zeta_domX = zeros(N,1);
zeta_domY = zeros(N,1);

k_modal_X = zeros(N,1);
k_modal_Y = zeros(N,1);

Phi_all = cell(N,1);     % store tracked mode shapes
S_all   = cell(N-1,1);   % store MAC matrices between consecutive points
phi_prev = [];

% --- Numerical derivative step (for Kg) ---
epsK = 1e-5;

% --- Sanity check storage (for research plots / QA logs) ---
symErrD  = zeros(N,1);
minEigD  = zeros(N,1);
condD    = zeros(N,1);
Pgrav    = zeros(N,1);
Kg_sens  = zeros(N,1);  % simple sensitivity metric
Kx_static = zeros(N,1);
Ky_static = zeros(N,1);

for p = 1:N

    q0 = Angcord(p,:).';

    % --- Inertia matrix (numeric, exact) ---
    D = inertia_numeric(q0,a,d,alpha,lc,lc_y,lc_z,m,Ixx,Iyy,Izz,Ixy,Ixz,Iyz);

    % --- Gravity stiffness Kg (numeric linearization) ---
    Kg = zeros(n);
    for i=1:n  % Numeric derivative (central difference)
        dq = zeros(n,1);     % A column vector with zeros.
        dq(i)=epsK;          % we apply the perturbation just to the joint i
        Gp = gravity_numeric(q0+dq,a,d,alpha,lc,lc_y,lc_z,m,g);
        Gm = gravity_numeric(q0-dq,a,d,alpha,lc,lc_y,lc_z,m,g);
        Kg(:,i) = (Gp-Gm)/(2*epsK);
    end

    Kt = Ks + Kg;

    % --- Undamped modal problem (eig + mode tracking with MAC) ---
    [phi, W2] = eig(Kt, D);
    wn_tmp = sqrt(real(diag(W2)));   % n×1

    % --- Sort modes by increasing natural frequency (important!) ---
    [wn_tmp, idx_sort] = sort(wn_tmp, 'ascend');
    phi = phi(:, idx_sort);

    % Normalize modes (for MAC comparison)
    for k = 1:n
        phi(:,k) = phi(:,k) / sqrt(phi(:,k)'*D*phi(:,k));
    end

    if p == 1
        perm = 1:n;
    else
        % MAC matrix (0..1) because columns are normalized:
        % MAC = |phi_prev' * phi|^2
        S = abs(phi_prev.' * phi).^2;
        perm = best_perm(S);
        S_all{p-1} = S;
    end

    % Reorder according to tracking permutation
    phi    = phi(:,perm);
    wn_tmp = wn_tmp(perm);





% ================================

    % Store
    wn(p,:) = wn_tmp.';
    Phi_all{p} = phi;
    phi_prev = phi;

    % =========================================================
    % STATIC CARTESIAN STIFFNESS CHECK
    % =========================================================
    J_static = jacobian_numeric(q0,a,d,alpha);
    Jv_static = J_static(1:3,:);

    Kcart_static = inv(Jv_static / Ks * Jv_static');

    Kx_static(p) = Kcart_static(1,1);
    Ky_static(p) = Kcart_static(2,2);



    % ============================
    % Damped state-space model
    % ============================
    Z = zeros(n);
    I = eye(n);

    A = [ Z         I ;
         -D\Kt    -D\Cs ];

    lam_all = eig(A);

% Select damped eigenvalue for each mode by matching imag(lambda) ~ wn_tmp(k)
lam_c = lam_all(imag(lam_all) > 0);
lam_sel = zeros(n,1);

for k = 1:n
    [~, idx_min] = min(abs(abs(imag(lam_c)) - wn_tmp(k)));
    lam_sel(k) = lam_c(idx_min);
end


lam_sel = lam_sel(:).';
wd(p,:)   = abs(imag(lam_sel));  % [rad/s]
zeta(p,:) = -real(lam_sel) ./ sqrt(real(lam_sel).^2 + imag(lam_sel).^2); % [-]
zeta(p,:) = max(0, min(zeta(p,:), 0.99));
% =========================================================
% Dominant modal parameters for SLD (X and Y directions)
% =========================================================

J = jacobian_numeric(q0,a,d,alpha);
Jv = J(1:3,:);

PhiX = zeros(1,n);
PhiY = zeros(1,n);
WeightX = zeros(1,n);
WeightY = zeros(1,n);

for k = 1:n
    phi_k = phi(:,k);
    v_cart = Jv * phi_k;

    PhiX(k) = v_cart(1);
    PhiY(k) = v_cart(2);

    WeightX(k) = (PhiX(k)^2) / (wn_tmp(k)^2);
    WeightY(k) = (PhiY(k)^2) / (wn_tmp(k)^2);
end

[~, idxX] = max(WeightX);
[~, idxY] = max(WeightY);

DomModeX(p) = idxX;
DomModeY(p) = idxY;

wn_domX(p) = wn_tmp(idxX);
wn_domY(p) = wn_tmp(idxY);

zeta_domX(p) = zeta(p,idxX);
zeta_domY(p) = zeta(p,idxY);

k_modal_X(p) = wn_domX(p)^2 / (PhiX(idxX)^2);
k_modal_Y(p) = wn_domY(p)^2 / (PhiY(idxY)^2);



    %% =====================================================================
    %% Consistency / sanity checks (stored + optional warnings)
    %% =====================================================================
    symErrD(p) = norm(D - D.', 'fro') / max(1, norm(D,'fro'));
    Ds = (D + D.')/2;                 % symmetric part (for numerical safety)
    ev = eig(Ds);
    minEigD(p) = min(real(ev));
    condD(p) = cond(D);

    % Potential energy (gravity)
    Pgrav(p) = potential_energy_numeric(q0,a,d,alpha,lc,lc_y,lc_z,m,g);

    % Simple Kg sensitivity metric (compare epsK vs 10*epsK)
    Kg2 = zeros(n);
    eps2 = 10*epsK;
    for i=1:n
        dq = zeros(n,1);
        dq(i)=eps2;
        Gp = gravity_numeric(q0+dq,a,d,alpha,lc,lc_y,lc_z,m,g);
        Gm = gravity_numeric(q0-dq,a,d,alpha,lc,lc_y,lc_z,m,g);
        Kg2(:,i) = (Gp-Gm)/(2*eps2);
    end
    Kg_sens(p) = norm(Kg - Kg2,'fro') / max(1, norm(Kg2,'fro'));

    % Optional: warning prints only if something looks bad
    if symErrD(p) > 1e-8
        fprintf('WARN: D not symmetric enough at p=%d (symErr=%g)\n', p, symErrD(p));
    end
    if minEigD(p) <= 0
        fprintf('WARN: D not positive definite at p=%d (minEig=%g)\n', p, minEigD(p));
    end
    if condD(p) > 1e12
        fprintf('WARN: D ill-conditioned at p=%d (cond=%g)\n', p, condD(p));
    end
    if Kg_sens(p) > 0.05
        fprintf('WARN: Kg sensitive to eps at p=%d (sens=%g)\n', p, Kg_sens(p));
    end

end


%% =========================================================================
%% MAC plots
%% =========================================================================
MAC_diag = zeros(N-1,n);
MAC_max  = zeros(N-1,n);

for p = 2:N
    S = S_all{p-1};  % MAC matrix (0..1)
    for k = 1:n
        MAC_diag(p-1,k) = S(k,k);
        MAC_max(p-1,k)  = max(S(k,:));
    end
end

figure; hold on;
colors = lines(n);
for k = 1:n
    plot(xpos(2:end), MAC_diag(:,k), 'LineWidth', 2, 'Color', colors(k,:));
end
xlabel('X (m)');
ylabel('MAC [-]');
title('MAC continuity (diagonal, tracked modes)');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6');
grid on; ylim([0 1.05]);

figure; hold on;
for k = 1:n
    plot(xpos(2:end), MAC_max(:,k), 'LineWidth', 2, 'Color', colors(k,:));
end
xlabel('X (m)');
ylabel('MAC_{max} [-]');
title('MAC max per mode (best match between consecutive points)');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6');
grid on; ylim([0 1.05]);

%% =========================================================================
%% Plot wn (undamped)
%% =========================================================================
figure; hold on;
colors = lines(n);
for i=1:n
    plot(xpos, wn(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
plot(xpos, Toolw*ones(size(xpos)), 'k--', 'LineWidth', 2);
xlabel('X (m)');
ylabel('\omega_n  [rad/s]');
title('Natural frequencies (undamped)');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','Tool line');
grid on;

%% =========================================================================
%% Plot wd (damped)
%% =========================================================================
figure; hold on;
for i=1:n
    plot(xpos, wd(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
plot(xpos, Toolw*ones(size(xpos)), 'k--', 'LineWidth', 2);
xlabel('X (m)');
ylabel('\omega_d  [rad/s]');
title('Damped natural frequencies (state-space)');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6','Tool line');
grid on;

%% =========================================================================
%% Plot damping ratios
%% =========================================================================
figure; hold on;
for i=1:n
    plot(xpos, zeta(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
xlabel('X (m)');
ylabel('\zeta [-]');
title('Modal damping ratio along trajectory');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6');
grid on;

%% =========================================================================
%% Plot sanity checks (optional but recommended for research)
%% =========================================================================
figure; plot(xpos, symErrD, 'LineWidth', 2); grid on;
xlabel('X (m)'); ylabel('||D-D^T||/||D|| [-]');
title('D symmetry error along trajectory');

figure; plot(xpos, minEigD, 'LineWidth', 2); grid on;
xlabel('X (m)'); ylabel('min eig(Ds) [-]');
title('Minimum eigenvalue of symmetric part of D');

figure; semilogy(xpos, condD, 'LineWidth', 2); grid on;
xlabel('X (m)'); ylabel('cond(D) [-]');
title('Condition number of D');

figure; plot(xpos, Pgrav, 'LineWidth', 2); grid on;
xlabel('X (m)'); ylabel('Potential energy P [J]');
title('Gravity potential energy along trajectory');

figure; plot(xpos, Kg_sens, 'LineWidth', 2); grid on;
xlabel('X (m)'); ylabel('Kg sensitivity metric [-]');
title('Kg sensitivity to eps (compare eps and 10*eps)');

%% =========================================================================
%% Visualización de modos
%% =========================================================================
p_vis  = 1;     % postura
scale  = 0.6;   % exageración modal

plotModalShapes(Angcord, Phi_all, a, d, alpha, p_vis, scale)

%% =========================================================================
%% Plot FRF magnitude
%% =========================================================================

plotFRF(wn,zeta,m,xpos);

%% =========================================================================
%% Plot Cartesian influence of the modes of vibration
%% =========================================================================

figure;
plot(xpos, DomModeX,'o-','LineWidth',2); hold on;
plot(xpos, DomModeY,'s-','LineWidth',2);
xlabel('X (m)');
ylabel('Mode index');
legend('X direction','Y direction');
title('Dominant mode along trajectory');
grid on;

figure;
plot(xpos, wn_domX,'LineWidth',2); hold on;
plot(xpos, wn_domY,'LineWidth',2);
xlabel('X (m)');
ylabel('\omega_n [rad/s]');
legend('X direction','Y direction');
title('Dominant natural frequency');
grid on;

figure;
plot(xpos, zeta_domX,'LineWidth',2); hold on;
plot(xpos, zeta_domY,'LineWidth',2);
xlabel('X (m)');
ylabel('\zeta [-]');
legend('X direction','Y direction');
title('Dominant damping ratio');
grid on;

figure;
plot(xpos, k_modal_X,'LineWidth',2); hold on;
plot(xpos, k_modal_Y,'LineWidth',2);
xlabel('X (m)');
ylabel('k_{modal} [N/m]');
legend('X direction','Y direction');
title('Equivalent Cartesian modal stiffness');
grid on;

figure;
plot(xpos, Kx_static,'LineWidth',2); hold on;
plot(xpos, Ky_static,'LineWidth',2);
xlabel('X (m)');
ylabel('Static Cartesian stiffness [N/m]');
legend('X static','Y static');
title('Static stiffness check');
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

O{1} = [0;0;0]; % Origin of base reference frame
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

epsG = 1e-6; % We define the perturbation used to apply the mathematical formula of the derivative

for k = 1:n
    dqk = zeros(n,1); % A column vector with zeros.
    dqk(k) = epsG; % we apply the perturbation just to the joint k

    Pp = 0; % Potential energy with positive perturbation
    T  = eye(4); % Transformation Matrix storage.
    for i = 1:n
        T = T * dh_num(a(i),alpha(i),d(i),q(i) + dqk(i)); % Transformation matrix with perturbation in the joint k.
        R = T(1:3,1:3);
        O = T(1:3,4);
        rcm = [lc(i); lc_y(i); lc_z(i)]; % Center of mass coordinates for each link.
        Ocm = O + R*rcm;
        Pp = Pp + m(i)*g*Ocm(3); % Potential Energy with Influence on the joint k. (positive perturbation)
    end

    Pm = 0; % Potential energy with negative perturbation
    T  = eye(4); % Transformation Matrix storage.
    for i = 1:n
        T = T * dh_num(a(i),alpha(i),d(i),q(i) - dqk(i)); % Transformation matrix with perturbation in the joint k.
        R = T(1:3,1:3);
        O = T(1:3,4);
        rcm = [lc(i); lc_y(i); lc_z(i)]; % Center of mass coordinates for each link.
        Ocm = O + R*rcm;
        Pm = Pm + m(i)*g*Ocm(3); % Potential Energy with Influence on the joint k. (negative perturbation)
    end

    G(k) = (Pp - Pm) / (2*epsG); % Derivative of the potential energy for each joint. (Column vector)
end

end

%% =========================================================================
%% Potential energy (gravity)
%% =========================================================================
function P = potential_energy_numeric(q,a,d,alpha,lc,lc_y,lc_z,m,g)

n = length(q);
P = 0;
T = eye(4);

for i = 1:n
    T = T * dh_num(a(i),alpha(i),d(i),q(i));
    R = T(1:3,1:3);
    O = T(1:3,4);
    rcm = [lc(i); lc_y(i); lc_z(i)];
    Ocm = O + R*rcm;
    P = P + m(i)*g*Ocm(3);
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

%% =========================================================================
%% Global assignment per step (permutation)
%% =========================================================================
function perm = best_perm(S)
n = size(S,1);
P = perms(1:n);              % 720 x 6
score = zeros(size(P,1),1);

for k = 1:size(P,1)
    p = P(k,:);
    score(k) = sum(diag(S(:,p)));
end

[~,idx] = max(score);
perm = P(idx,:);
end


%% =========================================================================
%% GEOMETRIC JACOBIAN (end-effector)
%% =========================================================================
function J = jacobian_numeric(q,a,d,alpha)

n = length(q);
T = eye(4);

O = cell(n+1,1);
z = cell(n+1,1);

O{1} = [0;0;0];
z{1} = [0;0;1];

for i=1:n
    T = T * dh_num(a(i),alpha(i),d(i),q(i));
    R = T(1:3,1:3);
    O{i+1} = T(1:3,4);
    z{i+1} = R*[0;0;1];
end

Oe = O{n+1};

Jv = zeros(3,n);
Jw = zeros(3,n);

for j=1:n
    Jv(:,j) = cross(z{j}, Oe - O{j});
    Jw(:,j) = z{j};
end

J = [Jv; Jw];  % 6×n
end



