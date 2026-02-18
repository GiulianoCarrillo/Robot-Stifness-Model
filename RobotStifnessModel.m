function RobotStifnessModel()

% =========================================================================
% RobotStifnessModel.m
%
% Dynamic stiffness + modal model for a 6-DOF industrial robot along a
% Cartesian trajectory.
%
% Pipeline (high level):
%   1) Define robot parameters: DH geometry, link COM, masses, inertias,
%      joint stiffness/damping.
%   2) Define Cartesian target points p [m] and constant TCP orientation Rot [rad].
%   3) Inverse kinematics (GetAngles): returns joint angles q [rad] for each target.
%      - Uses an analytic seed for joints 1–3 (wrist center geometry)
%      - Refines full 6-DOF solution with numerical IK (damped least squares).
%   4) Forward-kinematics sanity check: verifies FK(q) ≈ target p.
%   5) For each posture q:
%      - D(q): inertia matrix [kg·m²]
%      - G(q): gravity torque vector [N·m]
%      - Kg = dG/dq: gravity stiffness (linearization) [N·m/rad]
%      - Kt = Ks + Kg: total joint stiffness [N·m/rad]
%      - Undamped eigenproblem: Kt φ = ωn² D φ
%      - Mode tracking along the path using MAC (modal assurance criterion)
%      - Damped state-space eigenvalues → ωd and damping ratio ζ
%      - Dominant mode in Cartesian X/Y via Jacobian-based weighting
%
% Units:
%   - Length: m
%   - Angle: rad
%   - Mass: kg
%   - Inertia: kg·m²
%   - Joint stiffness Ks: N·m/rad
%   - Joint damping  Cs: N·m·s/rad
%   - Natural frequency ω: rad/s   (f = ω / 2π in Hz)
%   - Cartesian stiffness (modal/static): N/m
% =========================================================================
% Key computed outputs (per posture index p):
%   wn(p,k)      : undamped natural frequencies [rad/s]
%   wd(p,k)      : damped natural frequencies [rad/s] (from state-space)
%   zeta(p,k)    : damping ratios [-]
%   Phi_all{p}   : tracked joint-space mode shapes [6x6] (mass-normalized)
%   DomModeX/Y   : dominant mode index in Cartesian X / Y [-]
%   wn_domX/Y    : dominant undamped natural frequency in X / Y [rad/s]
%   zeta_domX/Y  : damping ratio of the dominant mode in X / Y [-]
%   k_modal_X/Y  : equivalent Cartesian modal stiffness in X / Y [N/m]
%   Kx_static/Y  : static Cartesian stiffness check in X / Y [N/m]
% =========================================================================

clear; clc;
pkg load control

%% =============================================================================
%% ------------------ Parameters -----------------------------------------------
%% =============================================================================

% Toolw is used as a reference excitation line in the frequency plots.

RPM = 2000;                % [rev/min]
z_teeth = 4;               % [-] number of teeth
Toolw = RPM*(2*pi/60)*z_teeth;    % [rad/s] tooth passing frequency

% DH convention (as used by dh_num):
% T_i = RotZ(theta_i) * TransZ(d_i) * TransX(a_i) * RotX(alpha_i)
% where theta_i = q(i). All kinematics/dynamics functions below assume this.


a     = [0.33 1.15 0.115 0 0 0];             % [m]
d     = [0.645 0 0 1.22 0 0.24];             % [m]
alpha = [pi/2 0 pi/2 pi/2 -pi/2 0];         % [rad]

% Link center-of-mass offsets in each link frame (DH frame i):
% r_cm_i = [lc(i); lc_y(i); lc_z(i)] [m]
% In this dataset only the x-offset is used (lc_y = lc_z = 0).

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

% Ks: diagonal joint stiffness matrix (torsional spring at each joint).
% Cs: diagonal joint viscous damping matrix (linearized joint damping).

Ks = diag([1409800 400760 935280 360000 370000 380000]);    % [N*m/rad]
Cs = diag([2e1  2e2  2e1  5e1  5e1  2e1 ]);     % [N*m*s/rad]
g  = 9.81;                                       % [m/s^2]

%% =============================================================================
%% Trajectory (q1..q6)
%% =============================================================================

% Trajectory definition:
%   p(i,:) = [x y z] target TCP position in base frame [m].
%   Here y is kept at 0 and z at 0.2 m, while x is swept from 2.6 → 1.0 m.

p = [2.6,0,0.2;
2.4,0,0.2;
2.2,0,0.2;
2.0,0,0.2;
1.8,0,0.2;
1.6,0,0.2;
1.4,0,0.2;
1.2,0,0.2;
1.0,0,0.2];   %[m]

p_traj = p;

% --- Orientation convention used in this project ---
% Rot = [psi, theta, phi] = [Rx, Ry, Rz] in radians (roll, pitch, yaw).
% The TCP rotation matrix is built as: R = Rz(phi) * Ry(theta) * Rx(psi)  (ZYX order).

Rot = [0,pi/2,-pi/2];


%% =============================================================================

xpos = p(:,1);

Angcord = GetAngles(a, d, alpha, p, Rot);

% Forward kinematics check:
%   For each IK solution Angcord(i,:), compute FK via dh_num and compare
%   with the requested target position p(i,:). This is the quickest way to
%   detect convention mismatches (DH / tool offset / orientation).

fprintf("\n=== TCP FK check ===\n");
for ii=1:size(Angcord,1)
    T = eye(4);
    for jj=1:6
        T = T * dh_num(a(jj), alpha(jj), d(jj), Angcord(ii,jj));
    end
    fk = T(1:3,4).';
    fprintf("i=%d target=[%.3f %.3f %.3f] fk=[%.3f %.3f %.3f]\n", ...
        ii, p(ii,1),p(ii,2),p(ii,3), fk(1),fk(2),fk(3));
end
fprintf("====================\n\n");

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

% === Modal constants at TCP for SLD (per posture, per mode) ===
PhiX_all = zeros(N,n);   % b_x = (Jv*phi)_x
PhiY_all = zeros(N,n);   % b_y = (Jv*phi)_y

uX = zeros(N,n);         % u_xj = b_xj^2  (units ~ 1/kg)
uY = zeros(N,n);         % u_yj = b_yj^2

uXY = zeros(N,n);        % cross term b_xj*b_yj (useful if someone wants Gxy)

k_eq_X = zeros(N,n);     % k_xj = wn^2 / u_xj  [N/m]
k_eq_Y = zeros(N,n);     % k_yj = wn^2 / u_yj  [N/m]

c_eq_X = zeros(N,n);     % c_xj = 2*zeta*wn / u_xj [N*s/m]
c_eq_Y = zeros(N,n);     % c_yj = 2*zeta*wn / u_yj [N*s/m]


for p = 1:N

% =====================================================================
% Per-posture loop
%   For each joint configuration q0 = Angcord(p,:), we build the joint-space
%   dynamic model and extract modal parameters.
%
%   Kt = Ks + Kg includes gravity stiffness (linearization of gravity torque).
%   Undamped modes come from:  Kt * phi = wn^2 * D * phi
%   where D is the inertia matrix.
% =====================================================================

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

    % Generalized eigenvalue problem:
    % eig(Kt, D) returns eigenvectors phi and eigenvalues W2 = wn^2.
    % The mode shapes are then mass-normalized using D for MAC tracking.

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

    % Damped second-order system (joint space):
    %   D*qdd + Cs*qd + Kt*q = 0
    % Written in state-space with x = [q; qd]:
    %   xdot = A*x  where
    %   A = [0 I; -D\Kt  -D\Cs]

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
% Modal constants at TCP for SLD (X/Y) + dominant indices (optional)
% =========================================================
J = jacobian_numeric(q0,a,d,alpha);
Jv = J(1:3,:);

% b = Jv * phi  (3x6).  With phi mass-normalized (phi'*D*phi = I),
% the collocated receptance in X/Y is:
% Gxx(ω) = Σ b_xj^2 / (wn_j^2 - ω^2 + 2i zeta_j wn_j ω)
% Gyy(ω) = Σ b_yj^2 / (wn_j^2 - ω^2 + 2i zeta_j wn_j ω)

B = Jv * phi;              % 3x6
PhiX = B(1,:);             % 1x6  (b_xj)
PhiY = B(2,:);             % 1x6  (b_yj)

PhiX_all(p,:) = PhiX;
PhiY_all(p,:) = PhiY;

% Modal constants for collocated FRFs
uX(p,:)  = PhiX.^2;         % u_xj = b_xj^2
uY(p,:)  = PhiY.^2;         % u_yj = b_yj^2
uXY(p,:) = PhiX.*PhiY;      % cross term if needed (Gxy)

% Avoid division by zero
epsu = 1e-18;
uX(p,:) = max(uX(p,:), epsu);
uY(p,:) = max(uY(p,:), epsu);

% Equivalent (per-mode) stiffness and damping in X/Y (handy to export)
wn_row = wn_tmp(:).';       % 1x6
z_row  = zeta(p,:);         % 1x6

k_eq_X(p,:) = (wn_row.^2) ./ uX(p,:);              % [N/m]
k_eq_Y(p,:) = (wn_row.^2) ./ uY(p,:);              % [N/m]
c_eq_X(p,:) = (2*z_row.*wn_row) ./ uX(p,:);        % [N*s/m]
c_eq_Y(p,:) = (2*z_row.*wn_row) ./ uY(p,:);        % [N*s/m]

% Optional: "dominant" by static-compliance contribution (NOT a SLD rule)
WeightX = uX(p,:) ./ (wn_row.^2);
WeightY = uY(p,:) ./ (wn_row.^2);

[~, idxX] = max(WeightX);
[~, idxY] = max(WeightY);

DomModeX(p) = idxX;
DomModeY(p) = idxY;

wn_domX(p) = wn_row(idxX);
wn_domY(p) = wn_row(idxY);

zeta_domX(p) = z_row(idxX);
zeta_domY(p) = z_row(idxY);

k_modal_X(p) = k_eq_X(p, idxX);
k_modal_Y(p) = k_eq_Y(p, idxY);



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


% =========================================================================
% ==== SLD OUTPUT: save/extrct each value per position =========
% =========================================================================

SLD = struct();

% trayectoria y eje de barrido
SLD.p_xyz   = p_traj;     % [N x 3]
SLD.xpos    = xpos;       % [N x 1]

% modos dominantes
SLD.DomModeX = DomModeX;  % [N x 1]
SLD.DomModeY = DomModeY;  % [N x 1]

% frecuencias dominantes (rad/s) y también en Hz
SLD.wn_domX_rad = wn_domX;                 % [N x 1]
SLD.wn_domY_rad = wn_domY;                 % [N x 1]
SLD.fn_domX_Hz  = wn_domX/(2*pi);          % [N x 1]
SLD.fn_domY_Hz  = wn_domY/(2*pi);          % [N x 1]

% amortiguamiento dominante
SLD.zeta_domX = zeta_domX;  % [N x 1]
SLD.zeta_domY = zeta_domY;  % [N x 1]

% rigidez modal equivalente cartesiana
SLD.k_modal_X = k_modal_X;  % [N x 1]  (N/m)
SLD.k_modal_Y = k_modal_Y;  % [N x 1]  (N/m)

% opcional: chequeo de rigidez estática
SLD.Kx_static = Kx_static;  % [N x 1]
SLD.Ky_static = Ky_static;  % [N x 1]

% 1) Dejarlo disponible en Workspace aunque esto sea una función:
assignin('base','SLD',SLD);

% 2) Guardarlo a disco (MAT):
save(fullfile(tempdir(),'SLD_results.mat'),'SLD');

% =========================================================================
% EXPORT for colleague (SLD): ONE ROW = (posture, mode)
% Use ';' delimiter so Excel (ES locale) doesn't break decimals
% =========================================================================

outDir = tempdir();

% -------- (A) Long table: posture+mode (recommended for SLD code) --------
outCsvLong = fullfile(outDir, 'SLD_modes_long.csv');
fid = fopen(outCsvLong, 'w');

fprintf(fid, ['xpos_m;px_m;py_m;pz_m;mode;', ...
              'wn_rad_s;wd_rad_s;zeta;', ...
              'bX; bY; uX_1_kg; uY_1_kg; uXY_1_kg;', ...
              'kX_N_m;kY_N_m;cX_Ns_m;cY_Ns_m\n']);

for i = 1:N
    for k = 1:n
        fprintf(fid, ['%.12g;%.12g;%.12g;%.12g;%d;', ...
                      '%.12g;%.12g;%.12g;', ...
                      '%.12g;%.12g;%.12g;%.12g;%.12g;', ...
                      '%.12g;%.12g;%.12g;%.12g\n'], ...
            xpos(i), p_traj(i,1), p_traj(i,2), p_traj(i,3), k, ...
            wn(i,k), wd(i,k), zeta(i,k), ...
            PhiX_all(i,k), PhiY_all(i,k), uX(i,k), uY(i,k), uXY(i,k), ...
            k_eq_X(i,k), k_eq_Y(i,k), c_eq_X(i,k), c_eq_Y(i,k));
    end
end

fclose(fid);
fprintf('Saved (SLD long table): %s\n', outCsvLong);

% -------- (B) Summary table: dominant X/Y per posture (for your report) --------
outCsvSummary = fullfile(outDir, 'SLD_summary_XY.csv');
fid = fopen(outCsvSummary, 'w');

fprintf(fid, ['xpos_m;px_m;py_m;pz_m;', ...
              'DomModeX;wnX_rad_s;zetaX;kX_N_m;', ...
              'DomModeY;wnY_rad_s;zetaY;kY_N_m\n']);

for i = 1:N
    fprintf(fid, ['%.12g;%.12g;%.12g;%.12g;', ...
                  '%d;%.12g;%.12g;%.12g;', ...
                  '%d;%.12g;%.12g;%.12g\n'], ...
        xpos(i), p_traj(i,1), p_traj(i,2), p_traj(i,3), ...
        DomModeX(i), wn_domX(i), zeta_domX(i), k_modal_X(i), ...
        DomModeY(i), wn_domY(i), zeta_domY(i), k_modal_Y(i));
end

fclose(fid);
fprintf('Saved (summary): %s\n', outCsvSummary);


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

% ========================================================================
% Local functions
%   The following helper functions implement kinematics and dynamics and
%   are written to be consistent with the DH transform dh_num().
% ========================================================================


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

function Angcoord = GetAngles(a, d, alpha, p, Rot)
% IK numéric consistent with dh_num + jacobian_numeric
%solve position + constant orientation for Rot for each point p(i,:)

    P = p;
    phi = Rot(3); theta = Rot(2); psi = Rot(1);

    % wished rotation (tu convención)
    Rz = [cos(phi),-sin(phi),0; sin(phi),cos(phi),0; 0,0,1];
    Ry = [cos(theta),0,sin(theta); 0,1,0; -sin(theta),0,cos(theta)];
    Rx = [1,0,0; 0,cos(psi),-sin(psi); 0,sin(psi),cos(psi)];
    Rdes = Rz*Ry*Rx;

    N = rows(P);
    Angcoord = zeros(N,6);


    q = zeros(6,1);

    for i=1:N
        pdes = P(i,:).';
        [q,info] = ik_dls(q, pdes, Rdes, a, d, alpha);

        Angcoord(i,:) = q.';

        if info.fail
            printf("WARN IK: point %d not converged. pos_err=%.3e ori_err=%.3e\n", ...
                i, info.pos_err, info.ori_err);
        end
    end
end


function [q,info] = ik_dls(q0, pdes, Rdes, a, d, alpha)
% Damped Least Squares IK usando Jacobiano geométrico (espacial) de tu DH

    q = q0;

    maxIter = 200;
    tol_p   = 1e-6;     % [m]
    tol_o   = 1e-6;     % [rad]
    lambda  = 1e-2;     % damping
    stepMax = 0.2;      % [rad] límite por iter (evita saltos)

    fail = false;

    for it=1:maxIter
        [pfk, Rfk] = fk_pos_rot(q, a, d, alpha);

        ep = pdes - pfk;                 % position error (base)
        % orientation error when usign log(Rfk' * Rdes) and transfer to base
        Rerr_body = Rfk.' * Rdes;
        eo_body   = rotvec_from_R(Rerr_body);  % in end effector frame
        eo        = Rfk * eo_body;             % to base frame (coherente con Jacobiano espacial)

        if norm(ep) < tol_p && norm(eo) < tol_o
            break;
        end

        e = [ep; eo];

        J = jacobian_numeric(q, a, d, alpha);   % 6x6 (tu función)

        % DLS: dq = (J'J + λ²I)^-1 J' e
        dq = (J.'*J + (lambda^2)*eye(6)) \ (J.'*e);

        % limit for stability
        ndq = norm(dq);
        if ndq > stepMax
            dq = dq * (stepMax/ndq);
        end

        q = q + dq;

        % wrap a [-pi,pi] to avoid problem with certain configuration
        q = atan2(sin(q), cos(q));
    end

    [pfk, Rfk] = fk_pos_rot(q, a, d, alpha);
    ep = pdes - pfk;
    Rerr_body = Rfk.' * Rdes;
    eo_body   = rotvec_from_R(Rerr_body);
    eo        = Rfk * eo_body;

    info.pos_err = norm(ep);
    info.ori_err = norm(eo);
    info.fail    = (info.pos_err > 1e-4 || info.ori_err > 1e-4);
end


function [p, R] = fk_pos_rot(q, a, d, alpha)
% FK consistent with dh_num

    T = eye(4);
    for j=1:6
        T = T * dh_num(a(j), alpha(j), d(j), q(j));
    end
    p = T(1:3,4);
    R = T(1:3,1:3);
end


function w = rotvec_from_R(R)
% Return rotation vector
% stable for small angles

    tr = trace(R);
    c = (tr - 1)/2;
    c = max(-1, min(1, c));
    ang = acos(c);

    if ang < 1e-12
        w = [0;0;0];
        return;
    end

    % axis = (1/(2sin(ang))) * vee(R - R')
    ax = (1/(2*sin(ang))) * [R(3,2)-R(2,3);
                             R(1,3)-R(3,1);
                             R(2,1)-R(1,2)];
    w = ang * ax;
end



