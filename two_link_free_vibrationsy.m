function two_link_free_vibrationsy()
pkg load symbolic
clear; clc;
warning('off', 'all');
%% =============================================================================
%% ------------------ Parameters -----------------------------------------------
%% =============================================================================

%% ------------------ Evaluation range -----------------------------------------

dp = 0.1;           % Interval of evaluation points.
p0 = [0.25,0.0];    % Initial Point
path = [0.5,0.0]; % Final point or path ([100,100; 100,140; 140,180])


%% ------------------ Physical parameters --------------------------------------
##syms t q1(t) q2(t)

% Geometry
l1  = 0.380;       % [m] length link 1
l2  = 0.380;        % [m] length link 2

q1_min = 0;  q1_max = 180;             % joint limits
q2_min = -180; q2_max = 350;



lc1 = 0.196;      % [m] COM link 1 from joint 1
lc2 = 0.176;      % [m] COM link 2 from joint 2

%q1 = deg2rad(20) ; % (Angle of first joint)
%q2 = deg2rad(30) ; % (Angle of second joint)

% Mass & inertia about COM (z out of plane)
%m1 = 2.845;        % [kg]
%m2 = 2.774;        % [kg]

m1 = 5.639;        % [kg]
m2 = 2.774;        % [kg]

%Iz1 = 0.073;      % [kg m^2] Inertia IZ taken on the reference system with origin in the CM of the link.
%Iz2 = 0.03460;      % [kg m^2] Inertia IZ taken on the reference system with origin in the CM of the link.

Iz1 = 0.0563;      % [kg m^2] Inertia IZ taken on the reference system with origin in the CM of the link.
Iz2 = 0.0147;      % [kg m^2] Inertia IZ taken on the reference system with origin in the CM of the link.

% Gravity
g  = 9.81;       % [m/s^2]

% Torsional springs (potential Energy: 0.5*(q - q0)'*Ks*(q - q0)).
k1 = 150;      % [N·m/rad]
k2 = 1100;      % [N·m/rad]

% Viscous damping.
c1 = 20 ;       % [N·m·s/rad]; set >0 if we want damping
c2 = 20;       % [N·m·s/rad]; set >0 if we want damping

%% =============================================================================
%% ------------------ Evaluation -----------------------------------------------
%% =============================================================================

[Angcord] = joint_calculator(l1,l2,p0,path,dp,q1_min,q2_min,q1_max,q2_max)
[Cs,Ks,D,Cc,G,kg] = Dynamic_Model();

%% =============================================================================
%% =============================================================================
%% ------------------ Functions -----------------------------------------------
%% =============================================================================

%% =============================================================================
%% 2R PLANAR ARM - Linearization.
%% =============================================================================
n = size(Angcord, 1);   % o n = rows(Angcord);

wn = zeros(n,2) ;

for i = 1:n;

  q1=deg2rad(Angcord(i,1)) ;
  q2=deg2rad(Angcord(i,2)) ;
  dq1=0 ;
  dq2=0 ;

  D_eval = double(subs(D,{"l1","l2","lc1","lc2","m1","m2","Iz1","Iz2","g","q1","q2"},{l1, l2, lc1, lc2, m1, m2, Iz1, Iz2, g, q1, q2}));
  G_eval = double(subs(G,{"l1","l2","lc1","lc2","m1","m2","Iz1","Iz2","g","q1","q2"},{l1, l2, lc1, lc2, m1, m2, Iz1, Iz2, g, q1, q2}));
  Cc_eval = double(subs(Cc,{"l1","l2","lc1","lc2","m1","m2","Iz1","Iz2","g","q1","q2","dq1","dq2"},{l1, l2, lc1, lc2, m1, m2, Iz1, Iz2, g, q1, q2, dq1, dq2}));
  Cs_eval = double(subs(Cs, {"c1", "c2"}, {c1, c2}));
  Ks_eval = double(subs(Ks, {"k1", "k2"}, {k1, k2}));
  kg_eval = double(subs(kg,{"l1","l2","lc1","lc2","m1","m2","Iz1","Iz2","g","q1","q2"},{l1, l2, lc1, lc2, m1, m2, Iz1, Iz2, g, q1, q2}));

  Kt= Ks_eval + kg_eval;

%% =============================================================================
%% Modal Parameters.
%% =============================================================================

  [eigVec, eigVal] = eig(Kt, D_eval);
  wn(i,:) = sqrt(diag(eigVal))   % rad/s

  % Modos propios no amortiguados
  phi = eigVec;    % cada columna es un modo

  % Normalizar cada modo para tener máx = 1 (solo para interpretación)
  for k = 1:2
    phi(:,k) = phi(:,k) / max(abs(phi(:,k)));
  end
  phi


%% =============================================================================
%% State-space formulation
%% =============================================================================

invD = inv(D_eval);

A11 = zeros(2);
A12 = eye(2);
A21 = -invD * Kt;
A22 = -invD * Cs_eval;

A = [A11 A12; A21 A22];

lambda = eig(A);

% Dampened frecuencies (rad/s)
wd = abs(imag(lambda))

% Modal dampening
zeta = -real(lambda) ./ sqrt(real(lambda).^2 + imag(lambda).^2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------- FRF en el TCP para configuración i --------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numerical Jacobian of TCP (for FRF in X-Y-Z)
% You only need the first 2 rows of the linear velocity Jacobian
J_eval = [ -l1*sin(q1) - l2*sin(q1+q2),   -l2*sin(q1+q2);
            l1*cos(q1) + l2*cos(q1+q2),    l2*cos(q1+q2) ];

% Frecuency range
w_min = 1;
w_max = 300;
dw = 1;
w = w_min:dw:w_max;
nw = length(w);

% Preallocate FRF in TCP (magnitude of Hxx → address X)
FRF_TCP_mag = zeros(1, nw);

for ii = 1:nw
    ww = w(ii);

    % FRF articular (2x2)
    H = inv(-ww^2 * D_eval + 1i*ww * Cs_eval + Kt);

    % FRF en TCP:  H_tcp = J * H * J'
    H_tcp = J_eval * H * J_eval';

    % Tomamos la FRF en la dirección X (Hxx)
    FRF_TCP_mag(ii) = abs(H_tcp(1,1));
endfor

% Save for graphing (3D SLD)
FRF_matrix(i, :) = FRF_TCP_mag;


end

%% =============================================================================
%% Chart Natural Frequencies
%% =============================================================================

% X coordinate
X = Angcord(:,3);

% Natural Frequencies
wn1 = wn(:,1);
wn2 = wn(:,2);

figure;
plot(X, wn1, 'b-', 'LineWidth', 2); hold on;
plot(X, wn2, 'r-', 'LineWidth', 2);

xlabel('X position [m]');
ylabel('Natural frequencies \omega_n [rad/s]');
title('Natural frequencies along trajectory');
legend('Mode 1', 'Mode 2');
grid on;

n;

end
function [Angcord] = joint_calculator(l1,l2,p0,path,dp,q1_min,q2_min,q1_max,q2_max)
%% =============================================================================
%% Joint calculator - Inverse Kinematic.
%% =============================================================================


%% ------------------ Trayectory -----------------------------------------------


L = sqrt(l1^2 + l2^2);                 % Max lenght

points = [p0; path];                    % Extreme points
all_points = [];                        % Empty matrix for all the points

for i = 1:(rows(points)-1)              % subdivision for each link of the path

    p1 = points(i,:);
    p2 = points(i+1,:);
    diff = p2 - p1;
    dist = norm(diff);                  % Distance of the path
    nn = max(1, ceil(dist/dp));         % Ensure that there is at least one subdivision
    t = linspace(0,1,nn+1)';             % Linspace(a, b, m) generates m numbers equally spaced between a and b (inclusive).
    seg = p1 + t .* diff;
    if i > 1
        seg = seg(2:end,:);             % avoid to duplicate points
    end
    all_points = [all_points; seg];
end

% --------------------- Angle joint calculation (just elbow-up) ----------------


num_points = rows(all_points);
anglesUp = zeros(num_points, 2);  % [q1 q2] in degrees

for i = 1:num_points
    x = all_points(i,1);
    z = all_points(i,2);
    r = sqrt(x^2 + z^2);

    if r <= L
        r2 = x^2 + z^2;
        cos_q2 = (r2 - l1^2 - l2^2) / (2*l1*l2);

        if abs(cos_q2) <= 1
            q2_rad = -acos(cos_q2);  % “elbow-up”
            q1_rad = atan2(z,x) - atan2(l2*sin(q2_rad), l1 + l2*cos(q2_rad));
            q1_deg = rad2deg(q1_rad);
            q2_deg = rad2deg(q2_rad);


            if q1_deg > q1_min && q1_deg < q1_max && ...
               q2_deg > q2_min && q2_deg < q2_max
                anglesUp(i,:) = [q1_deg, q2_deg];
            else
                anglesUp(i,:) = [NaN, NaN]; %
                disp("out of angle range.")
            end
        else
            anglesUp(i,:) = [NaN, NaN];
            disp("out of reach.")
        end
    else
        anglesUp(i,:) = [NaN, NaN];
        disp("out of max radius")
    end
end
% --- Resultados ---
anglesUp;
all_points;
Angcord = [anglesUp, all_points];
end
function [Cs,Ks,D,Cc,G,kg] = Dynamic_Model () ;

%% =============================================================================
%% 2R PLANAR ARM - Mathematical model.
%% =============================================================================

%%------------------------------------------------------------------------------
%% ------------------ Parameters -----------------------------------------------
%%------------------------------------------------------------------------------
syms t q1 q2 l1 l2 lc1 lc2 m1 m2 Iz1 Iz2 g k1 k2 c1 c2 dq1 dq2 ;


q = [q1; q2];
dq = [dq1; dq2];

n = length(q); % Number of generalized coordinates


%%------------------------------------------------------------------------------
%% ------------------ Reference systems ----------------------------------------
%%------------------------------------------------------------------------------


% ------------------ Inertial system -------------------------------------------
x0 = sym([1,0,0]) ;
y0 = sym([0,1,0]) ;
z0 = sym([0,0,1]) ;
O0 = sym([0,0,0]) ;

Rs0 = sym([x0;y0;z0]) ;

% ------------------ Joint 1 ---------------------------------------------------

% Reference system of link 1 is located on the center of joint 2 for convention.

O1 = [ l1*cos(q1) , l1*sin(q1) , sym(0) ] ; % Origin Joint 1
x1 = sym([1,0,0]) ;
y1 = sym([0,1,0]) ;
z1 = sym([0,0,1]) ;

Rs1 = [x1;y1;z1] ;


R10 = [ sym(cos(q1)), sym(- sin(q1)), sym(0) ; sym(sin(q1)), sym(cos(q1)), sym(0) ; sym(0),sym(0), sym(1) ];  % Rotational Matrix from RS1 TO 0 / P0 = R10*P1

% ------------------ Joint 2 ---------------------------------------------------

O2 = [ l1*cos(q1) + l2*cos(q2+q1) , l1*sin(q1) + l2*sin(q1+q2) , sym(0) ] ; % Origin Joint 2
x2 = sym([1,0,0]) ;
y2 = sym([0,1,0]) ;
z2 = sym([0,0,1]) ;

Rs2 = [x2;y2;z2] ;

R21 = [ sym(cos(q2)), sym(-sin(q2)), sym(0) ; sym(sin(q2)), sym(cos(q2)), sym(0) ; sym(0),sym(0), sym(1) ] ;

R20 = R10 * R21 ;

% ------------------ Cm1 -------------------------------------------------------

Om1 = [ lc1*cos(q1) , lc1*sin(q1) , sym(0) ] ; % Origin Cm1
xm1 = sym([1,0,0]) ;
ym1 = sym([0,1,0]) ;
zm1 = sym([0,0,1]) ;

Rsm1 = [xm1;ym1;zm1] ;

Rm110 = R10  ;% Rotational Matrix of CM1 to RS0

% ------------------ Cm2 -------------------------------------------------------

Om2 = [ l1*cos(q1) + lc2*cos(q2+q1) , l1*sin(q1) + lc2*sin(q1+q2) , sym(0) ] ; % Origin Cm2
xm2 = sym([1,0,0]) ;
ym2 = sym([0,1,0]) ;
zm2 = sym([0,0,1]) ;

Rsm2 = [xm2;ym2;zm2] ;

Rm220 = R20 ;

%%------------------------------------------------------------------------------
%% ------------------ Jacobians ------------------------------------------------
%%------------------------------------------------------------------------------

% The velocity of each point on the robot can be written as a function of the generalized coordinates and the geometric characteristics of the robot.

% J = [Ji, .... , Jn] ----> Ji= [Zi-1 x (On - Oi-1) ; Zi-1] with i as the number of the joint.

% Jacobian of center of mass 1: Jcm1 = [Jcm11  Jcm12]
% Jacobian of center of mass 2: Jcm2 = [Jcm21  Jcm22]

% where the first rows are tangencial velocity and the seocond ones returns the angular velocity

Jcm11 = [cross(Rs0(3,:),(Om1-O0)).'; Rs0(3,:).'] ;

Jcm12 = sym([0;0;0;0;0;0]) ; % there is not influence of the joint 2 in the movement of the point of cm1.

Jcm21 = [cross(Rs0(3,:),(Om2-O0)).';Rs0(3,:).'] ;
Jcm22 = [cross(Rs1(3,:),(Om2-O1)).';Rs1(3,:).'] ;

Jcm1= [Jcm11 , Jcm12] ;
Jcm2= [Jcm21 , Jcm22] ;

%%------------------------------------------------------------------------------
%% ------------------ Modal Matrix ---------------------------------------------
%%------------------------------------------------------------------------------

%% ------------------ Damping Matrix -------------------------------------------

Cs = diag([c1 c2]);

%% ------------------ Stifness Matrix ------------------------------------------

Ks = diag([k1 k2]);

% q0 = [q1; q2];   % set q0 = q* if we want neutral springs at posture q*

%% ------------------ Mass Matrix ----------------------------------------------

M= [m1 m2];

%% ------------------ Rotational Inertia Matrix --------------------------------

% By writing Ii over the principal axes of inertia of the link, the rotational inertia matrix around the z axis is:

I1 = [sym(0),sym(0),sym(0);sym(0),sym(0),sym(0);sym(0),sym(0),sym(Iz1)] ;

I2 = [sym(0),sym(0),sym(0);sym(0),sym(0),sym(0);sym(0),sym(0),sym(Iz2)] ;

%%------------------------------------------------------------------------------
%% ------------------ Matrix Valued functions-----------------------------------
%%------------------------------------------------------------------------------

%------------------ Inertia matrix D(q) ----------------------------------------

% D(q) = sum_i to n [ m_i * Jv_i^T * Jv_i + Jw_i^T * Ri0 * I_i * Ri0^T * Jw_i ]

D = m1*Jcm1(1:3,:).'*Jcm1(1:3,:)+Jcm1(4:6,:).'*Rm110*I1*Rm110.'*Jcm1(4:6,:)+m2*Jcm2(1:3,:).'*Jcm2(1:3, :)+Jcm2(4:6,:).'*Rm220*I2*Rm220.'*Jcm2(4:6,:)

%------------------ Coroilis Matriz Cc(q;qdot) ---------------------------------

Cc = sym(zeros(n,n));  % o sym(zeros(...)) si querés simbólico


for k = 1:n ;
  for j = 1:n ;
    for i = 1:n ;
     Cc(k,j) = Cc(k,j) + (1/2)*(diff(D(k,j), q(i))+ diff(D(k,i), q(j))-diff(D(i,j), q(k)))*dq(i);
    end ;
  end ;
end ;
Cc;


%------------------ Gravity Matrix G ---------------------------------

Pg1= g*M(1)*Om1(2) ;
Pg2= g*M(2)*Om2(2) ;
Pg = Pg1 + Pg2 ;

G = sym(zeros(n,1)) ;

for k = 1:n ;
  G(k)= G(k) + diff(Pg,q(k)) ;
end ;

kg = sym(zeros(n,n)) ;

for k = 1:n ;
  for j = 1:n ;
  kg(k,j)= kg(k,j) + diff(G(k),q(j)) ;
  end ;
end ;

end
