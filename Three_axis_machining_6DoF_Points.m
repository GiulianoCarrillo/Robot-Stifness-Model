pkg load symbolic
pkg load mapping
clear; clc;
warning('off', 'all');
%------------------Input TCP position and orientation---------------------------

%Note: Only points in XZ plane works:
P = [2000,0,2200;2000,0,2000;2000,0,1800;2000,0,1600;1000,0,2200;1000,0,2000;1000,0,1800;1000,0,1600];
phi = deg2rad(0); theta = deg2rad(90); psi = deg2rad(0);

%-------------------------Robot parameters--------------------------------------

# Link lengths:
a1 = 330; a2 = 1150; a3 = 115;
a4 = 0; a5 = 0; a6 = 0;

# Link offset
d1  = 645; d2 = 0; d3 = 0;
d4 = 1220; d5 = 0; d6 = 240;

# Link twist
alph1 = pi/2; alph2 = 0; alph3 = pi/2;
alph4 = pi/2; alph5 = -pi/2; alph6 = 0;

%-------------------------Rotation matrix of TCP--------------------------------

Rz = [cos(phi),-sin(phi),0;sin(phi),cos(phi),0;0,0,1];
Ry = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
Rx = [1,0,0;0,cos(psi),-sin(psi);0,sin(psi),cos(psi)];
R = Rz*Ry*Rx;

for i = 1:rows(P)
ox = P(i,1); oy = P(i,2); oz = P(i,3);
point = ["p" num2str(i)];
#-------------------------------------------------------------------------------
#------------STEP 0: Calculate wrist center position----------------------------
#-------------------------------------------------------------------------------

xc = ox-d6*R(1,3); yc = oy - d6*R(2,3); zc = oz - d6*R(3,3);

#-------------------------------------------------------------------------------
#------------STEP 1: Determine joint parameters for elbow manipulator:----------
#-------------------------------------------------------------------------------
# Geometric considerations:
L = sqrt(a3^2+d4^2);
r1 = sqrt(xc^2+yc^2)-a1;
r2 = zc-d1;
r3 = sqrt(r1^2+r2^2);
D1 = (L^2-r3^2-a2^2) / (-2*a2*r3);
D3 = (r3^2-L^2-a2^2) / (-2*a2*L);

# Helping angles:
phi1 = atan2(sqrt(1-D1^2),D1);
phi2 = atan2(r2,r1);
phi3 = atan2(sqrt(1-D3^2),D3);
delta = atan2(a3,d4); % Correction to get actual motor angle

# Elbow up:
q1 = atan2(yc,xc);
q2 = (phi2+phi1);
q3 = -(pi-phi3+delta);

#Convert to deg
q1_deg = -q1*(180/pi);
q2_deg = -q2*(180/pi);
q3_deg = -q3*(180/pi);

#Save in structure:
q.("q1") = q1_deg;
q.("q2") = q2_deg;
q.("q3") = q3_deg;

#-------------------------------------------------------------------------------
#-------------------------STEP 2: Evaluate matrix R03---------------------------
#-------------------------------------------------------------------------------

A1 = [cos(q1), 0, sin(q1), a1*cos(q1); sin(q1), 0, -cos(q1), a1*sin(q1); 0, 1, 0, d1; 0, 0, 0, 1];
A2 = [cos(q2), -sin(q2), 0, a2*cos(q2); sin(q2), cos(q2), 0, a2*sin(q2); 0, 0, 1, 0; 0, 0, 0, 1];
A3 = [cos(q3), 0, sin(q3), a3*cos(q3); sin(q3), 0, -cos(q3), a3*sin(q3); 0, 1, 0, 0; 0, 0, 0, 1];

T = A1*A2*A3;
r03 = T(1:3,1:3);
r03T = transpose(r03);
R36 = r03T; #Right hand side of equation should be: (R03)^T*R, but something is wrong here.

#-------------------------------------------------------------------------------
#----------------STEP 3: Solve R36 = (R03)^T*R (similar to euler angles)--------
#-------------------------------------------------------------------------------
# In this case we need to define the ranges of q4, q5 and q6:
eps = 1e-6;
cond1 = abs(R36(1,3)) > eps;
cond2 = abs(R36(2,3)) > eps;

# a) IF not BOTH R36(1,3) and R36(2,3) are 0 => sin(q5) neq 0
if or(cond1, cond2)
    q4 = zeros(1,2);
    q5 = zeros(1,2);
    q6 = zeros(1,2);

    q4(1,1) = atan2(R36(2,3),R36(1,3));
    q5(1,1) = atan2(sqrt(R36(1,3)^2+R36(2,3)^2),-R36(3,3));
    q6(1,1) = atan2(-R36(3,2),R36(3,1));
    q4(1,2) = q4(1,1)+pi;
    q5(1,2) = atan2(-sqrt(1-R36(3,3)^2),-R36(3,3));
    q6(1,2) = q6(1,1)+pi;

    q4_deg = q4*(180/pi);
    q5_deg = q5*(180/pi);
    q6_deg = q6*(180/pi);

    q.("q4") = q4_deg;
    q.("q5") = q5_deg;
    q.("q6") = q6_deg;

# b) IF R36(1,3) = R36(2,3) = 0 => R(3,3) = cos(5)= +/- 1
else
    # If R36(3,3) = -1 => cos(q5) = 1 and sin(q5) = 0:
  if R36(3,3) = 1
    q5 = 0;
    q4 = 0; #This can have infinitely many values - convention to choose =0
    q6 = atan2(R36(2,1),R36(1,1))

    # If R36(3,3) = 1 => cos(q5) = -1 and sin(q5) = 0:
  elseif R36(3,3) = -1
    q5 = pi;
    q4 = 0; #This can have infinitely many values - convention to choose =0
    q6 = -atan2(R36(2,1),R36(1,1))
  end
    q4_deg = q4*(180/pi);
    q5_deg = q5*(180/pi);
    q6_deg = q6*(180/pi);

    q.("q4") = q4_deg;
    q.("q5") = q5_deg;
    q.("q6") = q6_deg;

end

Q.(point) = q;

end
