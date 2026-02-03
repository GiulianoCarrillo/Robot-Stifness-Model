function [ Angcoord ] = GetAngles(a, d, alpha, p, Rot)
% =========================================================================
% GetAngles
%
%   Calculate joint parameters for a given position or set of positions
%
% INPUTS:
%   a,d,alpha : DH parameters [
%   p       : position of TCP frame origin [x1,y1,z1; ... ; xn,yn,zn]
%   Rot     : Orientation of TCP frame [Rx,Ry,Rz] in rad
%
% OUTPUTS:
%   AngCoord : Joint angles for each position in p
%   Pwrist : Position of wrist center for each position in p
%
% =========================================================================

%--------------------------------Initialize-------------------------------------
ai = a;
di = d;
alphi = alpha;
P = p;
phi = Rot(3); theta = Rot(2); psi = Rot(1);

%-------------------------Rotation matrix of TCP--------------------------------

Rz = [cos(phi),-sin(phi),0;sin(phi),cos(phi),0;0,0,1];
Ry = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
Rx = [1,0,0;0,cos(psi),-sin(psi);0,sin(psi),cos(psi)];
R = Rz*Ry*Rx;

%-------------------------Calculate joint angles--------------------------------

for i = 1:rows(P)
ox = P(i,1); oy = P(i,2); oz = P(i,3);
#-------------------------------------------------------------------------------
#------------STEP 0: Calculate wrist center position----------------------------
#-------------------------------------------------------------------------------

xc = ox-di(6)*R(1,3); yc = oy - di(6)*R(2,3); zc = oz - di(6)*R(3,3);
Pc(i,1) = xc; Pc(i,2) = yc; Pc(i,3) = zc;

#-------------------------------------------------------------------------------
#------------STEP 1: Determine joint parameters for elbow manipulator:----------
#-------------------------------------------------------------------------------
# Geometric considerations:
L = sqrt(ai(3)^2+di(4)^2);
r1 = sqrt(xc^2+yc^2)-ai(1);
r2 = zc-di(1);
r3 = sqrt(r1^2+r2^2);
D1 = (L^2-r3^2-ai(2)^2) / (-2*ai(2)*r3);
D3 = (r3^2-L^2-ai(2)^2) / (-2*ai(2)*L);

# Helping angles:
phi1 = atan2(sqrt(1-D1^2),D1);
phi2 = atan2(r2,r1);
phi3 = atan2(sqrt(1-D3^2),D3);
delta = atan2(ai(3),di(4)); % Correction to get actual motor angle

# Elbow up:
q1 = -atan2(yc,xc); q2 = -(phi2+phi1); q3 = (pi-phi3+delta);

#Save in matrix:
Q(i,1) = q1; Q(i,2) = q2; Q(i,3) = q3;

#-------------------------------------------------------------------------------
#-------------------------STEP 2: Evaluate matrix R03---------------------------
#-------------------------------------------------------------------------------

A1 = [cos(q1), 0, sin(q1); sin(q1), 0, -cos(q1); 0, 1, 0,];
A2 = [cos(q2), -sin(q2), 0; sin(q2), cos(q2), 0; 0, 0, 1];
A3 = [cos(q3), 0, sin(q3); sin(q3), 0, -cos(q3); 0, 1, 0];

r03 = A1*A2*A3;
r03T = transpose(r03);
R36 = r03T; #Right hand side of equation should be: (R03)^T*R, but something is wrong here.

#-------------------------------------------------------------------------------
#----------------STEP 3: Solve R36 = (R03)^T*R (similar to euler angles)--------
#-------------------------------------------------------------------------------
# In this case we need to define the ranges of q4, q5 and q6:
eps = 1e-6;
cond1 = (abs(R36(1,3)) > eps);
cond2 = (abs(R36(2,3)) > eps);

# a) IF not BOTH R36(1,3) and R36(2,3) are 0 => sin(q5) neq 0
  if or(cond1, cond2)
      q5(1,1) = atan2(sqrt(R36(1,3)^2+R36(2,3)^2),-R36(3,3));
      q4(1,1) = atan2(R36(2,3),R36(1,3));
      q5(1,1) = atan2(-sqrt(R36(1,3)^2+R36(2,3)^2),-R36(3,3));
      q6(1,1) = atan2(-R36(3,2),R36(3,1));
      q4(1,2) = q4(1,1)+pi;
      q5(1,2) = atan2(sqrt(1-R36(3,3)^2),-R36(3,3));
      q6(1,2) = q6(1,1)+pi;

      %Always choose q4 = q6 = 0 or 2pi (and not pi/-pi) and make sure -125 < q5 < 125
      checkq4 = (abs( abs(q4)-pi ) > 0.1);
      checkq5 = (abs(q5) < 125*pi/180);
      if checkq5
        id = find(checkq4);
        Q(i,4) = q4(1,id); Q(i,5) = q5(1,id); Q(i,6) = q6(1,id);
      else
        point = ["p" num2str(i)];
        printf("For point %d, q5 = %d or %d rad is out of range\n", [i,q5]);
        Q(i,:) = [];
      end

  # b) IF R36(1,3) = R36(2,3) = 0 => R(3,3) = cos(5)= +/- 1
  else
      # If R36(3,3) = -1 => cos(q5) = 1 and sin(q5) = 0:
    if R36(3,3) = 1
      q5 = 0;
      q4 = 0; #This can have infinitely many values - convention to choose =0
      q6 = atan2(R36(2,1),R36(1,1))

      Q(i,4) = q4; Q(i,5) = q5; Q(i,6) = q6;

      # If R36(3,3) = 1 => cos(q5) = -1 and sin(q5) = 0:
    elseif R36(3,3) = -1
      q5 = pi;
      q4 = 0; #This can have infinitely many values - convention to choose =0
      q6 = -atan2(R36(2,1),R36(1,1))

      Q(i,4) = q4; Q(i,5) = q5; Q(i,6) = q6;
    end

  end

end
Angcoord = Q;
Pwrist = Pc;
endfunction

