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

pw = [ox; oy; oz] - R * [0; 0; di(6)];
xc = pw(1);
yc = pw(2);
zc = pw(3);
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

#-------------------------------------------------------------------------------
#-------------------------STEP 2: Evaluate matrix R03 (DH-consistent)-----------
#-------------------------------------------------------------------------------

T03 = dh_num(ai(1), alphi(1), di(1), q1) * ...
      dh_num(ai(2), alphi(2), di(2), q2) * ...
      dh_num(ai(3), alphi(3), di(3), q3);

R03  = T03(1:3,1:3);
R36  = R03' * R;   % <-- ESTA es la ecuación correcta

#-------------------------------------------------------------------------------
#----------------STEP 3: Solve wrist angles from R36----------------------------
#-------------------------------------------------------------------------------

eps = 1e-9;
s = sqrt(R36(1,3)^2 + R36(2,3)^2);

if s > eps
    % Two solutions for the wrist (q5 positive/negative)
    q4a = atan2(R36(2,3), R36(1,3));
    q5a = atan2( s, -R36(3,3));
    q6a = atan2(-R36(3,2), R36(3,1));

    q4b = q4a + pi;
    q5b = atan2(-s, -R36(3,3));
    q6b = q6a + pi;

    % Choose one solution (keep your criteria, but now consistent)
    cand = [q4a q5a q6a;
            q4b q5b q6b];

    % example criterion: keep |q5| < 125deg and avoid q4 near +-pi if possible
    ok = (abs(cand(:,2)) < 125*pi/180);

    if any(ok)
        % if both ok, pick the one with q4 farthest from +-pi
        idx_ok = find(ok);
        [~,best] = max(abs(abs(cand(idx_ok,1)) - pi));
        pick = idx_ok(best);
        Q(i,4) = cand(pick,1);
        Q(i,5) = cand(pick,2);
        Q(i,6) = cand(pick,3);
    else
        printf("Point %d: wrist solutions out of q5 range\n", i);
        Q(i,4:6) = [NaN NaN NaN];
    end

else
    % Singularity: q5 ~ 0 or pi
    if R36(3,3) > 0
        q5 = 0;
        q4 = 0;
        q6 = atan2(R36(2,1), R36(1,1));
    else
        q5 = pi;
        q4 = 0;
        q6 = -atan2(R36(2,1), R36(1,1));
    end
    Q(i,4) = q4; Q(i,5) = q5; Q(i,6) = q6;
end
  end

end
Angcoord = Q;
Pwrist = Pc;
endfunction

