function plotFRF(wn,zeta,m,xpos);
% =========================================================================
% plotFRF
%
%   Plot the magnitude of the FRF vs. frequency [rad/s] and position [m]
%   FRF: g(s) = 1 / ((k/wn^2)*s^2+2*zeta*(k/wn)*s+k);
%   with s = jw, k = wn^2*m;
%
% INPUTS:
%   wn      : Natural frequencies
%   zeta    : damping factor
%   m       : masses
%   xpos    : positions
%
% OUTPUTS:
%   plot of FRF magnitude
%
% REQUIRED PACKAGES
%   control (type "pkg load control" in octave prompt
%
% =========================================================================

% Define transfer function variable:
s = tf('s');

% Define freq. upper limit:
ymax = ceil(max(max(wn))/100)*100+100;

% Calculate magnitude of FRF (total and individual modes):
for i = 1:length(xpos)
    k1 = wn(i,1)^2*m(1);
    k2 = wn(i,2)^2*m(2);
    k3 = wn(i,3)^2*m(3);
    k4 = wn(i,4)^2*m(4);
    k5 = wn(i,5)^2*m(5);
    k6 = wn(i,6)^2*m(6);
    g1 = 1 / ((k1/wn(i,1)^2)*s^2+2*zeta(i,1)*(k1/wn(i,1))*s+k1);
    g2 = 1 / ((k2/wn(i,2)^2)*s^2+2*zeta(i,2)*(k2/wn(i,2))*s+k2);
    g3 = 1 / ((k3/wn(i,3)^2)*s^2+2*zeta(i,3)*(k3/wn(i,3))*s+k3);
    g4 = 1 / ((k4/wn(i,4)^2)*s^2+2*zeta(i,4)*(k4/wn(i,4))*s+k4);
    g5 = 1 / ((k5/wn(i,5)^2)*s^2+2*zeta(i,5)*(k5/wn(i,5))*s+k5);
    g6 = 1 / ((k6/wn(i,6)^2)*s^2+2*zeta(i,6)*(k6/wn(i,6))*s+k6);
    g_tot = g1+g2+g3+g4+g5+g6;

    [mag1(:,i),pha1(:,i),w1(:,i)] = bode(g1);
    [mag2(:,i),pha2(:,i),w2(:,i)] = bode(g2);
    [mag3(:,i),pha3(:,i),w3(:,i)] = bode(g3);
    [mag4(:,i),pha4(:,i),w4(:,i)] = bode(g4);
    [mag5(:,i),pha5(:,i),w5(:,i)] = bode(g5);
    [mag6(:,i),pha6(:,i),w6(:,i)] = bode(g6);
    [mag_tot(:,i),pha_tot(:,i),w_tot(:,i)] = bode(g_tot);
    num = length(w_tot(:,i));
    num2 = length(w1);

    x1(:,i) = ones(1,num2)*xpos(i);
    x2(:,i) = ones(1,num2)*xpos(i);
    x3(:,i) = ones(1,num2)*xpos(i);
    x4(:,i) = ones(1,num2)*xpos(i);
    x5(:,i) = ones(1,num2)*xpos(i);
    x6(:,i) = ones(1,num2)*xpos(i);
    xtot(:,i) = ones(1,num)*xpos(i);
end

% =========================================================================
% Plot magnitude of total transfer fucntion H = sum(H_i):
% =========================================================================

% Full plot:
figure; hold on
plot3(xtot(:,1),w_tot(:,1),mag_tot(:,1),'color','blue');
plot3(xtot(:,2),w_tot(:,2),mag_tot(:,2),'color','blue');
plot3(xtot(:,3),w_tot(:,3),mag_tot(:,3),'color','blue');
plot3(xtot(:,4),w_tot(:,4),mag_tot(:,4),'color','blue');
plot3(xtot(:,5),w_tot(:,5),mag_tot(:,5),'color','blue');
plot3(xtot(:,6),w_tot(:,6),mag_tot(:,6),'color','blue');
plot3(xtot(:,7),w_tot(:,7),mag_tot(:,7),'color','blue');
plot3(xtot(:,8),w_tot(:,8),mag_tot(:,8),'color','blue');
plot3(xtot(:,9),w_tot(:,9),mag_tot(:,9),'color','blue');
ylim([0 ymax]);
%zlim([0 0.0005]);
xlabel('x [m]');
ylabel('\omega [rad/s]');
zlabel('Magnitude dB');
title('Magnitude of combined FRF');
set(gca, 'YScale', 'log'); % Log scale for frequency
view(3500,200);

% Plot zoom-in version
figure; hold on
plot3(xtot(:,1),w_tot(:,1),mag_tot(:,1),'color','blue');
plot3(xtot(:,2),w_tot(:,2),mag_tot(:,2),'color','blue');
plot3(xtot(:,3),w_tot(:,3),mag_tot(:,3),'color','blue');
plot3(xtot(:,4),w_tot(:,4),mag_tot(:,4),'color','blue');
plot3(xtot(:,5),w_tot(:,5),mag_tot(:,5),'color','blue');
plot3(xtot(:,6),w_tot(:,6),mag_tot(:,6),'color','blue');
plot3(xtot(:,7),w_tot(:,7),mag_tot(:,7),'color','blue');
plot3(xtot(:,8),w_tot(:,8),mag_tot(:,8),'color','blue');
plot3(xtot(:,9),w_tot(:,9),mag_tot(:,9),'color','blue');
ylim([0 ymax]);
zlim([0 0.000003]);
xlabel('x [m]');
ylabel('\omega [rad/s]');
zlabel('Magnitude dB');
title('Magnitude of combined FRF - zoom');
view(3500,200);
set(gca, 'YScale', 'log'); % Log scale for frequency

% =========================================================================
% Plot magnitude of FRF for individual modes
% =========================================================================

% Plot as is:
figure; hold on
for i = 1:length(xpos)
  plot3(x1(:,i),w1(:,i),mag1(:,i),'color','blue');
  plot3(x2(:,i),w2(:,i),mag2(:,i),'color','red');
  plot3(x3(:,i),w3(:,i),mag3(:,i),'color','black');
  plot3(x4(:,i),w4(:,i),mag4(:,i),'color','magenta');
  plot3(x5(:,i),w5(:,i),mag5(:,i),'color','cyan');
  plot3(x6(:,i),w6(:,i),mag6(:,i),'color','yellow');
end
ylim([0 ymax]);
xlabel('x [m]');
ylabel('\omega [rad/s]');
zlabel('Magnitude dB');
view(100,40);
%set(gca, 'YScale', 'log'); % Log scale for frequency
title('FRF - individual modes');
legend('Mode 1','Mode 2',' Mode 3', 'Mode 4','Mode 5','Mode 6','location','eastoutside');


%Plot normalized version
figure; hold on
for i = 1:length(xpos)
  norm1 = mag1(:,i)./max(mag1(:,i));
  norm2 = mag2(:,2)./max(mag2(:,2));
  norm3 = mag3(:,3)./max(mag3(:,3));
  norm4 = mag4(:,4)./max(mag4(:,4));
  norm5 = mag5(:,5)./max(mag5(:,5));
  norm6 = mag6(:,6)./max(mag6(:,6));
  plot3(x1(:,i),w1(:,i),norm1,'color','blue');
  plot3(x2(:,i),w2(:,i),norm2,'color','red');
  plot3(x3(:,i),w3(:,i),norm3,'color','black');
  plot3(x4(:,i),w4(:,i),norm4,'color','magenta');
  plot3(x5(:,i),w5(:,i),norm5,'color','cyan');
  plot3(x6(:,i),w6(:,i),norm6,'color','yellow');
end
ylim([0 ymax]);
xlabel('x [m]');
ylabel('\omega [rad/s]');
zlabel('Normalized Magnitude');
view(100,40);
title('FRF - individual modes, normalized');
legend('Mode 1','Mode 2',' Mode 3', 'Mode 4','Mode 5','Mode 6','location','eastoutside');
%set(gca, 'YScale', 'log'); % Log scale for frequency
end
