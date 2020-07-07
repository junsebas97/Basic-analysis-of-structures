%clc, clear all , close all
%{
This codes assess the stress and strain distribution in the cross section
of a circular or tubular bar subjected to constant torsion

MADE BY: Juan Sebastián Delgado Trujillo
%}

%% problem parameters:

T     = 94;          % [kN - m]
G     = 84e6;        % [kPa]
L     = 1.5;         % [m]
r_int = 2*0.0254;    % [m]
r_ext = 3*0.0254;    % [m]

% solid bar trial:
r_ext = 0.072;       r_int = 0;

%% main:

Ip = pi*(r_ext^4 - r_int^4)/2;    % polar inertia

phi   = (T*L)/(Ip*G);    % rotation of the final section  [rad]
theta =  phi/L;          % rotation angle per length unit [rad/m]

% equation of the stress and strain distribution in the section
tau_ro = @(ro) T*ro./Ip;
gam_ro = @(ro) tau_ro(ro)./G;

% to find the maxima and the minima they are evaluated in the surfaces
tau_max = tau_ro(r_ext);        gam_max = gam_ro(r_ext);
tau_min = tau_ro(r_int);        gam_min = gam_ro(r_int);

%% plots ans reports:

fprintf('\nrotation of the final section = %.3f°', phi*180/pi)
fprintf('\ntau_max = %.3f MPa', tau_max/1000)
fprintf('\ntau_min = %.3f MPa', tau_min/1000)
fprintf('\ngamma_max = %.3f',   gam_max)
fprintf('\ngamma_min = %.3f\n', gam_min)

% this block is for creating the circular domain
nR = linspace(r_int, r_ext) ;        nT = linspace(0,     2*pi);
[R, T] = meshgrid(nR,nT) ;
X = R.*cos(T);                       Y = R.*sin(T);

figure
hold on
plot(nR, tau_ro(nR), 'b-', -nR, tau_ro(nR), 'b-')
axis tight
grid on
title('Stress distribution')
ylabel('\tau (\rho)')
xlabel('\rho')

figure

subplot(1, 2, 1)
contour(X, Y, tau_ro(R), 100)
colorbar
axis equal tight
grid on
xlabel('z')
ylabel('y')
title('\tau (x, y, z)')

subplot(1, 2, 2)
contour(X, Y, gam_ro(R), 100)
colorbar
axis equal tight
grid on
xlabel('z')
ylabel('y')
title('\gamma (x, y, z)')