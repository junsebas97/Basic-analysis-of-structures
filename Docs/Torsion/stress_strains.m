clc, clear all , close all
%{
This codes plots the possible planes of failure and calcules the produced
distribution of normal strains in a tubular bar subjected to constant 
torsion

MADE BY: Juan Sebastián Delgado Trujillo
%}

%% problem parameters:

T     = 50;          % [kN - m]
E     = 210e6;       % [kPa]
miu   = 0.25;        % [-]
L     = 2.5;         % [m]
r_int = 1*0.0254;    % [m]
r_ext = 2.5*0.0254;  % [m]

%% main:

G = E/(2*(1 + miu));    % calculated shear modulus [kPa]

Ip  = pi*(r_ext^4 - r_int^4)/2;   % polar inertia [-]
phi = (T*L)/(Ip*G);               % rotation of the final section  [rad]

% equation of the shear stress and normal strain distribution in the cross
% section
tau_ro = @(ro) T*ro./Ip;
eps_ro = @(ro) tau_ro(ro)*(1 + miu)/E;

% the maximum stress (shear and normal) occur in the surface and have the
% same magnitud
tau_max = tau_ro(r_ext);
sig_max = tau_max;

%% plots ans reports:

fprintf('\nrotation of the final section = %.3f°', phi*180/pi)
fprintf('\ntau_max   = %.3f MPa', tau_max/1000)
fprintf('\nsigma_max = %.3f MPa\n', sig_max/1000)

% this block is for creating the circular domain
nR = linspace(r_int, r_ext) ;        nT = linspace(0,     2*pi);
[R, T] = meshgrid(nR,nT) ;
X = R.*cos(T);                       Y = R.*sin(T);

figure
contour(X, Y, eps_ro(R), 100)
colorbar
axis equal tight
grid on
xlabel('z')
ylabel('y')
title('\epsilon_{\theta} (x, y, z)')

% this block creates the cylinder
n_levels     = 25;
[Xc, Yc, Zc] = cylinder(r_ext, n_levels);      
temp = Zc(2, :);
for i = 2:n_levels
    Xc(i, :) = Xc(1, :);
    Yc(i, :) = Yc(1, :);
    Zc(i, :) = temp*(L*(i - 1)/(n_levels - 1));
end
temp = Xc;         Xc = Zc;        Zc = temp;
temp = ones(size(Zc));


figure

subplot(3, 1, 1)
hold on
surf(Xc, Yc, Zc)
colormap([0.85, 0.85, 0.85])
quiver3(Xc, Yc, Zc, temp*tau_max, -Zc*0, Yc*0, 0.5, 'LineWidth', 1.2)
grid on 
axis tight
view([0.05, -1, 0.75])
title 'First failure direction of \tau_{max}'
xlabel 'x(m)'
ylabel 'y(m)'
zlabel 'z(m)'

subplot(3, 1, 2)
hold on
surf(Xc, Yc, Zc)
colormap([0.85, 0.85, 0.85])
quiver3(Xc, Yc, Zc, temp*0, -Zc*tau_max, Yc*tau_max, 0.25, 'LineWidth', 1.2)
grid on 
axis tight
view([0.25, -1, 0.25])
title 'Second failure direction of \tau_{max}'
xlabel 'x(m)'
ylabel 'y(m)'
zlabel 'z(m)'

subplot(3, 1, 3)
hold on
surf(Xc, Yc, Zc)
colormap([0.85, 0.85, 0.85])
quiver3(Xc, Yc, Zc, temp*tau_max, -Zc*tau_max, Yc*tau_max, 1, 'LineWidth', 1.2)
grid on 
axis tight
view([0.15, -1, 0.05])
title 'ailure direction of \sigma_{\theta max}'
xlabel 'x(m)'
ylabel 'y(m)'
zlabel 'z(m)'