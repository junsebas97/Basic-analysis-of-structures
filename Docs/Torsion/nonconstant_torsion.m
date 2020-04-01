clc, clear all , close all
%{
This codes assess the stress and strain distribution in the cross section
of a circular or tubular bar subjected to non - constant torsion

MADE BY: Juan Sebastián Delgado Trujillo
%}

%% problem parameters:

MD    = 20;                 % [kN - m/m]
G     = 84e6;               % [kPa]
L     = 4;                  % [m]
r_ext = 3*0.0254;           % [m]
r_int = @(x) 0.0254*x/2;    % [m]

%% main:

Rm  = -MD*L;             % torsion moment reaction in the support
T_x = @(x) MD*x + Rm;    % internal torque

Ip = @(x) pi*(r_ext^4 - r_int(x).^4)/2;    % polar inertia

% the ratio of rotation and the rotation of the final section are
% calculated
theta = @(x) T_x(x)./(Ip(x).*G);
phi   = integral(theta, 0, L);

% equation of the stress and strain distribution in the section
tau   = @(x, ro) T_x(x).*ro./Ip(x);
gamma = @(x, ro) tau(x, ro)./G;

%% plots ans reports:

fprintf('\nrotation of the final section = %.3f°\n', phi*180/pi)

figure
hold on
eval = [0, L/2, 3*L/4];

for i = 1:3
    
    x_eval = eval(i);
    
    % this block is for creating the circular domain
    nR = linspace(r_int(x_eval), r_ext) ;        nT = linspace(0,     2*pi);
    [R, MD] = meshgrid(nR,nT) ;
    X = R.*cos(MD);                       Y = R.*sin(MD);

    subplot(3, 3, 3*i)
    hold on
    plot(nR, tau(x_eval, nR), 'b-', -nR, tau(x_eval, nR), 'b-')
    axis tight
    grid on
    title_text = strcat('Stress distribution in x = ', num2str(x_eval),' m');
    title(title_text)
    ylabel('\tau (\rho)')
    xlabel('\rho')

    subplot(3, 3, 3*i - 2)
    contour(X, Y, tau(x_eval, R), 100)
    colorbar
    axis equal tight
    grid on
    xlabel('z')
    ylabel('y')
    title_text = strcat('\tau (', num2str(x_eval),' m , y, z)');
    title(title_text)

    subplot(3, 3, 3*i - 1)
    contour(X, Y, gamma(x_eval, R), 100)
    colorbar
    axis equal tight
    grid on
    xlabel('z')
    ylabel('y')
    title_text = strcat('\gamma (', num2str(x_eval),' m , y, z)');
    title(title_text)
end