%{
This codes solves a simetric cable subjected to a uniform load in the
horizontal dimension. The deformed shape and tensions are assessed 

MADE BY: Juan Sebastián Delgado Trujillo
%}
clc, clear all, close all

%% problem parameters:

%h = 50;
h = 25;       % cable heigth [m]
L = 100;      % cable length [m]
w = 65;       % uniform load [kN/m]

L_t = L/2;    % length of the rigth section [m]

max_t  = 5000;    % failure tension of the cable [kN]

%% reactions:

Rh = (w*L_t^2)/(2*h);    % horizontal reaction at the support [kN]
Rv = w*L_t;              % vertical reaction at the support   [kN]

%% cable tension:

T_x   = Rh;                              % tension x-component [kN]
T_y   = @(x) w*x;                        % tension y-component [kN]
T_mag = @(x) sqrt(T_x^2 + T_y(x).^2);    % tension magnitude   [kN]

T_max = T_mag(L_t);     % maximum tension   [kN]

%% shape:

y = @(x) w/(2*T_x)*x.^2;

%% sections: 

% the analysis is extended to the left section beacause it follows the same
% rules and relations.

x_lim = [-L_t, L_t];    % [left section, rigth section]
                        % [  (-7.5, 0),     (0, 7.5)  ]
                        
%% plots:

fprintf('Support reactions:')
fprintf('\nRx = %.3f kN        Ry = %.3f kN\n', Rh, Rv)
fprintf('Maximum tension:')
fprintf('\nT max = %.3f kN\n', T_max)

figure
fmesh(@(u, v) u, @(u, v) y(u), @(u, v) T_mag(u), [x_lim(1), x_lim(2), 0, h], 'LineWidth', 2);
view([0, 0, 1])
colorbar
caxis([0, max_t])
axis equal
xlabel('x [m]')
ylabel('y [m]')
title 'Deformed cable and cable stress magnitude'

