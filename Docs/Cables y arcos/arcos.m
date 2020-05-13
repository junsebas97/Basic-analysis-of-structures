%{
This codes calcules the internal reactions of a three hinged arch with 
semicircular shape and subjected to a uniform distributed load  

MADE BY: Juan Sebastián Delgado Trujillo
%}
clc, close all, clear all

%% data:

L = 3;      % arch span (diameter) [m]
R = 1.5;    % arch heigth (radius) [m]
w = 50;     % distributed force    [kN/m]

l_sa = L/2;    % semi-arch span    [m]

%% equilibium:

P_eq = w*l_sa;    % equivalent point load in each section [kN]

syms Cx_s Cy_s

% with the moment equillibrium in the left support the Cx reaction is
% ESTIMATED
S_Ma = -P_eq*l_sa/2 + Cy_s*l_sa + Cx_s*R == 0;
Cx   = solve(S_Ma, Cx_s);

% then, the moment equilibrium is used to CALCULE the magnitude of the Cy
% reaction
S_Mb = P_eq*l_sa/2 + Cy_s*l_sa - Cx*R == 0;
Cy   = double(solve(S_Mb, Cy_s));

% when the magnitude of Cy is known the Cx magnitude is calculated
Cx = double(subs(Cx, Cy_s, Cy));

% finally, the forces equilibrium is used to assess the supports reactions
Ay = P_eq - Cy;                 By = P_eq + Cy;
Ax = Cx;                        Bx = Cx;

%% Arch shape:
syms x

% the shape of the arch is used to obtain the orientation of the cross-area
%y     = sqrt(R^2 - x^2);
y    = -2*x^2/3 + 1.5;
dydx  = diff(y);
theta = atan(dydx);

%% internal reactions:

% using the equilibrium in a segment of the right arch the internal forces
% are computed 

M_x = Cx*(1.5 - y) - w*x^2/2;    % bending moment [kN- m]

syms N

V1 = ( Cx  + N*cos(theta))/sin(theta);
V2 = (-w*x - N*sin(theta))/cos(theta);

N_x = solve(V1 == V2, N);
V_x = subs(V1, N, N_x);

