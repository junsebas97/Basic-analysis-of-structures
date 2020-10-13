clear all, close all, clc
%{
Code to calcule the stiffness matrix OR the nodal equivalent forces of a 2D
beam

MADE BY: Juan Sebastián Delgado Trujillo
NOTE: This code is based on
      https://github.com/diegoandresalvarez/elementosfinitos/blob/master/codigo/repaso_matricial/portico_2d/deduccion_K_y_fe_elemento_portico_2D/c1_deduccion_K_portico2D.m
%}
syms EI x L V(x) M(x) th(x) u(x) w
%% Stiffness matrix:

% boundary conditions to be analyzed:
case_1 = [1, 0, 0, 0];                       case_2 = [0, 1, 0, 0];
case_3 = [0, 0, 1, 0];                       case_4 = [0, 0, 0, 1];
cases  = [case_1; case_2; case_3; case_4];

% if the stiffness matrix will be assessed the no load is applied to the
% beam
q_x = 0;

%% Equivalent forces:
%{
q_x = w*x;

% the nodal forces are caluclated only in the undeformed beam:
case_1 = [0, 0, 0, 0];
cases  = [case_1];
%}

%% Main: 

n_cas = size(cases, 1);        % number of cases
Mat   = sym(zeros(4, n_cas));

% for each case:
for i = 1:n_cas
    
    % 1) the initial conditions are extracted
    cs = cases(i, :);
    
    % 2) with thos conditions the beam is solved:
    solu = dsolve(diff(V,  x) == q_x,  ... %-------------------------------
                  diff(M,  x) == V,    ... %             Beam
                  diff(th, x) == M/EI, ... %           equation
                  diff(u,  x) == th,   ... %-------------------------------
                  u(0) == cs(1), th(0) == cs(2), ...    Boundary conditions
                  u(L) == cs(3), th(L) == cs(4));%-------------------------
    
    % 3) the shear force and bending moment of the solution are considered
    V_x = solu.V;               M_x = solu.M;
    
    % 4) the equations are evaluated at the beam ends to obtain the
    % coefficients
    Mat(:, i) = [ subs(V_x, x, 0);     % Ry on the left support
                 -subs(M_x, x, 0);     % M  on the left support
                 -subs(V_x, x, L);     % Ry on the right support
                  subs(M_x, x, L)];    % M  on the right support
end
pretty(Mat)