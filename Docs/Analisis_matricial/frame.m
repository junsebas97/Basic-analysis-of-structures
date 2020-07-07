clear all, close all, clc
X = 1; Y = 2; Th = 3;
%{
his codes solves a frame using the stiffness method. The loads, the
geometry of the structure and the properties of the elements must be passed
to the code

MADE BY: Juan Sebastián Delgado Trujillo
%}

%% Problem's parameters:

% nodal coordinates (x[m], y[m], theta[rad]):
coord = [0,   0, 0;
         3,   0, 0;
         1.5, 2, 0;
         3,   2, 0];

% bars properties:
E_mod  = [210,  210, 210]*10^6;    % elastic modulus  [kPa]
base   = [0.1,  0.1, 0.1];         % section's base   [m2]
height = [0.25, 0.3, 0.4];         % section's height [m2]

% correspondence matrix:
LaG = [3, 4;
       2, 4;
       1, 3];

% loads:
% nodal loads[kN]
P = [0;      % Px_1
     0;      % Py_1
     0;      % M_1
     0;      % Px_2
     0;      % Py_2
     0;      % M_2
     0;      % Px_3
     0;      % Py_3
     0;      % M_3
     250;    % Px_4
     0;      % Py_4
     0];     % M_4
w = [0; 0; 100];                   % distributed loads [kN/m]
     
% supports:
cc = [1, 2, 3, 4, 5, 6];

n_nodos = size(coord, 1);    % number of nodes
n_el    = size(E_mod, 2);    % number of elements
n_gdl   = 3*n_nodos;         % number of degree of freedoms


gdl = [(1:3:n_gdl)', (2:3:n_gdl)', (3:3:n_gdl)'];    % degree of freedoms

dd  = setdiff((1:n_gdl), cc);          % free degree of freedom
%% plots:
figure
hold on
for el = 1:n_el
    plot(coord(LaG(el, :), X), coord(LaG(el, :), Y), 'k', 'LineWidth', 1.5)
end
title 'Portico'
axis equal
grid on

%% Stiffness matrix:

K    = zeros(n_gdl, n_gdl);
k_st = cell(el);
t_st = cell(el);
len  = cell(el);

% in each element of the truss:
for el = 1:n_el
    
    % 1)the coordinates, the section's lengths and the elasticity modulus
    % are extracted:
    x_i = coord(LaG(el, 1), X);        x_f = coord(LaG(el, 2), X);
    y_i = coord(LaG(el, 1), Y);        y_f = coord(LaG(el, 2), Y);
    b   = base(el);                    h   = height(el);
    E   = E_mod(el);
    
    % 2) the length and sin(a), cos(a), area and inertia are calculated:
    L = sqrt((x_f - x_i)^2 + (y_f - y_i)^2);
    c = (x_f - x_i)/L;
    s = (y_f - y_i)/L;
    A = b*h;
    I = b*h^3/12;
    
    % 3) The stiffness matrix in local coordinates and the transformation
    % matrix are assembled:
    
    k = [ E*A/L,            0,           0,-E*A/L,            0,           0;
              0, 12*E*I/(L^3), 6*E*I/(L^2),     0,-12*E*I/(L^3), 6*E*I/(L^2);
              0,  6*E*I/(L^2),     4*E*I/L,     0, -6*E*I/(L^2), 2*E*I/(L^2);
         -E*A/L,            0,           0, E*A/L,            0,           0;
              0,-12*E*I/(L^3),-6*E*I/(L^2),     0, 12*E*I/(L^3),-6*E*I/(L^2);
              0,  6*E*I/(L^2),     2*E*I/L,     0, -6*E*I/(L^2),     4*E*I/L];
          
    T = [ c, s, 0,  0, 0, 0;
         -s, c, 0,  0, 0, 0;
          0, 0, 1,  0, 0, 0;
          0, 0, 0,  c, s, 0;
          0, 0, 0, -s, c, 0;
          0, 0, 0,  0, 0, 1];
      
    % 4) it's transformed the stiffness matrix of the element to global
    % coordinates
    K_e = T'*k*T;
    
    % 5) the degrees of fredom of the element are extracted
    idx = gdl(LaG(el, :), :)';
    
    % 6) the elemental matrix is taken to the structure matrix
    K(idx, idx) = K(idx, idx) + K_e;
    
    % 7) the transformation and stiffness matrices and the length of the
    % element are stored to further steps
    k_st{el} = k;
    t_st{el} = T;
    len{el}  = L;
end

%% Nodal equivalent forces:

F    = zeros(n_gdl, 1);
f_st = cell(n_el);

% in each element:
for el = 1:n_el
    
    % 1) the distributed load, the transformation matrix and the length are
    % taken
    L = len{el};            dist = w(el);            T = t_st{el};
    
    % 2) the equivalent forces are calculated and transformed to global
    % coordinates
    f = [0; -dist*L/2; -dist*L^2/12; 0; -dist*L/2; dist*L^2/12];
    f_g = T'*f;
    
    % 3) these loads are assigned to global vector of equivalent forces
    idx    = reshape(gdl(LaG(el, :), :)', 6, 1);
    F(idx) = F(idx) + f_g;
    
    % 4) for future steps the nodal equivalent forces are stored
    f_st{el} = f;
end

% then, the nodal effective loads vector is computed
P_ef = P - F;

%% System solution:

% arrays are divided according to the degrees of freedom with unknown
% motion
Pc  = P_ef(dd);    
Kcc = K(cc, cc);    Kcd = K(cc, dd);    Kdc = K(dd, cc);    Kdd = K(dd, dd);
Uc  = zeros(length(cc), 1);
                
% then, displacement and the reactions are calculated:
Ud = inv(Kdd)*(Pc - Kdc*Uc);
Pd = Kcc*Uc + Kcd*Ud;

U     = zeros(n_gdl, 1);        U(cc) = Uc;        U(dd) = Ud;
P(cc) = Pd;

% the new shape of the structure is assessed
disp     = reshape(U, 3, n_nodos)';
scale    = 1000;
deformed = coord + disp*scale;

%% Internal forces:

% finally, the internal forces of each bar are calculated:
p = cell(n_el);
for el = 1:n_el    
    % nodal displacements of the bar:
    u_e = reshape(U(gdl(LaG(el, :), :))', 6, 1);
    
    % transformation matrix, stiffness matrix, and nodal equivalent forces:
    k = k_st{el};    T = t_st{el};    f = f_st{el};
    
    % force vector of the element:
    p{el} = k*T*u_e + f;
end

%% Results:

for el = 1:n_el
    b    = base(el);
    h    = height(el);
    E    = E_mod(el);
    dist = w(el);
    x1   = coord(LaG(el, 1), X);    y1 = coord(LaG(el, 1), Y);
    x2   = coord(LaG(el, 2), X);    y2 = coord(LaG(el, 2), Y);
    mov  = t_st{el}*reshape(U(gdl(LaG(el, :), :))', 6, 1);    
    
    c1_dibujar_barra_deformada_portico(b*h, E, b*h^3/12, x1, y1, x2, y2,...
    @(x) 0, @(x) dist, p{el}, mov, 1000, 0.01, 0.01, 0.01)
     
    fprintf("\nForce vector of the El.%d:\n", el)
    fprintf(['F_x1 = %1.3f\nF_y1 = %1.3f\nM_1  = %1.3f\n'...
             'F_x2 = %1.3f\nF_y2 = %1.3f\nM_2  = %1.3f\n'],...
             p{el}(1), p{el}(2), p{el}(3), p{el}(4), p{el}(5), p{el}(6))
end