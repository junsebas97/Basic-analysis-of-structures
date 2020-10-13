clear all, close all, clc
X = 1; Y = 2;
%{
This codes solves a truss using the stiffness method. The loads, the
geometry of the structure and the properties of the elements must be passed
to the code

MADE BY: Juan Sebastián Delgado Trujillo
%}

%% Problem's parameters:

% nodal coordinates (x[m], y[m]):
coord = [0,   0;
         2.5, 0;
         2.5, 3];

% bars properties:
E = [100, 150, 200]*10^6;    % elastic modulus  [kPa]
%A = [0.1, 0.4, 0.2];         % area of the bars [m2]
A = [0.15*0.15, 0.05*0.5, 0.05*0.5];         % area of the bars [m2]

% correspondence matrix:
LaG = [1, 2;
       2, 3;
       1, 3];

% loads:
P = [0; 0; 0; 0; 1e6*cosd(60); 1e6*sind(60)];

% supports:
cc = [1, 2, 3, 4];

n_nodos = size(coord, 1);    % number of nodes
n_el    = size(E, 2);        % number of elements
n_gdl   = 2*n_nodos;         % number of degree of freedoms


gdl = [(1:2:n_gdl)', (2:2:n_gdl)'];    % degree of freedoms

dd  = setdiff((1:n_gdl), cc);          % free degree of freedom

%% plots:
figure
hold on
for el = 1:n_el
    plot(coord(LaG(el, :), X), coord(LaG(el, :), Y), 'k', 'LineWidth', 1.5)
end
title 'Cercha'
axis equal
grid on

%% Stiffness matrix:

K = zeros(n_gdl, n_gdl);

% in each element of the truss:
for el = 1:n_el
    
    % 1)the coordinates are extracted:
    x_i = coord(LaG(el, 1), X);        x_f = coord(LaG(el, 2), X);
    y_i = coord(LaG(el, 1), Y);        y_f = coord(LaG(el, 2), Y);
    
    % 2) the length and sin(theta) - cos(theta) are calculated:
    L = sqrt((x_f - x_i)^2 + (y_f - y_i)^2);
    c = (x_f - x_i)/L;
    s = (y_f - y_i)/L;
    
    % 3) The stiffness matrix in local coordinates and the transformation
    % matrix are assembled:
    
    k = (E(el)*A(el)/L)*[ 1, 0, -1, 0;
                          0  0,  0, 0;
                         -1, 0,  1, 0;
                          0, 0,  0, 0];
    T = [ c, s,  0, 0;
         -s, c,  0, 0;
          0, 0,  c, s;
          0, 0, -s, c];
      
    % 4) it's transformed the stiffness matrix of the element to global
    % coordinates
    K_e = T'*k*T;
    
    % 5) the degrees of fredom of the element are extracted
    idx = gdl(LaG(el, :), :)';
    
    % 5) the elemental matrix is taken to the structure matrix
    K(idx, idx) = K(idx, idx) + K_e;
end

%% System solution:

% when there is the structure's stiffness equation it's divided to
% according to the restraints:
Pc = P(dd);    K_cc = K(cc, cc);    K_cd = K(cc, dd);
               K_dc = K(dd, cc);    K_dd = K(dd, dd);     Uc = zeros(4, 1);
                
% then, displacement and the reactions are calculated:
Ud = inv(K_dd)*(Pc - K_dc*Uc);
Pd = K_cc*Uc + K_cd*Ud;

U     = zeros(n_gdl, 1);        U(cc) = Uc;        U(dd) = Ud;
P(cc) = Pd;

% the new shape of the structure is assessed
disp     = reshape(U, 2, n_nodos)';
deformed = coord + disp;

%% Internal forces:

% finally, the internal forces of each bar are calculated:
p = zeros(4, n_el);
for el = 1:n_el
    
    % nodal displacements of the bar:
    u_e = reshape(U(gdl(LaG(el, :), :))', 4, 1);
    
    % transformation and stiffness matrices:
    x_i = coord(LaG(el, 1), X);        x_f = coord(LaG(el, 2), X);
    y_i = coord(LaG(el, 1), Y);        y_f = coord(LaG(el, 2), Y);
    
    L = sqrt((x_f - x_i)^2 + (y_f - y_i)^2);
    c = (x_f - x_i)/L;                 s = (y_f - y_i)/L;
    
    k = (E(el)*A(el)/L)*[ 1, 0, -1, 0;
                          0  0,  0, 0;
                         -1, 0,  1, 0;
                          0, 0,  0, 0];
    T = [ c, s,  0, 0;
         -s, c,  0, 0;
          0, 0,  c, s;
          0, 0, -s, c];
    
    % force vector of the element:
    p(:, el) = k*T*u_e;
end

%% Results:

figure
hold on
for el = 1:n_el
    plot(   coord(LaG(el, :), X),    coord(LaG(el, :), Y), 'k--')
    plot(deformed(LaG(el, :), X), deformed(LaG(el, :), Y), 'r-',...
         'LineWidth', 2)
end
title 'Deformed shape'
axis equal
grid on

for el = 1:n_el
    fprintf("\nForce vector of the El.%d:\n", el)
    p(:, el)
end