clc, clear all, close all

%{
This codes solves a composite rectangular beam of two materials. The
deflections and the normal stresses are calculated and plotted

MADE BY: Juan Sebastián Delgado Trujillo
%}

%% beam properties:

b    = 0.2;    % base                         [m]
h_m1 = 0.3;    % height of the upper material [m]
h_m2 = 0.1;    % height of the lower material [m]

lz = 6;       % span [m]

w  = -10;     % distributed load [kN/m]

E1 = 21e6;    % elastic modulus of the upper material [kPa]
E2 = 210e6;   % elastic modulus of the lower material [kPa]

h_tot = h_m1 + h_m2;    % beam's total height [m]
n     = E2/E1;          % modular ratio

%% neutral axis:
syms h1 y_var z

% the integral of the "y" in each material is formulated and the neutral
% axis rule is solved
int_1 = int(int(y_var, h1 - h_m1,  h1),  z, b, 0);
int_2 = int(int(y_var, h1 - h_tot, h1 - h_m1), z, b, 0);

d_tna = solve(int_1 + n*int_2 == 0);
d_tna = double(d_tna);

%% inertia:

% moment of inertia respect to the composite section neutral axis of each
% material
I_m1 = int(int(y_var^2, d_tna - h_m1,  d_tna),        z, 0, b);
I_m2 = int(int(y_var^2, d_tna - h_tot, d_tna - h_m1), z, 0, b);

It = I_m1 + n*I_m2;    % total inertia of the composite section
It = double(It);

%% Euler - Bernoulli beam equation:
syms M(x) V(x) th(x) u(x)

% EB equation is solved for a simply supported beam, i.e. no moments nor
% displacements in the support
solu = dsolve(diff(V(x),  x) == w,            ...
              diff(M(x),  x) == V(x),         ...
              diff(th(x), x) == M(x)/(E1*It), ...
              diff(u(x),  x) == th(x),        ...
              M(0) == 0, M(lz) == 0, u(0) == 0, u(lz) == 0);
          
M_x = matlabFunction(solu.M);
u_x = solu.u;

%% stress field:

% the stress in the transformed section are obtained, these stress match
% with the real stresses unchanged material but must be magnified in the
% transformed materials
sx_1 = @(x, y) -M_x(x)*y./It;
sx_2 = @(x, y) sx_1(x, y)*n;

% the stress are assesed in a cut in the mid of the span
y1 = linspace(d_tna - h_m1,  d_tna);            sx1_l2 = sx_1(lz/2, y1);
y2 = linspace(d_tna - h_tot, d_tna - h_m1);     sx2_l2 = sx_2(lz/2, y2);

%% plots:

figure()
fplot(u_x)
grid on
axis tight
xlim([0, lz])
xlabel('x [m]')
ylabel('v(x) [m]')
title 'Deflection'

figure
subplot(2, 1, 1)
hold on
fsurf(sx_1, [0, lz, d_tna - h_m1,  d_tna])
fsurf(sx_2, [0, lz, d_tna - h_tot, d_tna - h_m1])
view([0, 0, 1])
axis tight
xlabel('x [m]')
ylabel('\sigma_{x}(x, y) [kPa]')
ylim([d_tna - h_tot, d_tna])
title 'Normal stress in the composite beam'

subplot(2, 1, 2)
hold on
stem(y1, sx1_l2)
stem(y2, sx2_l2);
grid on
axis tight
ylabel('Y [m]')
xlabel('\sigma_{x}(l/2, y) [kPa]')
title 'Normal stress in x = L/2'