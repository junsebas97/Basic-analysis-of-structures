clc, clear all, close all
%{
This codes solves a rectangular beam subjected to biaxial bending due to an
inclined load. The sigma_x equation and the shear, moment, theta and
deflection diagrams in each direction are reported

MADE BY: Juan Sebastián Delgado Trujillo
%}

syms x My(x) Vy(x) thy(x) uy(x) Mz(x) Vz(x) thz(x) uz(x)
%% beam properties:

b = 0.3;    % base            [m]
h = 0.5;    % height          [m]
E = 210e6;  % elastic modulus [kPa]
lz = 6;     % span            [m]

Iy = b*h^3/12;   % inertia respect y axis [m4]
Iz = h*b^3/12;   % inertia respect z axis [m4]

%% loads:

q_x  = -2*x;    % magnitude of the distributed load [kN/m]
beta = 60;      % inclination of the load respect to z axis [degress]

% distributed loads in each dimention
q_y = q_x*sind(beta); 
q_z = q_x*cosd(beta);

%% Euler - Bernoulli beam equations:

% the EB equation is solved for a double-fixed supported beam, i.e. no 
% rotation nor displacements in the support

% EB for y-bending
solu_y = dsolve(diff(Vz(x),  x) == q_y,          ...
                diff(Mz(x),  x) == Vz(x),        ...
                diff(thz(x), x) == Mz(x)/(E*Iy), ...
                diff(uz(x),  x) == thz(x),       ...
                thz(0) == 0, thz(lz) == 0, uz(0) == 0, uz(lz) == 0);

Mz_x  = matlabFunction(solu_y.Mz);

% EB for z-bending
solu_z = dsolve(diff(Vy(x),  x) == q_z,          ...
                diff(My(x),  x) == Vy(x),        ...
                diff(thy(x), x) == My(x)/(E*Iz), ...
                diff(uy(x),  x) == thy(x),       ...
                thy(0) == 0, thy(lz) == 0, uy(0) == 0, uy(lz) == 0);
          
My_x  = matlabFunction(solu_z.My);

%% peprpendicular reactions:

%these reactions are not colinear so their resultant must be calculated
% with the algebra of vectors
V_x  = sqrt(solu_y.Vz^2  + solu_z.Vy^2);
M_x  = sqrt(solu_y.Mz^2  + solu_z.My^2);
th_x = sqrt(solu_y.thz^2 + solu_z.thy^2);
u_x  = sqrt(solu_y.uz^2  + solu_z.uy^2);

%% stress field:

% the normal stress field is obtained with the superposition method and
% it's evaluated in 4 representative ponits
sx = @(x, y, z) -(Mz_x(x)*y./Iy + My_x(x)*z./Iz);

sx_0  = @(y, z) sx(0,    y, z);
sx_l4 = @(y, z) sx(lz/4, y, z);
sx_l3 = @(y, z) sx(lz/3, y, z);
sx_l2 = @(y, z) sx(lz/2, y, z);

% the neural axis:
na    = @(x, z) -z.*(Mz_x(x)*Iz)/(My_x(x)*Iy);
na_0  = @(y) na(0,    y);
na_l4 = @(y) na(lz/4, y);
na_l3 = @(y) na(lz/3, y);
na_l2 = @(y) na(lz/2, y);

%% plots:

fprintf('\nShear force magnitude in biaxial bending:\n')
pretty(simplify(V_x))
fprintf('\nBending moment magnitude in biaxial bending:\n')
pretty(simplify(M_x))
fprintf('\nAngle magnitude in biaxial bending:\n')
pretty(simplify(th_x))
fprintf('\ntotal deflection in biaxial bending:\n')
pretty(simplify(u_x))

figure

subplot(2, 2, 1)
hold on
plot3([0, lz] , [0, 0], [0, 0])
fplot3(x, solu_y.Vz, solu_z.Vy, [0, lz], 'Linewidth', 2)
fplot3(x, @(x) 0,    solu_z.Vy, [0, lz], 'k--')
fplot3(x, solu_y.Vz, @(x) 0,    [0, lz], 'k--')
view([1, -1, 0.75])
grid on
axis tight
xlabel('x [m]')
xlim([0, lz])
ylabel('V(x) [kN]')
title 'Shear force'

subplot(2, 2, 2)
hold on
plot3([0, lz] , [0, 0], [0, 0])
fplot3(x, solu_y.Mz, solu_z.My, [0, lz], 'Linewidth', 2)
fplot3(x, @(x) 0,    solu_z.My, [0, lz], 'k--')
fplot3(x, solu_y.Mz, @(x) 0,    [0, lz], 'k--')
view([1, -1, 0.5])
grid on
axis tight
xlabel('x [m]')
xlim([0, lz])
ylabel('M(x) [kN - m]')
title 'Bending moment'

subplot(2, 2, 3)
hold on
plot3([0, lz] , [0, 0], [0, 0])
fplot3(x, solu_y.thz, solu_z.thy, [0, lz], 'Linewidth', 2)
fplot3(x, @(x) 0,     solu_z.thy, [0, lz], 'k--')
fplot3(x, solu_y.thz, @(x) 0,     [0, lz], 'k--')
view([0.8, -1, 0.2])
grid on
axis tight
xlabel('x [m]')
xlim([0, lz])
ylabel('\theta(x) [rad]')
title 'Rotation'

subplot(2, 2, 4)
hold on
plot3([0, lz] , [0, 0], [0, 0])
fplot3(x, solu_y.uz, solu_z.uy, [0, lz], 'Linewidth', 2)
fplot3(x, @(x) 0,    solu_z.uy, [0, lz], 'k--')
fplot3(x, solu_y.uz, @(x) 0,    [0, lz], 'k--')
view([1, -1, 0.25])
grid on
axis tight
xlabel('x [m]')
xlim([0, lz])
ylabel('u(x) [m]')
title 'Deflection'


figure
subplot(2, 2, 1)
hold on
fsurf(sx_0, [-b/2, b/2, -h/2, h/2])
fplot(na_0, [-b/2, b/2], 'k-', 'Linewidth', 2)
view([0, 0, 1])
axis tight
xlabel('z [m]')
ylabel('y [m]')
legend('\sigma_x', 'Neutral axis')
title '\sigma_{x} (0, y, z)'

subplot(2, 2, 2)
hold on
fsurf(sx_l4, [-b/2, b/2, -h/2, h/2])
fplot(na_l4, [-b/2, b/2], 'k-', 'Linewidth', 2)
view([0, 0, 1])
axis tight
xlabel('z [m]')
ylabel('y [m]')
legend('\sigma_x', 'Neutral axis')
title '\sigma_{x} (L/4, y, z)'

subplot(2, 2, 3)
hold on
fsurf(sx_l3, [-b/2, b/2, -h/2, h/2])
fplot(na_l3, [-b/2, b/2], 'k-', 'Linewidth', 2)
view([0, 0, 1])
axis tight
xlabel('z [m]')
ylabel('y [m]')
legend('\sigma_x', 'Neutral axis')
title '\sigma_{x} (L/3, y, z)'

subplot(2, 2, 4)
hold on
fsurf(sx_l2, [-b/2, b/2, -h/2, h/2])
fplot(na_l2, [-b/2, b/2], 'k-', 'Linewidth', 2)
view([0, 0, 1])
axis tight
xlabel('z [m]')
ylabel('y [m]')
legend('\sigma_x', 'Neutral axis')
title '\sigma_{x} (L/2, y, z)'