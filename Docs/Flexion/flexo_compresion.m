%{
This code solves a recatangular beam fixed in the base and free in the top
that is subjected simultaneously to bending and axial using the
superposition principle

Made by: Juan Sebastián Delgado Trujillo
%}

clc, clear all, close all
syms N(y) V(y) M(y) th(y) u(y) x y

%% Beam's properties:

E = 200e6;          % Young modulus [kPa]
b = 0.2;            % base          [m]
h = 0.35;           % height        [m]
L = 2;              % span          [m]

q_b = -5*y + 10;    % perpendicular force [kN/m]
q_a = -8;            % axial force         [kN/m]

I = (b*h^3)/12;     % moment of inertia [m4]
A = b*h;            % area              [m2]

%% bending:

% the bending "effect" is calculated using the equation of the beam and the
% respective conditions of the supports
sol = dsolve(diff(V, y) == -q_b,    ...
             diff(M, y) == V,       ...
             diff(th,y) == M/(E*I), ... 
             diff(u, y) == th,      ...
             u(0)  == 0,  th(0) == 0, M(L) == 0, V(L) == 0);

% the solution is stored in other variables to improve the readibility
V_b = sol.V;        M_b = sol.M;        th_b = sol.th;        u_b = sol.u;

s_bend   = -M_b*x/I;                       % axial stresses
tau_bend = (V_b/(2*I))*((h/2)^2 - x^2);    % shear stresses

%% axial:

% the axial "effect" is calculated with the bar equation and the fact that
% in the top there is no axial reaction:
N_a = dsolve(diff(N, y) == -q_a, N(L) == 0);

s_axial = N_a/A;   % axial stresses

%% simultaneous action:

% the variables that are present in both states are summed and the other
% variables are kept
tau_y = tau_bend;                          s_y = s_bend + s_axial;
V_y   = V_b;
M_y   = M_b;
th_y  = th_b;
u_y   = u_b;
N_y   = N_a;

% the stress distribution (normal with axial, normal without axial and
% shear) are estimated in the middle of the span
s_l2   = subs(s_y, y, L/2);
sb_l2  = subs(s_bend, y, L/2);
tau_l2 = subs(tau_y, y, L/2);
%% plots:

figure
title('Variables that dont change')
subplot(3, 4, 1)
hold on
fplot(V_y, [0, L], 'LineWidth', 2)
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
axis tight
grid on
xlabel('y [m]')
ylabel('Shear force [kN]')

subplot(3, 4, 2)
hold on
fplot(M_y, [0, L], 'LineWidth', 2)
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
axis tight
grid on
xlabel('y [m]')
ylabel('Bending moment [kN - m]')

subplot(3, 4, 3)
hold on
fplot(th_y, [0, L], 'LineWidth', 2)
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
camroll(90)
axis tight
grid on
xlabel('y [m]')
ylabel('\theta [rad]')

subplot(3, 4, 4)
hold on
fplot(u_y, [0, L], 'LineWidth', 2)
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
camroll(90)
axis tight
grid on
xlabel('y [m]')
ylabel('\sl u [m]')

subplot(3, 4, 6)
hold on
fplot(N_y, [0, L], 'LineWidth', 2)
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
camroll(90)
axis tight
grid on
xlabel('y [m]')
ylabel('\sl Axial reaction [kN]')


subplot(3, 4, 10)
fsurf(tau_y, [-h/2, h/2, 0, L])
view([0, 0, 1])
xlabel('y [m]')
ylabel('\tau (x, y) [kPa]')

subplot(3, 4, 11)
hold on
fplot(tau_l2, [-h/2, h/2], 'LineWidth', 2)
plot([-h/2, h/2], [0, 0], 'k', 'LineWidth', 2)
axis tight
grid on
ylabel('\tau_{y}(x, L/2) [kPa]')
xlabel('x [m]')


figure
subplot(1, 3, 1)
fsurf(s_y, [-h/2, h/2, 0, L])
view([0, 0, 1])
xlabel('x [m]')
ylabel('\sigma_{x}(x, y) [kPa]')
colorbar

subplot(1, 3, 2)
hold on
fplot(sb_l2, [-h/2, h/2], 'LineWidth', 2)
plot([-h/2, h/2], [0, 0], 'k', 'LineWidth', 2)
axis tight
grid on
ylabel('\sigma_{y}(x, L/2) [kPa]')
xlabel('x [m]')
title 'Normal stress distribution with bending'

subplot(1, 3, 3)
hold on
fplot(s_l2, [-h/2, h/2], 'LineWidth', 2)
plot([-h/2, h/2], [0, 0], 'k', 'LineWidth', 2)
axis tight
grid on
ylabel('\sigma_{y}(x, L/2) [kPa]')
xlabel('x [m]')
title 'Normal stress distribution with bending and axial'

fprintf('\nNormal stress distribution with bending\n')
pretty(sb_l2)
fprintf('\nNormal stress distribution with bending and axial\n')
pretty(s_l2)