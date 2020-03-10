%{
This code assess the diagrams of: shear force, bending moment, rotation
angle, deflection, maximum normal stress and maximum shear stress, of a
simply supported beam with rectangular constant cross section subjected to
a constant distributed load

Made by: Juan Sebastián Delgado Trujillo
%}

clc, clear all, close all

%% Beam's properties:

E = 200e6;          % Young modulus [kPa]
b = 0.2;            % base   [m]
h = 0.3;            % height [m]
I = (b*h^3)/12;     % moment of inertia [m4]
w = 20;             % distributed load  [kN/m]
l = 2;              % span [m]

%% Diagrams:
% analytic solutions -equations- of the differential equations system

V   = @(x) -w*x + w*l/2;
M   = @(x) (-w*x^2)/2 + w*l*x/2;
th  = @(x) ((-w*x^3)/6 + (w*l*x^2)/4 - (w*l^3)/24)/(E*I); 
def = @(x) ((-w*x^4)/24 + (w*l*x^3)/12 - (w*x*l^3)/24)/(E*I);

%% Stresses:

% normal stresses
sx_y   = @(x, y) -M(x)*y./I;
sx_min = @(x)    sx_y(x, b/2);        sx_max = @(x) sx_y(x, -b/2);
sx_l2  = @(y)    sx_y(l/2, y);

% shear stresses:
tx_y   = @(x, y) (V(x)/(2*I))*((h^2/4) - y.^2);
tx_max = @(x)    tx_y(x, 0);
tx_ap  = @(y)    tx_y(0, y);

%% Plots: 

figure
subplot(4, 1, 1)
fplot(V, [0, l])
axis tight
grid on
xlabel('x [m]')
ylabel('Shear force [kN]')

subplot(4, 1, 2)
fplot(M, [0, l])
axis tight
grid on
xlabel('x [m]')
ylabel('Bending moment [kN - m]')

subplot(4, 1, 3)
fplot(th, [0, l])
axis tight
grid on
xlabel('x [m]')
ylabel('\theta [rad]')

subplot(4, 1, 4)
hold on
fplot(def, [0, l])
plot([0, l], [0, 0], 'r--')
axis tight
grid on
xlabel('x [m]')
ylabel('\sl v [m]')



figure
subplot(3, 2, 1)
ezsurf(sx_y, [0, l, -h/2, h/2])
view([0, 0, 1])
xlabel('x [m]')
ylabel('\sigma_{x}(x, y) [kPa]')

subplot(3, 2, 3)
hold on
fplot(sx_min, [0, l])
fplot(sx_max, [0, l], 'go-')
axis tight
grid on
legend('\sigma_{x}min', '\sigma_{x}min', 'Location', 'best')
xlabel('x [m]')
ylabel('\sigma_{x}[kPa]')

subplot(3, 2, 5)
hold on
plot(sx_l2(linspace(-h/2, h/2)), linspace(-h/2, h/2))
axis tight
grid on
ylabel('y [m]')
xlabel('\sigma_{x}[kPa]')

subplot(3, 2, 2)
ezsurf(tx_y, [0, l, -h/2, h/2])
view([0, 0, 1])
xlabel('x [m]')
ylabel('\tau (x, y) [kPa]')

subplot(3, 2, 4)
hold on
fplot(tx_max, [0, l])
axis tight
grid on
xlabel('x [m]')
ylabel('\tau_{max} [kPa]')

subplot(3, 2, 6)
hold on
plot(tx_ap(linspace(-h/2, h/2)), linspace(-h/2, h/2))
axis tight
grid on
ylabel('y [m]')
xlabel('\tau_{x}[kPa]')
