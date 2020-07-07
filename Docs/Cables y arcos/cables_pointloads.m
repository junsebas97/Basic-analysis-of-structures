clc, clear all , close all
%{
This codes solves a cable subjected to 2 point loads when the angle at the
first support is known. Tensile forces and stresses in the cable are
calculate and its deformed shape ploted

MADE BY: Juan Sebastián Delgado Trujillo
%}

%% problem parameters:

A       = 2*0.0254;    % area of the cable       [m2]
theta_a = 30;          % angle of the cable at A [grad]
%theta_a = 60;
%theta_a = 15;

Fc      = 10;    % [kN]
Fd      = 15;    % [kN]

% known coordinates of the points [m]
xa = 0;        ya = 0;
xb = 8;        yb = 1;
xc = 2;
xd = 6;

%% equilibrium:

% in the support A the reactions are expressed in terms of the AC segment
% tension
syms Tac
Ra_y = Tac*sind(theta_a);
Ra_x = Tac*cosd(theta_a);

% the equilibrium in the whole structure is formulated to assess the
% tension in the AC segment
sum_mom = Fc*(xb - xc) + Fd*(xb - xd) - Ra_x*yb - Ra_y*(xb - xa) == 0;
Tac     = double(solve(sum_mom));

% the magnitude of the reactions at support A are numerically assessed
Ra_y = double(subs(Ra_y, Tac));
Ra_x = double(subs(Ra_x, Tac));

% the equillibrium at the point C let us assess the tension components of
% the CD segment, its tension magnitude and the segment orientation
Tcd_x = Tac*cosd(theta_a);
Tcd_y = Fc - Tac*sind(theta_a);

Tcd     = sqrt(Tcd_x^2 + Tcd_y^2);
theta_c = atan2d(Tcd_y, Tcd_x);

% equillibrium at the point D to CD segment
Tdb_x = Tcd_x;
Tdb_y = Fd + Tcd_y;

Tdb     = sqrt(Tdb_x^2 + Tdb_y^2);
theta_d = atan2d(Tdb_y, Tdb_x);

% reactions in the support B are obtained with equilibrium at that point
Rb_x = Tdb_x;        Rb_y = Tdb_y;

%% stresses:

% stress formula for constant axial force [kPa]
s_ac = Tac/A;       s_cd = Tcd/A;         s_db = Tdb/A;

%% deflections:

% using the orientation of the cables the vertical displacement at points
% where the forces act is computed
yc = -xc*tand(theta_a);          yd = yb - (xb - xd)*tand(theta_d);


%% results:

fprintf('\nReactions\n')
fprintf('Rax = %.3f kN    Ray = %.3f kN\n', Ra_x, Ra_y)
fprintf('Rbx = %.3f kN    Rby = %.3f kN\n', Rb_x, Rb_y)

fprintf('\nTension in the cables\n')
fprintf('Tac = %.3f kN    sigma ac = %.3f kPa\n', Tac, s_ac)
fprintf('Tcd = %.3f kN    sigma cd = %.3f kPa\n', Tcd, s_cd)
fprintf('Tdb = %.3f kN    sigma db = %.3f kPa\n', Tdb, s_db)

figure
hold on
plot([xa, xc], [ya, yc], 'r-', [xc, xd], [yc, yd], 'r-', [xd, xb], [yd, yb], 'r-')
grid on
axis equal tight
xlabel 'x [m]'
ylabel 'y [m]'
title 'Deformed cable'