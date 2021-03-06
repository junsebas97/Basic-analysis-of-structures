function c1_dibujar_barra_deformada_portico(A, E, I, x1, y1, x2, y2, ...
    qxloc,qyloc, qe, ae, esc_def, esc_faxial, esc_V, esc_M)
% Esta funcion dibuja el elemento de portico deformado junto con sus 
% respectivos diagramas de fuerza axial, fuerza cortante y momento flector.
%
% El diagrama de momento flector se grafica en el lado opuesto de la fibra
% a tracción
%
% PARAMETROS DE ENTRADA (junto con algunos ejemplos):
% A = area
% E = E
% I = Ix local
% (x1,y1) y (x2,y2) son las coordenadas de los puntos iniciales
% qxloc = @(x) x.^2; % carga en la dir. del eje x local (function handle)
% qyloc = @(x) 0;    % carga en la dir. del eje y local (function handle)
% qe = [ 0.01        % U1, V1, M1 reacciones del nodo 1 en coord. locales
%       -0.01
%        0.04
%       -0.01        % U2, V2, M2 reacciones del nodo 2 en coord. locales
%        0.02
%       -0.07 ];
% ae = [ 0.01        % u1, v1, t1 desplazamientos nodo 1 en coord. locales
%       -0.01
%        0.04
%       -0.01        % u2, v2, t2 desplazamientos nodo 2 en coord. locales
%        0.02
%       -0.07 ];
% esc_def    = 10;  % escalamiento de la deformada
% esc_faxial = 10;  % escalamiento del diagrama de axiales
% esc_V      = 10;  % escalamiento del diagrama de cortantes
% esc_M      = 10;  % escalamiento del diagrama de momentos

%% se definen algunas constantes
X = 1; Y = 2; X1  = 1; Y1 = 2; M1 = 3; X2 = 4; Y2 = 5; M2 = 6; 

%% resolver la ecuacion diferencial
npuntos = 1001;
xinit = linspace(0, hypot(x2-x1, y2-y1), npuntos);
sol   = bvpinit(xinit, zeros(6,1));
sol   = bvp5c(@ecuacion_diferencial, @condiciones_de_apoyo, sol);

%% Calculos intermedios
s     = sol.x;
axial = sol.y(6,:);          % Fuerza axial [kN]
V     = sol.y(4,:);          % Fuerza cortante [kN]
M     = sol.y(3,:);          % Momento flector [kN/m]
u     = sol.y(5,:);          % Desplazamiento horizontal de la viga [m]
v     = sol.y(1,:);          % Desplazamiento vertical de la viga [m]
%theta = atan(sol.y(2,:));   % Angulo de giro  [rad]

% rotacion de la solucion antes de dibujar
ang = atan2(y2-y1, x2-x1);
T   = [ cos(ang)  -sin(ang)    % matriz de rotacion
        sin(ang)   cos(ang) ];

%% Dibujar de deformada
figure(2)
hold on
pos = T*[ s + esc_def*u; esc_def*v ];

xx = pos(X,:) + x1;
yy = pos(Y,:) + y1;

plot([x1 x2], [y1 y2], 'b-', xx, yy, 'r-','LineWidth',2);

%% Dibujar los diagramas de fuerza axial 
figure(3)
hold on
pos = T*[ s; esc_faxial*axial ]; % escalamiento del diagrama

ss = pos(X,:) + x1;
aa = pos(Y,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 aa y2], 'r-','LineWidth',2);
text(ss(1),   aa(1),   num2str(-qe(X1)));
text(ss(end), aa(end), num2str(+qe(X2)));

%% Dibujar los diagramas de fuerza cortante
figure(4)
hold on
pos = T*[ s; esc_V*V ]; % escalamiento del diagrama

ss = pos(X,:) + x1;
vv = pos(Y,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 vv y2], 'r-','LineWidth',2);
text(ss(1),   vv(1),   num2str(+qe(Y1)));
text(ss(end), vv(end), num2str(-qe(Y2)));

%% Dibujar los diagramas de momento flector
figure(5)
hold on
pos = T*[ s; esc_M*M ]; % escalamiento del diagrama

ss = pos(X,:) + x1;
mm = pos(Y,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 mm y2], 'r-','LineWidth',2);
text(ss(1),   mm(1),   num2str(-qe(M1)));
text(ss(end), mm(end), num2str(+qe(M2)));
[minM,idminM] = min(M); text(ss(idminM), mm(idminM), num2str(minM));
[maxM,idmaxM] = max(M); text(ss(idmaxM), mm(idmaxM), num2str(maxM));

%% ------------------------------------------------------------------------
   function dydx = ecuacion_diferencial(x,y)
      % aqui se implementa la ecuacion diferencial para vigas de material
      % homogeneo y seccion transversal constante (A, E, I, qx, qy las 
      % provee la funcion exterior)
      %      d^4 v(x)
      % E I ---------- = q(x)
      %        dx^4
      %
      %      d^2 u(x)
      % A E ---------- = -b(x)
      %        dx^2

      dydx = zeros(6,1);
      %         y(1)          = v
      dydx(1) = y(2);       % = theta
      dydx(2) = y(3)/(E*I); % = M/(EI)
      dydx(3) = y(4);       % = V
      dydx(4) = qyloc(x);   % = qyloc
      dydx(5) = y(6)/(A*E); % = u
      dydx(6) = -qxloc(x);  % = faxial
   end

%% ------------------------------------------------------------------------

   function res = condiciones_de_apoyo(YL,YR)
      % condiciones de apoyo (cond. de frontera de la ecuacion diferencial)
      u1  = 1; v1 = 2; t1 = 3; u2 = 4; v2 = 5; t2   = 6;
      v_  = 1; t_ = 2; M_ = 3; V_ = 4; u_ = 5; fax_ = 6;
      res = [ % YL: apoyo izquierdo (LEFT), YR: apoyo derecho (RIGHT)
              YL(u_) - ae(u1)          % uloc(0)     = u1
              YL(v_) - ae(v1)          % vloc(0)     = v1
              YL(t_) - ae(t1)          % thetaloc(0) = t1
              YR(u_) - ae(u2)          % uloc(L)     = u2
              YR(v_) - ae(v2)          % vloc(L)     = v2
              YR(t_) - ae(t2) ];       % thetaloc(L) = t2
   end
end
%% bye, bye!
