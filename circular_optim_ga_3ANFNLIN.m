
%%
% AQUÍ VAMOS A INTRODUCIR UN TERCER ANILLO
% CON FUNCION NO LINEAL
% SIN ELEMENTO CENTRAL
%%

clear; clc; close all;

%%
f0 = 1.3e9;
lmbda = physconst('lightspeed')/f0;
fprintf('λ = %.4f m\n', lmbda);

h = 1000; % altura en Km
elevacion = deg2rad(20); % en grados
r = 6378; % radio tierra en Km
epsilon = 4.6;
lmbda_eff = lmbda/sqrt(epsilon);
fprintf('λ = %.4f m\n', lmbda_eff);

slant_range = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));
fprintf('slant range = %.4f m\n', slant_range);

%%
% Parámetros
R_E       = 6378;        % km, radio de la Tierra
slant_ref = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));        % km

% ecuación de ley de cosenos para φ_max
phi_max = acos( ((R_E+h)^2 + R_E^2 - slant_ref^2) ...
              /(2*R_E*(R_E+h)) );
phi_max_deg = phi_max*180/pi;   % ≃ 15.69°

N       = 2001;
phi      = linspace(-phi_max, +phi_max, N);  % en rad
% slant-range exacto en cada φ
d_phi    = sqrt((R_E+h)^2 + R_E^2 ...
              - 2*R_E*(R_E+h)*cos(phi));
% lo convertimos a off–axis θ real
theta_iso = atan2( R_E*sin(phi), (R_E+h)-R_E*cos(phi) );
% ganancia relativa (0 dB en φ=±φ_max)
GdB_iso   = 20*log10( d_phi ./ slant_ref );


% ángulos diagrama
peakAngle         = 54;
% peakAngle = atan2(R_E*sin(phi_max), (R_E+h) - R_E*cos(phi_max)) * 180/pi;
intermediateAngle = 70;  
outerAngle        = 100;  

% niveles
peak_dB         = 15;
intermediate_dB = -5;
outer_dB        = -5;

delta = 4;    % grados de “ancho” extra del pico

azimuth = -100:100;
Beam_dB = zeros(size(azimuth));

for k = 1:numel(azimuth)
  th = azimuth(k);  x = abs(th);
  if x <= (peakAngle-delta)
    % valle isoflux
    dBrel      = interp1(theta_iso*180/pi, GdB_iso, th, 'linear');
    Beam_dB(k) = peak_dB + dBrel;

  elseif x <= (peakAngle + delta)
   % Plateau plano en peak_dB
   Beam_dB(k) = peak_dB;

  elseif x <= intermediateAngle
    % half‐cosine de +15 dB → –5 dB
    t = (x-peakAngle)/(intermediateAngle-peakAngle);
    Beam_dB(k) = peak_dB + (intermediate_dB-peak_dB)*(0.5-0.5*cos(pi*t));

  elseif x <= outerAngle
    % half‐cosine de –5 dB → –5 dB (plano)
    t = (x-intermediateAngle)/(outerAngle-intermediateAngle);
    Beam_dB(k) = intermediate_dB + (outer_dB-intermediate_dB)*(0.5-0.5*cos(pi*t));

  else
    Beam_dB(k) = outer_dB;
  end
end

polOrder = 3;      % orden del polinomio
frameLen = 9;      % longitud de la ventana (impar, p.ej. 7, 9, 11)
Beam_dB_sg = sgolayfilt( Beam_dB, polOrder, frameLen );

Beam_sg = 10.^( Beam_dB_sg/20 );


Beam_d = Beam_sg;


figure;
plot(azimuth, mag2db(Beam_d), 'b', 'LineWidth',2);
xlabel('\theta (°)'); ylabel('Ganancia (dB)');
title('Patrón deseado (half-cosine)');
grid on; ylim([-15 20]);

%%
N1 = 4;% Elementos en anillo 1
N2 = 6;% Elementos en anillo 2
N3 = 8;% Elementos en anillo 3
N_total = N1 + N2 + N3;

%%
objfun = @(x) CostGA3ANFNLIN(x, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle);

%%
dim = 5 + 2*N_total;  % 5 variables geométricas + 2*16 = 35

x0 = zeros(dim,1);
x0(1) = 0.5*lmbda_eff; % r1 
x0(2) = 1.0*lmbda_eff; % r2
x0(3) = 1.5*lmbda_eff; % r3
x0(4) = 0; % desfaseRel inicial
x0(5) = 0; % desfaseRel1 inicial
x0(6 : 5+N_total) = 1; % parte real a 1
x0(6+N_total : 5+2*N_total) = 0;  % parte imaginaria a 0

%%
% Definir límites:
lb = zeros(dim,1);
ub = zeros(dim,1);

% Límites para r1
lb(1) = 0.45*lmbda_eff;
ub(1) = 0.55*lmbda_eff;
% Límites para r2
lb(2) = 0.9*lmbda_eff;
ub(2) = 1.1*lmbda_eff;
% Límites para r3
lb(3) = 1.45*lmbda_eff;
ub(3) = 2.05*lmbda_eff;
% Límites para el desfase relativo (en grados)
lb(4) = -90;
ub(4) = 90;
% Límites para el desfase relativo1 (en grados)
lb(5) = -90;
ub(5) = 90;
% Límites para la parte real de los pesos
lb(6 : 5+N_total) = -50;
ub(6 : 5+N_total) = 50;
% Límites para la parte imaginaria de los pesos
lb(6+N_total : 5+2*N_total) = -50;
ub(6+N_total : 5+2*N_total) = 50;

%%

mask_dB = -5;  % Umbral de la máscara en dB
nonlcon = @(x) my3ANFNLin(x, N1, N2, N3, azimuth, lmbda, intermediateAngle, mask_dB);


% Arrays donde se irán guardando ambos costes
centralHistory = [];
externalHistory = [];

ofcn = @(options,state,flag) costMonitor(...
    options, state, flag, ...
    N1, N2, N3, ...
    azimuth, Beam_d, ...
    lmbda_eff, intermediateAngle);

options = optimoptions('ga', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'ConstraintTolerance',1e-3, ...
    'CrossoverFraction',0.7, ...
    'PopulationSize',300, ...
    'MaxGenerations',800, ...
    'PlotFcn','gaplotbestf', ...
    'OutputFcn',ofcn);

[x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, nonlcon, options);

% dos historiales
figure; hold on;
gens = 0:length(centralHistory)-1;
plot(gens, centralHistory, '-b','LineWidth',2);
plot(gens, externalHistory, '-r','LineWidth',2);
legend('Central cost','External cost','Location','best');
xlabel('Generation'); ylabel('Cost');
title('Evolution Cost Central vs External');
grid on;


r1_opt         = x_opt(1);
r2_opt         = x_opt(2);
r3_opt         = x_opt(3);
desfaseRel_opt = x_opt(4);
desfaseRel1_opt = x_opt(5);

w_re_opt = x_opt(6 : 5+N_total).';
w_im_opt = x_opt(6+N_total : 5+2*N_total).';
w_opt = w_re_opt + 1i*w_im_opt;


%% 
% centro:
anglesRing1_opt = (0:N1-1)*(360/N1);  
anglesRing2_opt = (0:N2-1)*(360/N2) + desfaseRel_opt;
anglesRing3_opt = (0:N3-1)*(360/N3) + desfaseRel1_opt;

posRing1_opt = [r1_opt*cosd(anglesRing1_opt); r1_opt*sind(anglesRing1_opt)];
posRing2_opt = [r2_opt*cosd(anglesRing2_opt); r2_opt*sind(anglesRing2_opt)];
posRing3_opt = [r3_opt*cosd(anglesRing3_opt); r3_opt*sind(anglesRing3_opt)];
elementPos_opt = [posRing1_opt, posRing2_opt, posRing3_opt];

% matriz de apuntamiento
stvmat_opt = zeros(N_total, length(azimuth));
for i = 1:length(azimuth)
    phi = azimuth(i);
    stvmat_opt(:, i) = exp(1i * 2*pi * ( elementPos_opt(1,:)'*cosd(phi) + elementPos_opt(2,:)'*sind(phi) ) / lmbda_eff );
end

pattern_synth_opt = abs(w_opt'*stvmat_opt);

fprintf('Radio first ring = %.4f m\n', r1_opt);
fprintf('Radio second ring = %.4f m\n', r2_opt);
fprintf('Radio third ring = %.4f m\n', r3_opt);
fprintf('Phase shift second ring = %.4f º\n', desfaseRel_opt);
fprintf('Phase shift third ring = %.4f º\n', desfaseRel1_opt);

% patrón deseado vs. sintetizado
figure;
plot(azimuth, mag2db([Beam_d; pattern_synth_opt])', 'LineWidth',2);
legend('Desired','Synthesized','Location','Best');
xlabel('\theta (°)'); ylabel('Gain (dB)');
title('Optimized Pattern with GA');
grid on; ylim([-15 20]);

%%
figure;
scatter(elementPos_opt(1,:), elementPos_opt(2,:), 60, 'filled', 'b');
hold on;
thetaC = linspace(0,2*pi,200);
plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--');  % círculo del anillo 1
plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--');  % círculo del anillo 2
plot(r3_opt*cos(thetaC), r3_opt*sin(thetaC), 'k--');  % círculo del anillo 3
axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)');
title('Optimized Array Positions (Circular, 16 elements)');

%% —– Gráfico numerado con leyenda de 20 elementos —–
figure; hold on;

hC1 = plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
hC2 = plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
hC3 = plot(r3_opt*cos(thetaC), r3_opt*sin(thetaC), 'k--', 'HandleVisibility','off');

% 2) Para cada elemento creamos su propio scatter, con DisplayName = "n: módulo ∠ fase"
hElem = gobjects(N_total,1);
for k = 1:N_total
    xk = elementPos_opt(1,k);
    yk = elementPos_opt(2,k);
    lbl = sprintf('%d: %.2f ∠ %.1f°', k, abs(w_opt(k)), angle(w_opt(k))*180/pi);

    % scatter individual con DisplayName
    hElem(k) = scatter(xk, yk, 80, 'b', 'filled', ...
                       'DisplayName', lbl);

    % el número en rojo encima
    text(xk, yk, sprintf('%d',k), ...
        'FontSize',8,'FontWeight','bold','Color','r', ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)');
title('Array Elements Numbered and Weights');


legend(hElem, 'Location','eastoutside', 'FontSize',8);



%%
function [state,options,optchanged] = costMonitor( ...
        options, state, flag, ...
        N1, N2, N3, ...
        azimuth, Beam_d, ...
        lmbda_eff, intermediateAngle)
    optchanged = false;
    persistent centralHistory externalHistory
    switch flag
      case 'init'
        centralHistory = [];
        externalHistory = [];
      case 'iter'
        % 1) mejor individuo
        [~, idx] = min(state.Score);
        xBest = state.Population(idx, :).';
        % 2) evalúo componentes
        [~, cC, cE] = CostGA3ANFNLIN( xBest, ...
                          N1, N2, N3, ...
                          azimuth, Beam_d, ...
                          lmbda_eff, intermediateAngle);
        centralHistory(end+1)  = cC;
        externalHistory(end+1) = cE;
      case 'done'
        % vuelco al workspace
        assignin('base','centralHistory', centralHistory);
        assignin('base','externalHistory', externalHistory);
    end
end

%%
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------

%%
% AQUÍ VAMOS A INTRODUCIR UN TERCER ANILLO
% CON FUNCIÓN NO LINEAL
% CON ELEMENTO CENTRAL
%%

% clear; clc; close all;
% 
% %%
% f0 = 1.3e9;
% lmbda = physconst('lightspeed')/f0;
% fprintf('λ = %.4f m\n', lmbda);
% 
% h = 1000; % altura en Km
% elevacion = deg2rad(35); % en grados
% r = 6378; % radio tierra en Km
% epsilon = 4.6;
% lmbda_eff = lmbda/sqrt(epsilon);
% fprintf('λ = %.4f m\n', lmbda_eff);
% 
% slant_range = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));
% fprintf('slant range = %.4f m\n', slant_range);
% 
% %% 
% 
% % Parámetros
% R_E       = 6378;        % km, radio de la Tierra
% slant_ref = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));        % km
% 
% % ecuación de ley de cosenos para φ_max
% phi_max = acos( ((R_E+h)^2 + R_E^2 - slant_ref^2) ...
%               /(2*R_E*(R_E+h)) );
% phi_max_deg = phi_max*180/pi;   % ≃ 15.69°
% 
% N       = 2001;
% phi      = linspace(-phi_max, +phi_max, N);  % en rad
% % slant-range exacto en cada φ
% d_phi    = sqrt((R_E+h)^2 + R_E^2 ...
%               - 2*R_E*(R_E+h)*cos(phi));
% % lo convertimos a off–axis θ real
% theta_iso = atan2( R_E*sin(phi), (R_E+h)-R_E*cos(phi) );
% % ganancia relativa (0 dB en φ=±φ_max)
% GdB_iso   = 20*log10( d_phi ./ slant_ref );
% 
% 
% % ángulos diagrama
% peakAngle         = 54;
% % peakAngle = atan2(R_E*sin(phi_max), (R_E+h) - R_E*cos(phi_max)) * 180/pi;
% intermediateAngle = 70;  
% outerAngle        = 100;  
% 
% % niveles
% peak_dB         = 15;
% intermediate_dB = -5;
% outer_dB        = -5;
% 
% delta = 4;    % grados de “ancho” extra del pico
% 
% azimuth = -100:100;
% Beam_dB = zeros(size(azimuth));
% 
% for k = 1:numel(azimuth)
%   th = azimuth(k);  x = abs(th);
%   if x <= (peakAngle-delta)
%     % valle isoflux
%     dBrel      = interp1(theta_iso*180/pi, GdB_iso, th, 'linear');
%     Beam_dB(k) = peak_dB + dBrel;
% 
%   elseif x <= (peakAngle + delta)
%    % Plateau plano en peak_dB
%    Beam_dB(k) = peak_dB;
% 
%   elseif x <= intermediateAngle
%     % half‐cosine de +15 dB → –5 dB
%     t = (x-peakAngle)/(intermediateAngle-peakAngle);
%     Beam_dB(k) = peak_dB + (intermediate_dB-peak_dB)*(0.5-0.5*cos(pi*t));
% 
%   elseif x <= outerAngle
%     % half‐cosine de –5 dB → –5 dB (plano)
%     t = (x-intermediateAngle)/(outerAngle-intermediateAngle);
%     Beam_dB(k) = intermediate_dB + (outer_dB-intermediate_dB)*(0.5-0.5*cos(pi*t));
% 
%   else
%     Beam_dB(k) = outer_dB;
%   end
% end
% 
% polOrder = 3;      % orden del polinomio
% frameLen = 9;      % longitud de la ventana (impar)
% Beam_dB_sg = sgolayfilt( Beam_dB, polOrder, frameLen );
% 
% Beam_sg = 10.^( Beam_dB_sg/20 );
% 
% Beam_d = Beam_sg;
% 
% figure;
% plot(azimuth, mag2db(Beam_d), 'b', 'LineWidth',2);
% xlabel('\theta (°)'); ylabel('Ganancia (dB)');
% title('Patrón deseado (half-cosine)');
% grid on; ylim([-15 20]);
% 
% %%
% N0 = 1; % Elemento central
% N1 = 4;% Elementos en anillo 1
% N2 = 6;% Elementos en anillo 2
% N3 = 8;% Elementos en anillo 3
% N_total = N0 + N1 + N2 + N3;
% 
% 
% %%
% objfun = @(x) CostGA3ANFNLIN(x, N0, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle);
% 
% %%
% dim = 5 + 2*N_total;  % 5 variables geométricas + 2*16 = 35
% 
% x0 = zeros(dim,1);
% x0(1) = 0.5*lmbda_eff; % r1 
% x0(2) = 1.0*lmbda_eff; % r2
% x0(3) = 1.5*lmbda_eff; % r3
% x0(4) = 0; % desfaseRel inicial
% x0(5) = 0; % desfaseRel1 inicial
% x0(6 : 5+N_total) = 1; % parte real a 1
% x0(6+N_total : 5+2*N_total) = 0;  % parte imaginaria a 0
% 
% %%
% % Definir límites:
% lb = zeros(dim,1);
% ub = zeros(dim,1);
% 
% % Límites para r1
% lb(1) = 0.45*lmbda_eff;
% ub(1) = 0.55*lmbda_eff;
% % Límites para r2
% lb(2) = 0.9*lmbda_eff;
% ub(2) = 1.1*lmbda_eff;
% % Límites para r3
% lb(3) = 1.45*lmbda_eff;
% ub(3) = 2.05*lmbda_eff;
% % Límites para el desfase relativo (en grados)
% lb(4) = -90;
% ub(4) = 90;
% % Límites para el desfase relativo1 (en grados)
% lb(5) = -90;
% ub(5) = 90;
% % Límites para la parte real de los pesos
% lb(6 : 5+N_total) = -50;
% ub(6 : 5+N_total) = 50;
% % Límites para la parte imaginaria de los pesos
% lb(6+N_total : 5+2*N_total) = -50;
% ub(6+N_total : 5+2*N_total) = 50;
% 
% %%
% mask_dB = -5;  % Umbral de la máscara en dB
% nonlcon = @(x) my3ANFNLin(x, N0, N1, N2, N3, azimuth, lmbda, intermediateAngle, mask_dB);
% 
% 
% % Arrays donde se irán guardando ambos costes
% centralHistory = [];
% externalHistory = [];
% 
% ofcn = @(options,state,flag) costMonitor(...
%     options, state, flag, ...
%     N0, N1, N2, N3, ...
%     azimuth, Beam_d, ...
%     lmbda_eff, intermediateAngle);
% 
% options = optimoptions('ga', ...
%     'Display','iter', ...
%     'FunctionTolerance',1e-8, ...
%     'ConstraintTolerance',1e-3, ...
%     'CrossoverFraction',0.4, ...
%     'PopulationSize',300, ...
%     'MaxGenerations',800, ...
%     'PlotFcn','gaplotbestf', ...
%     'OutputFcn',ofcn);
% 
% [x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, nonlcon, options);
%
% % dos historiales
% figure; hold on;
% gens = 0:length(centralHistory)-1;
% plot(gens, centralHistory, '-b','LineWidth',2);
% plot(gens, externalHistory, '-r','LineWidth',2);
% legend('Central cost','External cost','Location','best');
% xlabel('Generation'); ylabel('Cost');
% title('Evolution Cost Central vs External');
% grid on;
% 
% 
% r1_opt         = x_opt(1);
% r2_opt         = x_opt(2);
% r3_opt         = x_opt(3);
% desfaseRel_opt = x_opt(4);
% desfaseRel1_opt = x_opt(5);
% 
% w_re_opt = x_opt(6 : 5+N_total).';
% w_im_opt = x_opt(6+N_total : 5+2*N_total).';
% w_opt = w_re_opt + 1i*w_im_opt;
% 
% 
% %% 
% % centro:
% pos0 = [0;0];
% 
% anglesRing1_opt = (0:N1-1)*(360/N1);  
% anglesRing2_opt = (0:N2-1)*(360/N2) + desfaseRel_opt;
% anglesRing3_opt = (0:N3-1)*(360/N3) + desfaseRel1_opt;
% 
% posRing1_opt = [r1_opt*cosd(anglesRing1_opt); r1_opt*sind(anglesRing1_opt)];
% posRing2_opt = [r2_opt*cosd(anglesRing2_opt); r2_opt*sind(anglesRing2_opt)];
% posRing3_opt = [r3_opt*cosd(anglesRing3_opt); r3_opt*sind(anglesRing3_opt)];
% elementPos_opt = [pos0, posRing1_opt, posRing2_opt, posRing3_opt];
% 
% 
% % matriz de apuntamiento
% stvmat_opt = zeros(N_total, length(azimuth));
% for i = 1:length(azimuth)
%     phi = azimuth(i);
%     stvmat_opt(:, i) = exp(1i * 2*pi * ( elementPos_opt(1,:)'*cosd(phi) + elementPos_opt(2,:)'*sind(phi) ) / lmbda_eff );
% end
% 
% pattern_synth_opt = abs(w_opt'*stvmat_opt);
% 
% fprintf('Radio first ring = %.4f m\n', r1_opt);
% fprintf('Radio second ring = %.4f m\n', r2_opt);
% fprintf('Radio third ring = %.4f m\n', r3_opt);
% fprintf('Phase shift second ring = %.4f º\n', desfaseRel_opt);
% fprintf('Phase shift third ring = %.4f º\n', desfaseRel1_opt);
% 
% 
% % patrón deseado vs. sintetizado
% figure;
% plot(azimuth, mag2db([Beam_d; pattern_synth_opt])', 'LineWidth',2);
% legend('Desired','Synthesized','Location','Best');
% xlabel('\theta (°)'); ylabel('Gain (dB)');
% title('Optimized Pattern with GA');
% grid on; ylim([-15 20]);
% 
% %%
% figure;
% scatter(elementPos_opt(1,:), elementPos_opt(2,:), 60, 'filled', 'b');
% hold on;
% thetaC = linspace(0,2*pi,200);
% plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--');  % círculo del anillo 1
% plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--');  % círculo del anillo 2
% plot(r3_opt*cos(thetaC), r3_opt*sin(thetaC), 'k--');  % círculo del anillo 3
% axis equal; grid on;
% xlabel('X (m)'); ylabel('Y (m)');
% title('Optimized Array Positions (Circular, 16 elements)');
% 
% %%
% figure; hold on;
% 
% hC1 = plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% hC2 = plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% hC3 = plot(r3_opt*cos(thetaC), r3_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% 
% % 2) Para cada elemento creamos su propio scatter, con DisplayName = "n: módulo ∠ fase"
% hElem = gobjects(N_total,1);
% for k = 1:N_total
%     xk = elementPos_opt(1,k);
%     yk = elementPos_opt(2,k);
%     lbl = sprintf('%d: %.2f ∠ %.1f°', k, abs(w_opt(k)), angle(w_opt(k))*180/pi);
% 
%     % scatter individual con DisplayName
%     hElem(k) = scatter(xk, yk, 80, 'b', 'filled', ...
%                        'DisplayName', lbl);
% 
%     % opcional: el número en rojo encima
%     text(xk, yk, sprintf('%d',k), ...
%         'FontSize',8,'FontWeight','bold','Color','r', ...
%         'HorizontalAlignment','center','VerticalAlignment','middle');
% end
% 
% axis equal; grid on;
% xlabel('X (m)'); ylabel('Y (m)');
% title('Array Elements Numbered and Weights');
% 
% legend(hElem, 'Location','eastoutside', 'FontSize',8);
% 
% 
% 
% %%
% function [state,options,optchanged] = costMonitor( ...
%         options, state, flag, ...
%         N0, N1, N2, N3, ...
%         azimuth, Beam_d, ...
%         lmbda_eff, intermediateAngle)
%     optchanged = false;
%     persistent centralHistory externalHistory
%     switch flag
%       case 'init'
%         centralHistory = [];
%         externalHistory = [];
%       case 'iter'
%         % 1) mejor individuo
%         [~, idx] = min(state.Score);
%         xBest = state.Population(idx, :).';
%         % 2) evalúo componentes
%         [~, cC, cE] = CostGA3ANFNLIN( xBest, ...
%                           N0, N1, N2, N3, ...
%                           azimuth, Beam_d, ...
%                           lmbda_eff, intermediateAngle);
%         centralHistory(end+1)  = cC;
%         externalHistory(end+1) = cE;
%       case 'done'
%         % vuelco al workspace
%         assignin('base','centralHistory', centralHistory);
%         assignin('base','externalHistory', externalHistory);
%     end
% end