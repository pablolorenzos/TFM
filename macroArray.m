clear; clc; close all;

%% Parámetros iniciales
f0 = 1.3e9;
lmbda = physconst('lightspeed')/f0;
fprintf('λ = %.4f m\n', lmbda);

h = 1000; % Altura en Km
elevacion = deg2rad(20); % Elevación en radianes (20°)
r = 6378; % Radio de la Tierra en Km
epsilon = 4.6;
lmbda_eff = lmbda/sqrt(epsilon);
fprintf('λ efectiva = %.4f m\n', lmbda_eff);

slant_range = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));
fprintf('slant range = %.4f m\n', slant_range);

%% Patrón deseado (half-cosine)
% Parámetros
R_E       = 6378;        % km, radio de la Tierra
slant_ref = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));        % km, rango en los picos

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
% peakAngle         = 54;
peakAngle = atan2(R_E*sin(phi_max), (R_E+h) - R_E*cos(phi_max)) * 180/pi;
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
frameLen = 9;      % longitud de la ventana (impar)
Beam_dB_sg = sgolayfilt( Beam_dB, polOrder, frameLen );

Beam_sg = 10.^( Beam_dB_sg/20 );

Beam_d = Beam_sg;

figure;
plot(azimuth, mag2db(Beam_d), 'b', 'LineWidth',2);
xlabel('\theta (°)'); ylabel('Ganancia (dB)');
title('Patrón deseado (half-cosine)');
grid on; ylim([-15 20]);

%% Definición del Macro Array
numRings = 15;  % Ahora 20 anillos
% Radios de cada anillo (múltiplos de 0.1*λ_eff)
ringRadii = (1:numRings) * 0.1 * lmbda_eff;
% Número de elementos en cada anillo: primer anillo 6, y se añaden 2 elementos por cada anillo superior
ringElements = 6 + 2*(0:numRings-1);
N_total = sum(ringElements);
fprintf('Total de elementos: %d\n', N_total);

% Cálculo de las posiciones de los elementos y asignación de índices de anillo
elementPos = [];
ringIndex = [];  % Vector que indica a qué anillo pertenece cada elemento
for i = 1:numRings
    N_i = ringElements(i);
    angles = (0:N_i-1) * (360/N_i);  % Ángulos en grados, uniformemente espaciados
    x_ring = ringRadii(i) * cosd(angles);
    y_ring = ringRadii(i) * sind(angles);
    elementPos = [elementPos, [x_ring; y_ring]];
    ringIndex = [ringIndex, repmat(i, 1, N_i)];
end

%% Definición de la máscara binaria inicial (on/off) para cada elemento
% Al quitar la restricción, inicializamos la máscara con todos los elementos encendidos.
b0 = ones(N_total,1);

%% Optimización con GA (desfases por anillo, pesos complejos y máscara on/off)
% El vector de optimización contendrá:
%   - numRings variables para los desfases de cada anillo (en grados)
%   - 2*N_total variables para los pesos complejos (parte real e imaginaria)
%   - N_total variables binarias para la activación (on/off) de cada elemento
dim = numRings + 2 * N_total + N_total;  % numRings + 3*N_total

% Vector inicial:
%   - Desfases: inicializados a 0
%   - Pesos: parte real = 1 y parte imaginaria = 0
%   - Máscara: según b0 (todos encendidos)
x0 = [zeros(numRings,1); ones(N_total,1); zeros(N_total,1); b0];

% Límites:
lb = zeros(dim,1);
ub = zeros(dim,1);
% Límites para desfases (en grados)
lb(1:numRings) = 0;
ub(1:numRings) = 0;
% Límites para la parte real de los pesos
lb(numRings+1 : numRings+N_total) = -5;
ub(numRings+1 : numRings+N_total) = 5;
% Límites para la parte imaginaria de los pesos
lb(numRings+N_total+1 : numRings+2*N_total) = -8;
ub(numRings+N_total+1 : numRings+2*N_total) = 8;
% Límites para la máscara binaria: deben estar entre 0 y 1
lb(numRings+2*N_total+1 : end) = 0;
ub(numRings+2*N_total+1 : end) = 1;

% -------------------------------------------------------------
% A = [zeros(1, numRings+2*N_total),  ones(1, N_total)];
% b = N_total/2;
% -------------------------------------------------------------

% Especificar que las últimas N_total variables son enteras (0 o 1)
intcon = (numRings + 2*N_total+1) : dim;
objfun = @(x) macroArrayCostGA(x, elementPos, ringIndex, azimuth, Beam_d, lmbda_eff);

options = optimoptions('ga', ...
    'Display','iter', ...
    'FunctionTolerance',1e-16, ...
    'ConstraintTolerance',1e-10, ...
    'CrossoverFraction',0.7, ...
    'PopulationSize',300, ...
    'MaxGenerations',2, ...
    'PlotFcn', @gaplotbestf);

[x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [], options);

%% Extraer y aplicar las variables optimizadas
% Desfase para cada anillo:
phase_offsets_opt = x_opt(1:numRings);  % en grados
% Pesos complejos optimizados:
w_re_opt = x_opt(numRings+1 : numRings+N_total).';
w_im_opt = x_opt(numRings+N_total+1 : numRings+2*N_total).';
w_raw_opt = w_re_opt + 1i*w_im_opt;
% Máscara binaria optimizada:
b_opt = x_opt(numRings+2*N_total+1 : end);

% Aplicar el desfase correspondiente a cada elemento según su anillo y la máscara
w_opt = zeros(1, N_total);
for idx = 1:N_total
    ring_num = ringIndex(idx);
    w_opt(idx) = b_opt(idx) * w_raw_opt(idx) * exp(1i*deg2rad(phase_offsets_opt(ring_num)));
end

%% Síntesis del patrón utilizando la configuración optimizada
stvmat_opt = zeros(N_total, length(azimuth));
for i = 1:length(azimuth)
    phi = azimuth(i);
    stvmat_opt(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda_eff );
end

pattern_synth_opt = abs(w_opt * stvmat_opt);

figure;
plot(azimuth, mag2db(Beam_d), 'b', 'LineWidth',2);
hold on;
plot(azimuth, mag2db(pattern_synth_opt), 'r--', 'LineWidth',2);
legend('Desired','Synthesized','Location','Best');
xlabel('\theta (°)'); ylabel('Gain (dB)');
title('Optimized pattern with GA');
grid on; ylim([-15 20]);

%% Visualización de la geometría del Macro Array (distinguiendo encendidos y apagados)
figure;
% Determinar índices de elementos activados y desactivados según la máscara optimizada
active_indices = find(b_opt >= 0.5);
inactive_indices = find(b_opt < 0.5);

% Graficar elementos encendidos con marcadores llenos
scatter(elementPos(1,active_indices), elementPos(2,active_indices), 60, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Encendidos');
hold on;
% Graficar elementos apagados con marcadores no llenos (por ejemplo, círculos rojos)
scatter(elementPos(1,inactive_indices), elementPos(2,inactive_indices), 60, 'o', 'MarkerEdgeColor', 'r', 'LineWidth',1.5, 'DisplayName', 'Apagados');

% Dibujar los contornos de cada anillo
thetaC = linspace(0, 2*pi, 200);
for i = 1:numRings
    plot(ringRadii(i)*cos(thetaC), ringRadii(i)*sin(thetaC), 'k--');
end

axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)');
title('Positions in the Macro Array (ON vs OFF)');
legend('show');
