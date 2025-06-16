
%%
%%
%%
%%
% RESTRICCIONES 1/2/3/4/5/6/7
clear; clc; close all;
%% Script optimización de Macro Array sin desfases por anillo
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
% peakAngle         = 54;
peakAngle = atan2(R_E*sin(phi_max), (R_E+h) - R_E*cos(phi_max)) * 180/pi;
intermediateAngle = 70;  
outerAngle        = 100;  

% niveles
peak_dB         = 15;
intermediate_dB = -5;
outer_dB        = -5;

delta = 4;    % grados de “ancho” extra del pico

azimuth = -90:90;
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
%% Definición del Macro Array
numRings = 11;  % Ahora 7 anillos
% Radios: primer anillo en 0.5·λ_eff, separación 0.25·λ_eff
ringRadii = (0.5 + (0:numRings-1)*0.25) * lmbda_eff;
% Número de elementos por anillo: [4,8,16,16,16,16,32]
% ringElements = [4, 8, repmat(16,1,numRings-3), 32];
ringElements = [4, 0, 8, 0, 8, 0, 8, 0, 16, 0, 16];
N_total = sum(ringElements);
fprintf('Total de elementos: %d\n', N_total);

nBits      = 6;
Nlevels    = 2^nBits;
phaseStep  = 360 / Nlevels;


% Cálculo de posiciones: anillos concéntricos, elementos radiales
elementPos = [];
ringIndex = [];
for i = 1:numRings
    N_i = ringElements(i);
    angles = (0:N_i-1) * (360 / N_i);  % grados uniformes
    x_ring = ringRadii(i) * cosd(angles);
    y_ring = ringRadii(i) * sind(angles);
    elementPos = [elementPos, [x_ring; y_ring]];
    ringIndex = [ringIndex, repmat(i, 1, N_i)];
end

%% Máscara binaria inicial (todos encendidos)
b0 = ones(N_total,1);

%% Optimización con GA (pesos complejos, máscara on/off, activación de anillos y fases cuantizadas por elemento)
% Restricción elementos activos: mínimo y máximo
Mmin_el = 13;
Mmax_el = 25;
% Restricción anillos activos: mínimo y máximo
Mmin_an = 3;
Mmax_an = 4;
% Variables de fase discretas por elemento (6 bits → 64 niveles)
% Número de variables de fase: N_total

% Dimensión del vector de optimización:
% [Re(w) (N_total); Im(w) (N_total); b (N_total); z (numRings); n (N_total)]
dim = 3*N_total + numRings + N_total;

% Vector inicial x0:
x0 = [ ...
    ones(N_total,1);      % Re(w)
    zeros(N_total,1);     % Im(w)
    ones(N_total,1);      % b (todos encendidos inicial)
    ones(numRings,1);     % z (todos anillos activos inicial)
    zeros(N_total,1)      % n (fase mínima = -180°)
];

% Cotas:
lb = [ ...
    -ones(N_total,1);     % Re(w) >= -1
    -ones(N_total,1);   % Im(w) >= -1
    zeros(N_total,1);     % b >= 0
    zeros(numRings,1);    % z >= 0
    zeros(N_total,1)      % n >= 0
];
ub = [ ...
    2*ones(N_total,1);    % Re(w) <= 2
    2*ones(N_total,1);    % Im(w) <= 2
    ones(N_total,1);      % b <= 1
    ones(numRings,1);     % z <= 1
    (Nlevels-1)*ones(N_total,1)    % n <= 63 (índices de fase)
];

% Índices de variables enteras:
intcon = [ ...
    (2*N_total+1):(3*N_total), ...           % b
    (3*N_total+1):(3*N_total+numRings), ...   % z
    (3*N_total+numRings+1):dim ...            % n
];

% Construir Aineq y bineq
% 1) sum(b) >= Mmin_el  -> -sum(b) <= -Mmin_el
%                Re(Im)   b          z              fases
A1 = [zeros(1,2*N_total), -ones(1,N_total), zeros(1,numRings), zeros(1,N_total)];
b1 = -Mmin_el;

% 2) sum(b) <= Mmax_el
A2 = [zeros(1,2*N_total),  ones(1,N_total), zeros(1,numRings), zeros(1,N_total)];
b2 =  Mmax_el;

% 3) sum(z) >= Mmin_an -> -sum(z) <= -Mmin_an
A3 = [zeros(1,2*N_total), zeros(1,N_total), -ones(1,numRings), zeros(1,N_total)];
b3 = -Mmin_an;

% 4) sum(z) <= Mmax_an
A4 = [zeros(1,2*N_total), zeros(1,N_total),  ones(1,numRings), zeros(1,N_total)];
b4 =  Mmax_an;

% 5) enlace de activación: b_j <= z_ring(j)
A5 = zeros(N_total, dim);
b5 = zeros(N_total,1);
for j = 1:N_total
    ring_j = ringIndex(j);
    % b_j - z_k <= 0
    A5(j, 2*N_total + j)      = 1;    % b_j
    A5(j, 3*N_total + ring_j) = -1;   % -z_ring(j)
end

% 6) Restricción de separación mínima entre elementos activos en cada anillo
threshold = 0.25 * lmbda_eff;  % ajuste de separación deseada (p.ej. 0.4·λ_eff)
A6 = [];
b6 = [];
for ring_i = 1:numRings
    N_i = ringElements(ring_i);
    radius_i = ringRadii(ring_i);
    % Longitud de arco entre vecinos (un paso angular)
    arcStep = 2*pi * radius_i / N_i;
    % Número de pasos contiguos que quedan dentro del threshold
    maxSteps = floor(threshold / arcStep);
    if maxSteps >= 1
        idx_ring = find(ringIndex == ring_i);
        for m = 1:N_i
            for step = 1:maxSteps
                j1 = idx_ring(m);
                % Índice del vecino a 'step' posiciones adelante (circular)
                j2 = idx_ring(mod(m-1 + step, N_i) + 1);
                % b_j1 + b_j2 <= 1
                Arow = zeros(1, dim);
                Arow(2*N_total + j1) = 1;
                Arow(2*N_total + j2) = 1;
                A6 = [A6; Arow];
                b6 = [b6; 1];
            end
        end
    end
end

% 6) Restricción de separación mínima entre anillos
thresholdR = 0.4 * lmbda_eff;
A_ring = [];
b_ring = [];
% Los z_k están en las columnas 3*N_total+1 … 3*N_total+numRings
baseZ = 3*N_total;
for i = 1:numRings-1
    % separación radial
    if (ringRadii(i+1) - ringRadii(i)) <= thresholdR
        row = zeros(1, dim);
        % coeficientes para z_i y z_{i+1}
        row(baseZ + i)     = 1;
        row(baseZ + i + 1) = 1;
        A_ring = [A_ring; row];
        b_ring = [b_ring; 1];
    end
end

% para evitar optimizar todo, solo lo q no b=0
% Prealocar
lb_w = -1; ub_w = 2;       % para w_re y w_im
A7 = zeros(5*N_total, dim);
b7 = zeros(5*N_total,1);

M_w    = ub_w;           % 50
M_phase = Nlevels-1;     % 63

for j = 1:N_total
  row0 = (j-1)*5;

  % 1)  w_re_j ≤  M_w * b_j
  A7(row0+1, j            ) =  1;     % coef en w_re_j
  A7(row0+1, 2*N_total+j  ) = -M_w;   % –M_w·b_j

  % 2) -w_re_j ≤  M_w * b_j
  A7(row0+2, j            ) = -1;
  A7(row0+2, 2*N_total+j  ) = -M_w;

  % 3)  w_im_j ≤  M_w * b_j
  A7(row0+3, N_total+j       ) =  1; 
  A7(row0+3, 2*N_total+j     ) = -M_w;

  % 4) -w_im_j ≤  M_w * b_j
  A7(row0+4, N_total+j       ) = -1; 
  A7(row0+4, 2*N_total+j     ) = -M_w;

  % 5)  n_j ≤ M_phase * b_j
  %    columna de fases empieza en 3*N_total+numRings
  col_n = 3*N_total + numRings + j;
  A7(row0+5, col_n        ) =  1;
  A7(row0+5, 2*N_total+j  ) = -M_phase;
end

external_idx = find(abs(azimuth) > intermediateAngle);
B = 10^(-5/20)/sqrt(2);             % umbral lineal para -5 dB
S = numel(external_idx);    % número de ángulos externos

Aext = zeros(S, N_total);
for i = 1:S
    phi = azimuth(external_idx(i));
    sv  = exp(1i*2*pi*(elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi))/lmbda_eff);
    Aext(i,:) = sv.';       % 1×N_total
end

% Partir en real e imag
Are =  real(Aext);  % S×N
Aim =  imag(Aext);  % S×N

A_pat = zeros(4*S, dim);
b_pat =  B   * ones(4*S,1);  

% 1)  Re(w^H a)  ≤ B
%     sum_j [ Are(i,j)*w_re(j) - Aim(i,j)*w_im(j) ] ≤ B
A_pat(1:S, 1:N_total)       =  Are;
A_pat(1:S, N_total+1:2*N_total) = -Aim;

% 2) -Re(w^H a) ≤ B
A_pat(S+1:2*S, 1:N_total)       = -Are;
A_pat(S+1:2*S, N_total+1:2*N_total) =  Aim;

% 3)  Im(w^H a) ≤ B
%     sum_j [ Aim(i,j)*w_re(j) + Are(i,j)*w_im(j) ] ≤ B
A_pat(2*S+1:3*S, 1:N_total)     =  Aim;
A_pat(2*S+1:3*S, N_total+1:2*N_total) =  Are;

% 4) -Im(w^H a) ≤ B
A_pat(3*S+1:4*S, 1:N_total)     = -Aim;
A_pat(3*S+1:4*S, N_total+1:2*N_total) = -Are;

% Consolidar todas las restricciones lineales
Aineq = [ A1; A2; A3; A4; A5; A6; A_ring; A7; A_pat ];
bineq = [ b1; b2; b3; b4; b5; b6; b_ring; b7; b_pat ];

%% GA
% Llamada al GA con todas las restricciones
targetFun = @(x) macroArrayCostGAPRUEBA( ...
    x(1:3*N_total), ...       % Re(w), Im(w), b
    x(3*N_total+1:3*N_total+numRings), ... % z (solo restricciones)
    x(3*N_total+numRings+1:end), ...     % n (índices de fase)
    elementPos, ringIndex, azimuth, Beam_d, lmbda_eff, intermediateAngle);

options = optimoptions('ga', ...
    'Display', 'iter', ...
    'FunctionTolerance', 1e-8, ...
    'ConstraintTolerance', 0.3, ...
    'CrossoverFraction', 0.5, ...
    'PopulationSize', 400, ...
    'MaxGenerations', 1000, ...
    'EliteCount', 25, ...
    'MutationFcn', @mutationadaptfeasible, ...
    'SelectionFcn', @selectiontournament, ...
    'PlotFcn', @gaplotbestf);

[x_opt, fval, exitflag, output] = ga(targetFun, dim, Aineq, bineq, [], [], lb, ub, [], intcon, options);

%% Extracción y aplicación de resultados
w_re_opt = x_opt(1:N_total).';
w_im_opt = x_opt(N_total+1:2*N_total).';
w_raw_opt = w_re_opt + 1i*w_im_opt;
b_opt     = x_opt(2*N_total+1:3*N_total).';

%% FASES
% Extracción de fases discretas por elemento:
n_opt = x_opt(3*N_total+numRings+1 : 3*N_total+numRings+N_total);
% Convertir índices en grados:
phase_elem = 0 + n_opt * phaseStep;  % vector 1×N_total

% Sólo los elementos activos
active_idx   = find(b_opt >= 0.5);           % índices de elementos encendidos
phase_active = phase_elem(active_idx);       % sus fases

% Mostrar en consola
fprintf('Fase por elemento activo (índice, fase °):\n');
disp([active_idx(:), phase_active(:)]);

% Dibujo de tallo sólo para los activos
figure;
stem(active_idx, phase_active, 'filled');
xlabel('Índice de elemento');
ylabel('Fase (°)');
title('Phases of elements (6 bits)');
grid on;
%%
% Aplicar fase por elemento a los pesos finales:
w_opt = b_opt .* w_raw_opt;   % 1×N_total
for idx = 1:N_total
    w_opt(idx) = w_opt(idx) * exp(1i * deg2rad(phase_elem(idx)));
end
w_opt = w_opt(:).';  % fila
w_opt = w_opt(:).';  % asegurar fila

%% Síntesis del patrón utilizando la configuración optimizada
stvmat_opt = zeros(N_total, length(azimuth));
for i = 1:length(azimuth)
    phi = azimuth(i);
    stvmat_opt(:, i) = exp(1i * 2*pi * (elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi)) / lmbda_eff);
end
pattern_synth_opt = abs(w_opt * stvmat_opt);  % 1×length(azimuth)

figure;
plot(azimuth, mag2db(Beam_d), 'b', 'LineWidth',2);
hold on;
plot(azimuth, mag2db(pattern_synth_opt), 'r--', 'LineWidth',2);
legend('Deseado','Sintetizado','Location','Best');
xlabel('\theta (°)'); ylabel('Gain (dB)');
title('Optimized pattern');
grid on; ylim([-15 20]);

%% Visualización de la geometría del Macro Array (encendidos vs apagados)
figure;
active_idx   = find(b_opt >= 0.5);
inactive_idx = find(b_opt < 0.5);
scatter(elementPos(1,active_idx), elementPos(2,active_idx), 60, 'filled', 'MarkerFaceColor', 'b'); hold on;
scatter(elementPos(1,inactive_idx), elementPos(2,inactive_idx), 60, 'o', 'MarkerEdgeColor', 'r', 'LineWidth',1.5);

thetaC = linspace(0,2*pi,200);
for i = 1:numRings
    plot(ringRadii(i)*cos(thetaC), ringRadii(i)*sin(thetaC), 'k--');
end
axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)');
title('Macro Array (ON vs OFF)');
legend('ON','OFF');

pat_dB = mag2db(pattern_synth_opt);

fileID = fopen('resultado_cuant.txt','w');
fprintf(fileID,'%% theta_deg   G_dBi\n');
for k = 1:numel(azimuth)
    fprintf(fileID,'%8.3f   %8.4f\n', azimuth(k), pat_dB(k));
end
fclose(fileID);





%%
%%
%%
%%
% RESTRICCIONES 1/2/3/4/5/7

% % Script optimización de Macro Array sin desfases por anillo
% clear; clc; close all;
% 
% %%
% f0 = 1.3e9;
% lmbda = physconst('lightspeed')/f0;
% fprintf('λ = %.4f m\n', lmbda);
% 
% h = 1000; % altura en Km
% elevacion = deg2rad(20); % en grados
% r = 6378; % radio tierra en Km
% epsilon = 4.6;
% lmbda_eff = lmbda/sqrt(epsilon);
% fprintf('λ_eff = %.4f m\n', lmbda_eff);
% 
% slant_range = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));
% fprintf('slant range = %.4f m\n', slant_range);
% 
% %%
% azimuth = -90:90;
% 
% % Niveles en dB para el patrón:
% peak_dB         = 15;  
% valley_dB       = peak_dB - 20*log10(slant_range/h);  
% intermediate_dB = -5; 
% outer_dB        = -5;  
% 
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
% % tus ángulos “mágicos” (siguen intactos)
% % peakAngle         = 54;
% peakAngle = atan2(R_E*sin(phi_max), (R_E+h) - R_E*cos(phi_max)) * 180/pi;
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
% azimuth = -90:90;
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
% frameLen = 9;      % longitud de la ventana (impar, p.ej. 7, 9, 11)
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
% 
% numRings = 11;  % Ahora 7 anillos
% % Radios: primer anillo en 0.5·λ_eff, separación 0.25·λ_eff
% ringRadii = (0.5 + (0:numRings-1)*0.25) * lmbda_eff;
% % Número de elementos por anillo: [4,8,16,16,16,16,32]
% % ringElements = [4, 8, repmat(16,1,numRings-3), 32];
% ringElements = [4, 0, 8, 0, 8, 0, 8, 0, 16, 0, 16];
% N_total = sum(ringElements);
% fprintf('Total de elementos: %d\n', N_total);
% 
% 
% % Cálculo de posiciones: anillos concéntricos, elementos radiales
% elementPos = [];
% ringIndex = [];
% for i = 1:numRings
%     N_i = ringElements(i);
%     angles = (0:N_i-1) * (360 / N_i);  % grados uniformes
%     x_ring = ringRadii(i) * cosd(angles);
%     y_ring = ringRadii(i) * sind(angles);
%     elementPos = [elementPos, [x_ring; y_ring]];
%     ringIndex = [ringIndex, repmat(i, 1, N_i)];
% end
% 
% %% Máscara binaria inicial (todos encendidos)
% b0 = zeros(N_total,1);
% 
% %% Optimización con GA (pesos complejos, máscara on/off, activación de anillos y fases cuantizadas por elemento)
% % Restricción elementos activos: mínimo y máximo
% Mmin_el = 13;
% Mmax_el = 25;
% % Restricción anillos activos: mínimo y máximo
% Mmin_an = 3;
% Mmax_an = 4;
% 
% % Variables de fase discretas por elemento (6 bits → 64 niveles)
% % Número de variables de fase: N_total
% 
% % Dimensión del vector de optimización:
% % [Re(w) (N_total); Im(w) (N_total); b (N_total); z (numRings); n (N_total)]
% dim = 3*N_total + numRings;
% 
% % Vector inicial x0:
% x0 = [ ...
%     ones(N_total,1);      % Re(w)
%     zeros(N_total,1);     % Im(w)
%     zeros(N_total,1);      % b (todos encendidos inicial)
%     ones(numRings,1);     % z (todos anillos activos inicial)
% ];
% 
% % Cotas:
% lb = [ ...
%     -ones(N_total,1);     % Re(w) >= -1
%     -ones(N_total,1);   % Im(w) >= -1
%     zeros(N_total,1);     % b >= 0
%     zeros(numRings,1);    % z >= 0
% ];
% ub = [ ...
%     2*ones(N_total,1);    % Re(w) <= 2
%     2*ones(N_total,1);    % Im(w) <= 2
%     ones(N_total,1);      % b <= 1
%     ones(numRings,1);     % z <= 1
% ];
% 
% % Índices de variables enteras:
% intcon = [ ...
%     (2*N_total+1):(3*N_total), ...           % b
%     (3*N_total+1):(3*N_total+numRings), ...   % z
% ];
% 
% % Construir Aineq y bineq
% % 1) sum(b) >= Mmin_el  -> -sum(b) <= -Mmin_el
% %                Re(Im)   b          z              fases
% A1 = [zeros(1,2*N_total), -ones(1,N_total), zeros(1,numRings)];
% b1 = -Mmin_el;
% 
% % 2) sum(b) <= Mmax_el
% A2 = [zeros(1,2*N_total),  ones(1,N_total), zeros(1,numRings)];
% b2 =  Mmax_el;
% 
% % 3) sum(z) >= Mmin_an -> -sum(z) <= -Mmin_an
% A3 = [zeros(1,2*N_total), zeros(1,N_total), -ones(1,numRings)];
% b3 = -Mmin_an;
% 
% % 4) sum(z) <= Mmax_an
% A4 = [zeros(1,2*N_total), zeros(1,N_total),  ones(1,numRings)];
% b4 =  Mmax_an;
% 
% % 5) enlace de activación: b_j <= z_ring(j)
% A5 = zeros(N_total, dim);
% b5 = zeros(N_total,1);
% for j = 1:N_total
%     ring_j = ringIndex(j);
%     % b_j - z_k <= 0
%     A5(j, 2*N_total + j)      = 1;    % b_j
%     A5(j, 3*N_total + ring_j) = -1;   % -z_ring(j)
% end
% 
% % 6) Restricción de separación mínima entre elementos activos en cada anillo
% threshold = 0.25 * lmbda_eff;  % ajuste de separación deseada (p.ej. 0.4·λ_eff)
% A6 = [];
% b6 = [];
% for ring_i = 1:numRings
%     N_i = ringElements(ring_i);
%     radius_i = ringRadii(ring_i);
%     % Longitud de arco entre vecinos (un paso angular)
%     arcStep = 2*pi * radius_i / N_i;
%     % Número de pasos contiguos que quedan dentro del threshold
%     maxSteps = floor(threshold / arcStep);
%     if maxSteps >= 1
%         idx_ring = find(ringIndex == ring_i);
%         for m = 1:N_i
%             for step = 1:maxSteps
%                 j1 = idx_ring(m);
%                 % Índice del vecino a 'step' posiciones adelante (circular)
%                 j2 = idx_ring(mod(m-1 + step, N_i) + 1);
%                 % b_j1 + b_j2 <= 1
%                 Arow = zeros(1, dim);
%                 Arow(2*N_total + j1) = 1;
%                 Arow(2*N_total + j2) = 1;
%                 A6 = [A6; Arow];
%                 b6 = [b6; 1];
%             end
%         end
%     end
% end
% 
% % 6) Restricción de separación mínima entre anillos
% thresholdR = 0.4 * lmbda_eff;
% A_ring = [];
% b_ring = [];
% % Los z_k están en las columnas 3*N_total+1 … 3*N_total+numRings
% baseZ = 3*N_total;
% for i = 1:numRings-1
%     % separación radial
%     if (ringRadii(i+1) - ringRadii(i)) <= thresholdR
%         row = zeros(1, dim);
%         % coeficientes para z_i y z_{i+1}
%         row(baseZ + i)     = 1;
%         row(baseZ + i + 1) = 1;
%         A_ring = [A_ring; row];
%         b_ring = [b_ring; 1];
%     end
% end
% 
% % para evitar optimizar todo, solo lo q no b=0
% % Prealocar
% lb_w = -1; ub_w = 2;       % para w_re y w_im
% A7 = zeros(4*N_total, dim);
% b7 = zeros(4*N_total,1);
% 
% M_w    = ub_w;           % 50
% 
% for j = 1:N_total
%   row0 = (j-1)*4;
% 
%   % 1)  w_re_j ≤  M_w * b_j
%   A7(row0+1, j            ) =  1;     % coef en w_re_j
%   A7(row0+1, 2*N_total+j  ) = -M_w;   % –M_w·b_j
% 
%   % 2) -w_re_j ≤  M_w * b_j
%   A7(row0+2, j            ) = -1;
%   A7(row0+2, 2*N_total+j  ) = -M_w;
% 
%   % 3)  w_im_j ≤  M_w * b_j
%   A7(row0+3, N_total+j       ) =  1; 
%   A7(row0+3, 2*N_total+j     ) = -M_w;
% 
%   % 4) -w_im_j ≤  M_w * b_j
%   A7(row0+4, N_total+j       ) = -1; 
%   A7(row0+4, 2*N_total+j     ) = -M_w;
% end
% 
% external_idx = find(abs(azimuth) > intermediateAngle);
% B = 10^(-5/20)/sqrt(2);             % umbral lineal para -5 dB
% S = numel(external_idx);    % número de ángulos externos
% 
% Aext = zeros(S, N_total);
% for i = 1:S
%     phi = azimuth(external_idx(i));
%     sv  = exp(1i*2*pi*(elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi))/lmbda_eff);
%     Aext(i,:) = sv.';       % 1×N_total
% end
% 
% % Partir en real e imag
% Are =  real(Aext);  % S×N
% Aim =  imag(Aext);  % S×N
% 
% A_pat = zeros(4*S, dim);
% b_pat =  B   * ones(4*S,1);  
% 
% % 1)  Re(w^H a)  ≤ B
% %     sum_j [ Are(i,j)*w_re(j) - Aim(i,j)*w_im(j) ] ≤ B
% A_pat(1:S, 1:N_total)       =  Are;
% A_pat(1:S, N_total+1:2*N_total) = -Aim;
% 
% % 2) -Re(w^H a) ≤ B
% A_pat(S+1:2*S, 1:N_total)       = -Are;
% A_pat(S+1:2*S, N_total+1:2*N_total) =  Aim;
% 
% % 3)  Im(w^H a) ≤ B
% %     sum_j [ Aim(i,j)*w_re(j) + Are(i,j)*w_im(j) ] ≤ B
% A_pat(2*S+1:3*S, 1:N_total)     =  Aim;
% A_pat(2*S+1:3*S, N_total+1:2*N_total) =  Are;
% 
% % 4) -Im(w^H a) ≤ B
% A_pat(3*S+1:4*S, 1:N_total)     = -Aim;
% A_pat(3*S+1:4*S, N_total+1:2*N_total) = -Are;
% 
% % Consolidar todas las restricciones lineales
% Aineq = [ A1; A2; A3; A4; A5; A6; A_ring; A7; A_pat ];
% bineq = [ b1; b2; b3; b4; b5; b6; b_ring; b7; b_pat ];
% 
% %% GA
% % Llamada al GA con todas las restricciones
% targetFun = @(x) macroArrayCostGAPRUEBA( ...
%     x(1:3*N_total), ...       % Re(w), Im(w), b
%     x(3*N_total+1:3*N_total+numRings), ... % z (solo restricciones)
%     elementPos, ringIndex, azimuth, Beam_d, lmbda_eff, intermediateAngle);
% 
% % options = optimoptions('ga', ...
% %     'Display','iter', ...
% %     'UseParallel',true,...
% %     'UseVectorized',false,...
% %     'FunctionTolerance',1e-12, ...
% %     'ConstraintTolerance',1e-1, ...
% %     'CrossoverFraction',0.7, ...
% %     'PopulationSize',500, ...
% %     'MaxGenerations',600, ...
% %     'PlotFcn',@gaplotbestf);
% % % 'OutputFcn',@myOutputFcn,...
% 
% options = optimoptions('ga', ...
%     'Display', 'iter', ...
%     'FunctionTolerance', 1e-8, ...
%     'ConstraintTolerance', 0.3, ...
%     'CrossoverFraction', 0.5, ...
%     'PopulationSize', 400, ...
%     'MaxGenerations', 1000, ...
%     'EliteCount', 25, ...
%     'MutationFcn', @mutationadaptfeasible, ...
%     'SelectionFcn', @selectiontournament, ...
%     'PlotFcn', @gaplotbestf);
% % 'UseParallel', true, ...
% % 'HybridFcn', @fmincon, ...
% 
% [x_opt, fval, exitflag, output] = ga(targetFun, dim, Aineq, bineq, [], [], lb, ub, [], intcon, options);
% % [x_opt, fval, exitflag, output] = ga(targetFun, dim, Aineq, bineq, [], [], lb, ub, [], options);
% 
% %% Extracción y aplicación de resultados
% w_re_opt = x_opt(1:N_total).';
% w_im_opt = x_opt(N_total+1:2*N_total).';
% w_raw_opt = w_re_opt + 1i*w_im_opt;
% b_opt     = x_opt(2*N_total+1:3*N_total).';
% 
% 
% %%
% % Aplicar fase por elemento a los pesos finales:
% w_opt = b_opt .* w_raw_opt;   % 1×N_total
% w_opt = w_opt(:).';  % fila
% w_opt = w_opt(:).';  % asegurar fila
% 
% %% Síntesis del patrón utilizando la configuración optimizada
% stvmat_opt = zeros(N_total, length(azimuth));
% for i = 1:length(azimuth)
%     phi = azimuth(i);
%     stvmat_opt(:, i) = exp(1i * 2*pi * (elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi)) / lmbda_eff);
% end
% pattern_synth_opt = abs(w_opt * stvmat_opt);  % 1×length(azimuth)
% 
% figure;
% plot(azimuth, mag2db(Beam_d), 'b', 'LineWidth',2);
% hold on;
% plot(azimuth, mag2db(pattern_synth_opt), 'r--', 'LineWidth',2);
% legend('Desired','Synthesized','Location','Best');
% xlabel('\theta (°)'); ylabel('Gain (dB)');
% title('Optimized pattern');
% grid on; ylim([-15 20]);
% 
% %% Visualización de la geometría del Macro Array (encendidos vs apagados)
% figure;
% active_idx   = find(b_opt >= 0.5);
% inactive_idx = find(b_opt < 0.5);
% scatter(elementPos(1,active_idx), elementPos(2,active_idx), 60, 'filled', 'MarkerFaceColor', 'b'); hold on;
% scatter(elementPos(1,inactive_idx), elementPos(2,inactive_idx), 60, 'o', 'MarkerEdgeColor', 'r', 'LineWidth',1.5);
% 
% thetaC = linspace(0,2*pi,200);
% for i = 1:numRings
%     plot(ringRadii(i)*cos(thetaC), ringRadii(i)*sin(thetaC), 'k--');
% end
% axis equal; grid on;
% xlabel('X (m)'); ylabel('Y (m)');
% title('Macro Array (ON vs OFF)');
% legend('ON','OFF');
% 
% pat_dB = mag2db(pattern_synth_opt);
% 
% fileID = fopen('resultado_2.txt','w');
% fprintf(fileID,'%% theta_deg   G_dBi\n');
% for k = 1:numel(azimuth)
%     fprintf(fileID,'%8.3f   %8.4f\n', azimuth(k), pat_dB(k));
% end
% fclose(fileID);