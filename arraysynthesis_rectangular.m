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
    % valle isoflux geométrico
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

%%
M = 10; 
N = 10;
N_total = M * N;

d = 2*lmbda; % Espaciado

% malla
[x_grid, y_grid] = meshgrid(0:N-1, 0:M-1);
% Posiciones de los elementos en metros
elementPos = [x_grid(:)'; y_grid(:)'] * d;  % 2 x 16

nvars  = 2*N_total;

% Cálculo del steering vector para cada ángulo de azimut
numAngles = length(azimuth);
stvmat = zeros(N_total, numAngles);
for i = 1:numAngles
    phi = azimuth(i);
    stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)' * cosd(phi) + elementPos(2,:)' * sind(phi) ) );
end

w_i_re = ones(N_total, 1);
w_i_im = zeros(N_total, 1);
x_ini = [w_i_re; w_i_im];

Beam_d_col = Beam_d(:);

% Límites de la parte real e imaginaria en [-10, 10]
lb = -10*ones(2*N_total, 1);
ub = 10*ones(2*N_total, 1);

% función objetivo
objfun = @(x) norm( ...
    abs( stvmat.' * ( ( x(1:N_total) + 1i*x(N_total+1:end) ) .') ) ...  
    - Beam_d_col );

options = optimoptions('ga', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'ConstraintTolerance',1e-3, ...
    'CrossoverFraction',0.7, ...
    'PopulationSize',100, ...
    'MaxGenerations',200, ...
    'PlotFcn', @gaplotbestf);

x_opt = ga(objfun, nvars, [], [], [], [], lb, ub, [], options);

w_opt = x_opt(1:N_total) + 1i*x_opt(N_total+1:end);

pattern_synth = abs(w_opt * stvmat);

figure;
plot(azimuth, mag2db([Beam_d; pattern_synth])','LineWidth',2);
legend('desired','synthesized','Location','Best');
xlabel('\theta (°)'); ylabel('Gain (dB)');
title('Isoflux pattern');
grid on; ylim([-10 20]);
