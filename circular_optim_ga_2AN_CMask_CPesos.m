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
N1 = 8;
N2 = 8;
N_total = N1 + N2;


objfun = @(x) CostGA2ANCmaskCPesos(x, N1, N2, azimuth, Beam_d, lmbda, intermediateAngle);

%%
dim = 3 + 2*N_total;

x0 = zeros(dim,1);
x0(1) = 0.5*lmbda;
x0(2) = 1.0*lmbda;
x0(3) = 0;
x0(4 : 3+N_total) = 1;
x0(4+N_total : 3+2*N_total) = 0;

% Límites:
lb = zeros(dim,1);
ub = zeros(dim,1);

% Para r1:
lb(1) = 0.45*lmbda;
ub(1) = 0.55*lmbda;
% Para r2:
lb(2) = 0.9*lmbda;
ub(2) = 1.1*lmbda;
% Para desfaseRel:
lb(3) = -60;
ub(3) = 60;
% Para los pesos:
lb(4 : 3+N_total) = -1;
ub(4 : 3+N_total) = 2;
lb(4+N_total : 3+2*N_total) = -1;
ub(4+N_total : 3+2*N_total) = 2;


%%

options = optimoptions('ga', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'ConstraintTolerance',1e-3, ...
    'CrossoverFraction',0.7, ...
    'PopulationSize',300, ...
    'MaxGenerations',300, ...
    PlotFcn='gaplotbestf');

[x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [], options);


%% 
r1_opt = x_opt(1);
r2_opt = x_opt(2);
desfaseRel_opt = x_opt(3);

w_re_opt = x_opt(4 : 3+N_total).';
w_im_opt = x_opt(4+N_total : 3+2*N_total).';
w_opt = w_re_opt + 1i*w_im_opt;

anglesRing1_opt = (0:N1-1)*(360/N1);
anglesRing2_opt = (0:N2-1)*(360/N2) + desfaseRel_opt;

posRing1_opt = [r1_opt*cosd(anglesRing1_opt); r1_opt*sind(anglesRing1_opt)];
posRing2_opt = [r2_opt*cosd(anglesRing2_opt); r2_opt*sind(anglesRing2_opt)];
elementPos_opt = [posRing1_opt, posRing2_opt];

stvmat_opt = zeros(N_total, length(azimuth));
for i = 1:length(azimuth)
    phi = azimuth(i);
    stvmat_opt(:, i) = exp(1i * 2*pi * ( elementPos_opt(1,:)'*cosd(phi) + ...
                                           elementPos_opt(2,:)'*sind(phi) ) / lmbda );
end

w_opt = w_opt(:);
pattern_synth_opt = abs(w_opt' * stvmat_opt);

fprintf('phase shift = %.4f º\n', desfaseRel_opt);
fprintf('Radius Ring 1 = %.4f m\n', r1_opt);
fprintf('Radius Ring 2 = %.4f m\n', r2_opt);

%%
figure;
plot(azimuth, mag2db([Beam_d; pattern_synth_opt])', 'LineWidth',2);
legend('Desired','Synthesized','Location','Best');
xlabel('\theta (°)'); ylabel('Gain (dB)');
title('Optimized Pattern with GA and Linear Constraints');
grid on; ylim([-10 20]);

figure;
scatter(elementPos_opt(1,:), elementPos_opt(2,:), 60, 'filled', 'b');
hold on;
thetaC = linspace(0,2*pi,200);
plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--');
plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--');
axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)');
title('Optimized Array Positions (Circular, 16 elements)');
