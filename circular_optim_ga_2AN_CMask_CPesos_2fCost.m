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


%%
N1 = 8;% Elementos en anillo 1
N2 = 8;% Elementos en anillo 2
N_total = N1 + N2;

%%
objfun = @(x) CostGA2ANCmaskCPesos2FCost(x, N1, N2, azimuth, Beam_d, lmbda_eff, intermediateAngle);

%%
dim = 3 + 2*N_total;  % 3 variables geométricas + 2*16 = 35

x0 = zeros(dim,1);
x0(1) = 0.5*lmbda_eff; % r1 
x0(2) = 1.0*lmbda_eff; % r2
x0(3) = 0; % desfaseRel inicial
x0(4 : 3+N_total) = 1; % parte real a 1
x0(4+N_total : 3+2*N_total) = 0;  % parte imaginaria a 0

% Definir límites:
lb = zeros(dim,1);
ub = zeros(dim,1);

% Límites para r1
lb(1) = 0.45*lmbda_eff;
ub(1) = 0.55*lmbda_eff;
% Límites para r2
lb(2) = 0.9*lmbda_eff;
ub(2) = 1.1*lmbda_eff;
% Límites para el desfase relativo (en grados)
lb(3) = -60;
ub(3) = 60;
% Límites para la parte real de los pesos
lb(4 : 3+N_total) = -1;
ub(4 : 3+N_total) = 2;
% Límites para la parte imaginaria de los pesos
lb(4+N_total : 3+2*N_total) = -1;
ub(4+N_total : 3+2*N_total) = 2;

%%
options = optimoptions('ga', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8,...
    'ConstraintTolerance',1e-3,...
    'CrossoverFraction',0.7,...
    'PopulationSize',300, ...
    'MaxGenerations',300,...
    PlotFcn='gaplotbestf');

[x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [], options);

r1_opt         = x_opt(1);
r2_opt         = x_opt(2);
desfaseRel_opt = x_opt(3);

w_re_opt = x_opt(4 : 3+N_total).';
w_im_opt = x_opt(4+N_total : 3+2*N_total).';
w_opt = w_re_opt + 1i*w_im_opt;


%% 
anglesRing1_opt = (0:N1-1)*(360/N1);  
anglesRing2_opt = (0:N2-1)*(360/N2) + desfaseRel_opt;

posRing1_opt = [r1_opt*cosd(anglesRing1_opt); r1_opt*sind(anglesRing1_opt)];
posRing2_opt = [r2_opt*cosd(anglesRing2_opt); r2_opt*sind(anglesRing2_opt)];
elementPos_opt = [posRing1_opt, posRing2_opt];

% matriz de apuntamiento
stvmat_opt = zeros(N_total, length(azimuth));
for i = 1:length(azimuth)
    phi = azimuth(i);
    stvmat_opt(:, i) = exp(1i * 2*pi * ( elementPos_opt(1,:)'*cosd(phi) + elementPos_opt(2,:)'*sind(phi) ) / lmbda_eff );
end

pattern_synth_opt = abs(w_opt'*stvmat_opt);

fprintf('phase shift = %.4f º\n', desfaseRel_opt);
fprintf('Radius Ring 1 = %.4f m\n', r1_opt);
fprintf('Radius Ring 2 = %.4f m\n', r2_opt);

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
axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)');
title('Optimized Array Positions (Circular, 16 elements)');



%% PARTE 2: SIMULACIÓN DE LA HUELLA (PFD) EN TIERRA
%% Parámetros
txPower = 10; % Watts

slant_range = 2121;
off_nadir_lim = 54.33;

satAlt = 1000e3;
satLat  = 40.4; % Latitud (aprox. Madrid)
satLon  = -3.7; % Longitud (aprox. Madrid)

minElev = 20; % Elevación mínima de diseño (°)
maxOffNadir = 90 - minElev; % Off-nadir máximo permitido (°)

azimuth = -90:90;                     
pattern_interp = @(angle) interp1(azimuth, pattern_synth_opt, angle, 'linear', 'extrap');

%% Estaciones base
stations = { 'Cadiz', 36.5, -6.3;
             'Tarifa', 36, -5.6;
             'Barcelona', 41.4, 2.2;
             'La Coruna', 43.3, -8.4;
             'Ponferrada', 42.9, -6.7;
             'Madrid', 40.4, -3.7;
             'Salamanca', 41, -5.96;
             'Albacete', 39, -1.85;
             'Paris', 48.8, 2.35;
             'Londres', 51.7, 0};
% NO ES PONFERRADA COMO TAL PERO ESTA CERCA (ES JUSTO EL LIMITE ENTRE
% ASTURIAS Y LEON)

numStations = size(stations,1);
PFD_stations = zeros(numStations,1);

fprintf('--- Estaciones Base ---\n');
for i = 1:numStations
    gsLat = stations{i,2};
    gsLon = stations{i,3};

    dSurf = haversine(satLat, satLon, gsLat, gsLon);
    R = sqrt(dSurf^2 + satAlt^2);
    elev = atand(satAlt / dSurf);
    offNadir = min(90 - elev, maxOffNadir);

    G = pattern_interp(offNadir);
    EIRP = txPower * G;
    PFD_stations(i) = EIRP / (4 * pi * R^2);

    PFD_dBW = 10 * log10(PFD_stations(i));

    fprintf('%s: off-nadir = %.2f° , PFD = %.2e W/m², = %.2f dBW, distancia = %.2f km\n', ...
        stations{i,1}, offNadir, PFD_stations(i), PFD_dBW, R/1000);
end

%% Calcular la PFD en una malla fina (región de interés)
dLat = 0.1; % Paso en latitud
dLon = 0.1; % Paso en longitud
latVec = 25:dLat:57;  
lonVec = -25:dLon:17;
[lonGrid, latGrid] = meshgrid(lonVec, latVec);

PFD_map = zeros(size(latGrid));
for idx = 1:numel(latGrid)
    PFD_map(idx) = computePFD(latGrid(idx), lonGrid(idx), ...
                              satLat, satLon, satAlt, ...
                              txPower, pattern_interp, maxOffNadir);
end

PFD_map_dBW = 10 * log10(PFD_map);

% puntos dentro del área de cobertura
groundRange = slant_range * sind(off_nadir_lim);  % km

distances = arrayfun(@(lat, lon) haversine(satLat, satLon, lat, lon), latGrid, lonGrid);
mask = distances <= groundRange * 1000;  % Solo puntos dentro del ground range

latInside = latGrid(mask);
lonInside = lonGrid(mask);
PFD_map_dBW_inside = PFD_map_dBW(mask);

%% Visualización
figure
ax = geoaxes;
geobasemap(ax, 'satellite');
geolimits(ax, [25 60], [-25 20]);
hold(ax, 'on');

R_earth_km = 6378;
angularRadius = rad2deg(groundRange / R_earth_km);  % grados
[latCircle, lonCircle] = scircle1(satLat, satLon, angularRadius, R_earth_km, 'degrees');
geoplot(ax, latCircle, lonCircle, 'r', 'LineWidth', 2);

% "Mapa de calor" de la PFD en dBW
geoscatter(ax, latInside, lonInside, 100, PFD_map_dBW_inside, ...
           'filled', 's', 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.2);
colormap(ax, jet);
cb = colorbar(ax);
cb.Label.String = 'PFD (dBW)';

%estaciones base en el mapa
for i = 1:numStations
    gsLat = stations{i,2};
    gsLon = stations{i,3};
    geoplot(ax, gsLat, gsLon, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'red');
    text(gsLat + 0.3, gsLon, stations{i,1}, 'Color', 'yellow', 'FontWeight', 'bold');
end

title(ax, 'PFD en Tierra (dBW), Estaciones Base y Límite de Cobertura');
hold(ax, 'off');

%% Funciones locales

function val = computePFD(lat, lon, satLat, satLon, satAlt, txPower, pattern_interp, maxOffNadir)
    d = haversine(satLat, satLon, lat, lon);
    R = sqrt(d^2 + satAlt^2);
    elev = atand(satAlt / d);
    offNadir = min(90 - elev, maxOffNadir);
    G = pattern_interp(offNadir);
    EIRP = txPower * G;
    val = EIRP / (4 * pi * R^2);
end

function d = haversine(lat1, lon1, lat2, lon2)
    R = 6378e3;
    dLat = deg2rad(lat2 - lat1);
    dLon = deg2rad(lon2 - lon1);
    a = sin(dLat/2)^2 + cosd(lat1) * cosd(lat2) * sin(dLon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    d = R * c;
end
