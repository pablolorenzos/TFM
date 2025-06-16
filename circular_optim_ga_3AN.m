%%
% En este código hay varios ejemplos
% SIN ELEMENTO CENTRAL
% CON ELEMENTO CENTRAL
% OPTIMIZACIÓN ELEMNETOS SIN ELEMENTO CENTRAL
% OPTIMIZACION ELEMENTON CON ELEMENTO CENTRAL
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


% peakAngle         = 54;
peakAngle = atan2(R_E*sin(phi_max), (R_E+h) - R_E*cos(phi_max)) * 180/pi;
fprintf('peak angle for pattern = %.3f º\n', peakAngle);
intermediateAngle = 70; 
outerAngle        = 100;  

% niveles
peak_dB         = 15;
intermediate_dB = peak_dB - 20;
outer_dB        = peak_dB - 20;

delta = 1;    % grados de “ancho” extra del pico

azimuth = -100:100;
Beam_dB = zeros(size(azimuth));

for k = 1:numel(azimuth)
  th = azimuth(k);  x = abs(th);
  if x <= (peakAngle-delta)
    % --- valle isoflux geométrico ---
    % interpolas la curva “exacta” y la elevas al nivel del pico
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
xlabel('\theta (°)'); ylabel('Gain (dB)');
title('Desired pattern');
grid on; ylim([-15 20]);

%%
N1 = 4;% Elementos en anillo 1
N2 = 6;% Elementos en anillo 2
N3 = 8;% Elementos en anillo 3
N_total = N1 + N2 + N3;

%%
objfun = @(x) CostGA3AN(x, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle);

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
lb(6+N_total : 3+5*N_total) = -50;
ub(6+N_total : 3+5*N_total) = 50;

%%


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
    'PopulationSize',500, ...
    'MaxGenerations',1750, ...
    'PlotFcn','gaplotbestf', ...
    'OutputFcn',ofcn);


[x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [], options);

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

% 1) Dibujamos los círculos de anillo SÓLO para referencia, SIN que salgan en la leyenda
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
        [~, cC, cE] = CostGA3AN( xBest, ...
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

%% PARTE 2: SIMULACIÓN DE LA HUELLA (PFD) EN TIERRA
%% Parámetros
txPower = 10; % Watts

% off_nadir_lim = 54.33;
off_nadir_lim = peakAngle;

satAlt = h;
satLat  = 40.4; % Latitud (aprox. Madrid)
satLon  = -3.7; % Longitud (aprox. Madrid)

% minElev = 20; % Elevación mínima de diseño (°)
minElev = rad2deg(elevacion);
maxOffNadir = 90 - minElev; % Off-nadir máximo permitido (°)
                   
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

% ------------------------------------- curva roja (pq no se me superpone)
R_earth_km = r;
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
% azimuth = -90:90;
% 
% % Niveles en dB para el patrón:
% peak_dB         = 15;  
% valley_dB       = peak_dB - 20*log10(slant_range/h);  
% intermediate_dB = -5; 
% outer_dB        = -5;  
% 
% % Ángulos en grados
% peakAngle         = 20;  
% intermediateAngle = 40;  
% outerAngle        = 80;  
% 
% % Parámetros
% R_E       = 6378;        % km, radio de la Tierra
% slant_ref = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));        % km, rango en los picos
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
% objfun = @(x) CostGA3AN(x, N0, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle);
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
% lb(4) = 0;
% ub(4) = 0;
% % Límites para el desfase relativo1 (en grados)
% lb(5) = 0;
% ub(5) = 0;
% % Límites para la parte real de los pesos
% lb(6 : 5+N_total) = -50;
% ub(6 : 5+N_total) = 50;
% % Límites para la parte imaginaria de los pesos
% lb(6+N_total : 3+5*N_total) = -50;
% ub(6+N_total : 3+5*N_total) = 50;
% 
% %%
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
%     'MaxGenerations',2, ...
%     'PlotFcn','gaplotbestf', ...
%     'OutputFcn',ofcn);
% 
% [x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [], options);
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
% %% Gráfico numerado con leyenda
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
% % 3) Ahora la leyenda la creamos sólo con los handles de los elementos
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
%         [~, cC, cE] = CostGA3AN( xBest, ...
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
% 
% 
% %% PARTE 2: SIMULACIÓN DE LA HUELLA (PFD) EN TIERRA
% %% Parámetros
% txPower = 10; % Watts
% 
% slant_range = 2121;
% off_nadir_lim = 54.33;
% 
% satAlt = 1000e3;
% satLat  = 40.4; % Latitud (aprox. Madrid)
% satLon  = -3.7; % Longitud (aprox. Madrid)
% 
% minElev = 20; % Elevación mínima de diseño (°)
% maxOffNadir = 90 - minElev; % Off-nadir máximo permitido (°)
% 
% azimuth = -90:90;                     
% pattern_interp = @(angle) interp1(azimuth, pattern_synth_opt, angle, 'linear', 'extrap');
% 
% %% Estaciones base
% stations = { 'Cadiz', 36.5, -6.3;
%              'Tarifa', 36, -5.6;
%              'Barcelona', 41.4, 2.2;
%              'La Coruna', 43.3, -8.4;
%              'Ponferrada', 42.9, -6.7;
%              'Madrid', 40.4, -3.7;
%              'Salamanca', 41, -5.96;
%              'Albacete', 39, -1.85;
%              'Paris', 48.8, 2.35;
%              'Londres', 51.7, 0};
% % NO ES PONFERRADA COMO TAL PERO ESTA CERCA (ES JUSTO EL LIMITE ENTRE
% % ASTURIAS Y LEON)
% 
% numStations = size(stations,1);
% PFD_stations = zeros(numStations,1);
% 
% fprintf('--- Estaciones Base ---\n');
% for i = 1:numStations
%     gsLat = stations{i,2};
%     gsLon = stations{i,3};
% 
%     dSurf = haversine(satLat, satLon, gsLat, gsLon);
%     R = sqrt(dSurf^2 + satAlt^2);
%     elev = atand(satAlt / dSurf);
%     offNadir = min(90 - elev, maxOffNadir);
% 
%     G = pattern_interp(offNadir);
%     EIRP = txPower * G;
%     PFD_stations(i) = EIRP / (4 * pi * R^2);
% 
%     PFD_dBW = 10 * log10(PFD_stations(i));
% 
%     fprintf('%s: off-nadir = %.2f° , PFD = %.2e W/m², = %.2f dBW, distancia = %.2f km\n', ...
%         stations{i,1}, offNadir, PFD_stations(i), PFD_dBW, R/1000);
% end
% 
% %% Calcular la PFD en una malla fina (región de interés)
% dLat = 0.1; % Paso en latitud
% dLon = 0.1; % Paso en longitud
% latVec = 25:dLat:57;  
% lonVec = -25:dLon:17;
% [lonGrid, latGrid] = meshgrid(lonVec, latVec);
% 
% PFD_map = zeros(size(latGrid));
% for idx = 1:numel(latGrid)
%     PFD_map(idx) = computePFD(latGrid(idx), lonGrid(idx), ...
%                               satLat, satLon, satAlt, ...
%                               txPower, pattern_interp, maxOffNadir);
% end
% 
% PFD_map_dBW = 10 * log10(PFD_map);
% 
% % puntos dentro del área de cobertura
% groundRange = slant_range * sind(off_nadir_lim);  % km
% 
% distances = arrayfun(@(lat, lon) haversine(satLat, satLon, lat, lon), latGrid, lonGrid);
% mask = distances <= groundRange * 1000;  % Solo puntos dentro del ground range
% 
% latInside = latGrid(mask);
% lonInside = lonGrid(mask);
% PFD_map_dBW_inside = PFD_map_dBW(mask);
% 
% %% Visualización
% figure
% ax = geoaxes;
% geobasemap(ax, 'satellite');
% geolimits(ax, [25 60], [-25 20]);
% hold(ax, 'on');
% 
% R_earth_km = 6378;
% angularRadius = rad2deg(groundRange / R_earth_km);  % grados
% [latCircle, lonCircle] = scircle1(satLat, satLon, angularRadius, R_earth_km, 'degrees');
% geoplot(ax, latCircle, lonCircle, 'r', 'LineWidth', 2);
% 
% % "Mapa de calor" de la PFD en dBW
% geoscatter(ax, latInside, lonInside, 100, PFD_map_dBW_inside, ...
%            'filled', 's', 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.2);
% colormap(ax, jet);
% cb = colorbar(ax);
% cb.Label.String = 'PFD (dBW)';
% 
% %estaciones base en el mapa
% for i = 1:numStations
%     gsLat = stations{i,2};
%     gsLon = stations{i,3};
%     geoplot(ax, gsLat, gsLon, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'red');
%     text(gsLat + 0.3, gsLon, stations{i,1}, 'Color', 'yellow', 'FontWeight', 'bold');
% end
% 
% title(ax, 'PFD en Tierra (dBW), Estaciones Base y Límite de Cobertura');
% hold(ax, 'off');
% 
% %% Funciones locales
% 
% function val = computePFD(lat, lon, satLat, satLon, satAlt, txPower, pattern_interp, maxOffNadir)
%     d = haversine(satLat, satLon, lat, lon);
%     R = sqrt(d^2 + satAlt^2);
%     elev = atand(satAlt / d);
%     offNadir = min(90 - elev, maxOffNadir);
%     G = pattern_interp(offNadir);
%     EIRP = txPower * G;
%     val = EIRP / (4 * pi * R^2);
% end
% 
% function d = haversine(lat1, lon1, lat2, lon2)
%     R = 6378e3;
%     dLat = deg2rad(lat2 - lat1);
%     dLon = deg2rad(lon2 - lon1);
%     a = sin(dLat/2)^2 + cosd(lat1) * cosd(lat2) * sin(dLon/2)^2;
%     c = 2 * atan2(sqrt(a), sqrt(1 - a));
%     d = R * c;
% end
% 


%%
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
% OPTIMIZAMOS NUMERO DE ELEMENTOS CON ELEMENTO CENTRAL

%%

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
% fprintf('λ = %.4f m\n', lmbda_eff);
% 
% slant_range = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));
% fprintf('slant range = %.4f m\n', slant_range);
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
% Beam_d = Beam_sg;
% 
% figure;
% plot(azimuth, mag2db(Beam_d), 'b', 'LineWidth',2);
% xlabel('\theta (°)'); ylabel('Ganancia (dB)');
% title('Patrón deseado (half-cosine)');
% grid on; ylim([-15 20]);
% 
% %%
% 
% N0min = 0;
% N0max = 1;
% N1min = 2;
% N1max = 4;
% N2min = 4;
% N2max = 8;
% N3min = 4;
% N3max = 10;
% 
% % Initial guess for the number of elements in each ring
% N0 = 0;
% N1 = 4; % Initial guess for ring 1
% N2 = 6; % Initial guess for ring 2
% N3 = 8; % Initial guess for ring 3
% 
% % Total number of elements (initial guess)
% N_total = N0 + N1 + N2 + N3;
% max_N_total = N0max + N1max + N2max + N3max;
% 
% 
% %%
% objfun = @(x) CostGA3AN(x, azimuth, Beam_d, lmbda_eff, intermediateAngle);
% 
% %%
% dim = 9 + 2*max_N_total;  % 5 variables geométricas + 2*16 = 35
% 
% x0 = zeros(dim,1);
% x0(1) = 0.5*lmbda_eff; % r1 
% x0(2) = 1.0*lmbda_eff; % r2
% x0(3) = 1.5*lmbda_eff; % r3
% x0(4) = 0; % desfaseRel inicial for ring 2
% x0(5) = 0; % desfaseRel inicial for ring 3
% x0(6) = N0;
% x0(7) = N1; % Initial guess for N1
% x0(8) = N2; % Initial guess for N2
% x0(9) = N3; % Initial guess for N3
% x0(10 : 9+max_N_total) = 1; % Real part of weights
% x0(10+max_N_total : 9+2*max_N_total) = 0; % Imaginary part of weights
% 
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
% % Bounds for N1, N2, N3
% lb(6) = N0min; 
% ub(6) = N0max;
% lb(7) = N1min; 
% ub(7) = N1max;
% lb(8) = N2min; 
% ub(8) = N2max;
% lb(9) = N3min; 
% ub(9) = N3max;
% % Bounds for real and imaginary parts of weights
% lb(10 : 9+max_N_total) = -50;
% ub(10 : 9+max_N_total) = 50;
% lb(10+max_N_total : 9+2*max_N_total) = -50;
% ub(10+max_N_total : 9+2*max_N_total) = 50;
% 
% % Specify integer constraints for N1, N2, N3
% IntCon = 6:9;

% %%
% 
% % Arrays donde se irán guardando ambos costes
% centralHistory = [];
% externalHistory = [];
% 
% 
% ofcn = @(options,state,flag) costMonitor(...
%     options, state, flag, ...
%     azimuth, Beam_d, ...
%     lmbda_eff, intermediateAngle);
% 
% options = optimoptions('ga', ...
%     'Display','iter', ...
%     'FunctionTolerance',1e-8, ...
%     'ConstraintTolerance',1e-3, ...
%     'CrossoverFraction',0.7, ...
%     'PopulationSize',300, ...
%     'MaxGenerations',800, ...
%     'PlotFcn','gaplotbestf', ...
%     'OutputFcn',ofcn);
% 
% [x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [],IntCon, options);
% 
% 
% figure; hold on;
% gens = 0:length(centralHistory)-1;
% plot(gens, centralHistory, '-b','LineWidth',2);
% plot(gens, externalHistory, '-r','LineWidth',2);
% legend('Central cost','External cost','Location','best');
% xlabel('Generation'); ylabel('Coste');
% title('Evolution Cost Central vs External');
% grid on;
% 
% 
% % Extract results
% r1_opt = x_opt(1);
% r2_opt = x_opt(2);
% r3_opt = x_opt(3);
% desfaseRel_opt = x_opt(4);
% desfaseRel1_opt = x_opt(5);
% N0_opt = round(x_opt(6));
% N1_opt = round(x_opt(7));
% N2_opt = round(x_opt(8));
% N3_opt = round(x_opt(9));
% 
% N_total_opt = N0_opt + N1_opt + N2_opt + N3_opt;
% 
% w_re_opt = x_opt(10 : 9+N_total_opt).';
% w_im_opt = x_opt(10+N_total_opt : 9+2*N_total_opt).';
% w_opt = w_re_opt + 1i*w_im_opt;
% 
% 
% %% 
% % centro:
% 
% if N0_opt > 0
%    pos0 = [0; 0];
% else
%    pos0 = [];
% end
% 
% anglesRing1_opt = (0:N1_opt-1)*(360/N1_opt);  
% anglesRing2_opt = (0:N2_opt-1)*(360/N2_opt) + desfaseRel_opt;
% anglesRing3_opt = (0:N3_opt-1)*(360/N3_opt) + desfaseRel1_opt;
% 
% posRing1_opt = [r1_opt*cosd(anglesRing1_opt); r1_opt*sind(anglesRing1_opt)];
% posRing2_opt = [r2_opt*cosd(anglesRing2_opt); r2_opt*sind(anglesRing2_opt)];
% posRing3_opt = [r3_opt*cosd(anglesRing3_opt); r3_opt*sind(anglesRing3_opt)];
% elementPos_opt = [pos0, posRing1_opt, posRing2_opt, posRing3_opt];
% 
% 
% % matriz de apuntamiento
% stvmat_opt = zeros(N_total_opt, length(azimuth));
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
% %% Gráfico con leyenda
% figure; hold on;
% 
% hC1 = plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% hC2 = plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% hC3 = plot(r3_opt*cos(thetaC), r3_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% 
% % 2) Para cada elemento creamos su propio scatter, con DisplayName = "n: módulo ∠ fase"
% hElem = gobjects(N_total_opt,1);
% for k = 1:N_total_opt
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
% % 3) Ahora la leyenda la creamos sólo con los handles de los elementos
% legend(hElem, 'Location','eastoutside', 'FontSize',8);
% 
% 
% 
% %%
% function [state,options,optchanged] = costMonitor( ...
%         options, state, flag, ...
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
%         [~, cC, cE] = CostGA3AN( xBest, ...
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
% 
% 
% 
% %% PARTE 2: SIMULACIÓN DE LA HUELLA (PFD) EN TIERRA
% %% Parámetros
% txPower = 10; % Watts
% 
% slant_range = 2121;
% off_nadir_lim = 54.33;
% 
% satAlt = 1000e3;
% satLat  = 40.4; % Latitud (aprox. Madrid)
% satLon  = -3.7; % Longitud (aprox. Madrid)
% 
% minElev = 20; % Elevación mínima de diseño (°)
% maxOffNadir = 90 - minElev; % Off-nadir máximo permitido (°)
% 
% azimuth = -90:90;                     
% pattern_interp = @(angle) interp1(azimuth, pattern_synth_opt, angle, 'linear', 'extrap');
% 
% %% Estaciones base
% stations = { 'Cadiz', 36.5, -6.3;
%              'Tarifa', 36, -5.6;
%              'Barcelona', 41.4, 2.2;
%              'La Coruna', 43.3, -8.4;
%              'Ponferrada', 42.9, -6.7;
%              'Madrid', 40.4, -3.7;
%              'Salamanca', 41, -5.96;
%              'Albacete', 39, -1.85;
%              'Paris', 48.8, 2.35;
%              'Londres', 51.7, 0};
% % NO ES PONFERRADA COMO TAL PERO ESTA CERCA (ES JUSTO EL LIMITE ENTRE
% % ASTURIAS Y LEON)
% 
% numStations = size(stations,1);
% PFD_stations = zeros(numStations,1);
% 
% fprintf('--- Estaciones Base ---\n');
% for i = 1:numStations
%     gsLat = stations{i,2};
%     gsLon = stations{i,3};
% 
%     dSurf = haversine(satLat, satLon, gsLat, gsLon);
%     R = sqrt(dSurf^2 + satAlt^2);
%     elev = atand(satAlt / dSurf);
%     offNadir = min(90 - elev, maxOffNadir);
% 
%     G = pattern_interp(offNadir);
%     EIRP = txPower * G;
%     PFD_stations(i) = EIRP / (4 * pi * R^2);
% 
%     PFD_dBW = 10 * log10(PFD_stations(i));
% 
%     fprintf('%s: off-nadir = %.2f° , PFD = %.2e W/m², = %.2f dBW, distancia = %.2f km\n', ...
%         stations{i,1}, offNadir, PFD_stations(i), PFD_dBW, R/1000);
% end
% 
% %% Calcular la PFD en una malla fina (región de interés)
% dLat = 0.1; % Paso en latitud
% dLon = 0.1; % Paso en longitud
% latVec = 25:dLat:57;  
% lonVec = -25:dLon:17;
% [lonGrid, latGrid] = meshgrid(lonVec, latVec);
% 
% PFD_map = zeros(size(latGrid));
% for idx = 1:numel(latGrid)
%     PFD_map(idx) = computePFD(latGrid(idx), lonGrid(idx), ...
%                               satLat, satLon, satAlt, ...
%                               txPower, pattern_interp, maxOffNadir);
% end
% 
% PFD_map_dBW = 10 * log10(PFD_map);
% 
% % puntos dentro del área de cobertura
% groundRange = slant_range * sind(off_nadir_lim);  % km
% 
% distances = arrayfun(@(lat, lon) haversine(satLat, satLon, lat, lon), latGrid, lonGrid);
% mask = distances <= groundRange * 1000;  % Solo puntos dentro del ground range
% 
% latInside = latGrid(mask);
% lonInside = lonGrid(mask);
% PFD_map_dBW_inside = PFD_map_dBW(mask);
% 
% %% Visualización
% figure
% ax = geoaxes;
% geobasemap(ax, 'satellite');
% geolimits(ax, [25 60], [-25 20]);
% hold(ax, 'on');
% 
% R_earth_km = 6378;
% angularRadius = rad2deg(groundRange / R_earth_km);  % grados
% [latCircle, lonCircle] = scircle1(satLat, satLon, angularRadius, R_earth_km, 'degrees');
% geoplot(ax, latCircle, lonCircle, 'r', 'LineWidth', 2);
% 
% % "Mapa de calor" de la PFD en dBW
% geoscatter(ax, latInside, lonInside, 100, PFD_map_dBW_inside, ...
%            'filled', 's', 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.2);
% colormap(ax, jet);
% cb = colorbar(ax);
% cb.Label.String = 'PFD (dBW)';
% 
% %estaciones base en el mapa
% for i = 1:numStations
%     gsLat = stations{i,2};
%     gsLon = stations{i,3};
%     geoplot(ax, gsLat, gsLon, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'red');
%     text(gsLat + 0.3, gsLon, stations{i,1}, 'Color', 'yellow', 'FontWeight', 'bold');
% end
% 
% title(ax, 'PFD en Tierra (dBW), Estaciones Base y Límite de Cobertura');
% hold(ax, 'off');
% 
% %% Funciones locales
% 
% function val = computePFD(lat, lon, satLat, satLon, satAlt, txPower, pattern_interp, maxOffNadir)
%     d = haversine(satLat, satLon, lat, lon);
%     R = sqrt(d^2 + satAlt^2);
%     elev = atand(satAlt / d);
%     offNadir = min(90 - elev, maxOffNadir);
%     G = pattern_interp(offNadir);
%     EIRP = txPower * G;
%     val = EIRP / (4 * pi * R^2);
% end
% 
% function d = haversine(lat1, lon1, lat2, lon2)
%     R = 6378e3;
%     dLat = deg2rad(lat2 - lat1);
%     dLon = deg2rad(lon2 - lon1);
%     a = sin(dLat/2)^2 + cosd(lat1) * cosd(lat2) * sin(dLon/2)^2;
%     c = 2 * atan2(sqrt(a), sqrt(1 - a));
%     d = R * c;
% end


%%
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------

% OPTIMIZAMOS NUMERO DE ELEMENTOS SIN ELEMENTO CENTRAL


% %%
% 
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
% fprintf('λ = %.4f m\n', lmbda_eff);
% 
% slant_range = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));
% fprintf('slant range = %.4f m\n', slant_range);
% 
% %%
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
% % ángulos diegrama
% peakAngle         = 54;
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
% azimuth = -100:100;
% Beam_dB = zeros(size(azimuth));
% 
% for k = 1:numel(azimuth)
%   th = azimuth(k);  x = abs(th);
%   if x <= (peakAngle-delta)
%     % --- valle isoflux geométrico ---
%     % interpolas la curva “exacta” y la elevas al nivel del pico
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
% 
% N1min = 2;
% N1max = 4;
% N2min = 4;
% N2max = 6;
% N3min = 4;
% N3max = 8;
% 
% % Initial guess for the number of elements in each ring
% N1 = 4; % Initial guess for ring 1
% N2 = 6; % Initial guess for ring 2
% N3 = 8; % Initial guess for ring 3
% 
% % Total number of elements (initial guess)
% N_total = N1 + N2 + N3;
% max_N_total = N1max + N2max + N3max;
% 
% 
% %%
% objfun = @(x) CostGA3AN(x, azimuth, Beam_d, lmbda_eff, intermediateAngle);
% 
% %%
% dim = 8 + 2*max_N_total;  % 5 variables geométricas + 2*16 = 35
% 
% x0 = zeros(dim,1);
% x0(1) = 0.5*lmbda_eff; % r1 
% x0(2) = 1.0*lmbda_eff; % r2
% x0(3) = 1.5*lmbda_eff; % r3
% x0(4) = 0; % desfaseRel inicial for ring 2
% x0(5) = 0; % desfaseRel inicial for ring 3
% x0(6) = N1; % Initial guess for N1
% x0(7) = N2; % Initial guess for N2
% x0(8) = N3; % Initial guess for N3
% x0(9 : 8+max_N_total) = 1; % Real part of weights
% x0(9+max_N_total : 8+2max_*N_total) = 0; % Imaginary part of weights
% 
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
% % Bounds for N1, N2, N3
% lb(6) = N1min; 
% ub(6) = N1max;
% lb(7) = N2min; 
% ub(7) = N2max;
% lb(8) = N3min; 
% ub(8) = N3max;
% % Bounds for real and imaginary parts of weights
% lb(9 : 8+max_N_total) = -50;
% ub(9 : 8+max_N_total) = 50;
% lb(9+max_N_total : 8+2*max_N_total) = -50;
% ub(9+max_N_total : 8+2*max_N_total) = 50;
% 
% % Specify integer constraints for N1, N2, N3
% IntCon = 6:8;
% 
% %%
% 
% % Arrays donde se irán guardando ambos costes
% centralHistory = [];
% externalHistory = [];
% 
% 
% ofcn = @(options,state,flag) costMonitor(...
%     options, state, flag, ...
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
% [x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [],IntCon, options);
% 
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
% % Extract results
% r1_opt = x_opt(1);
% r2_opt = x_opt(2);
% r3_opt = x_opt(3);
% desfaseRel_opt = x_opt(4);
% desfaseRel1_opt = x_opt(5);
% N1_opt = round(x_opt(6));
% N2_opt = round(x_opt(7));
% N3_opt = round(x_opt(8));
% 
% N_total_opt = N1_opt + N2_opt + N3_opt;
% 
% w_re_opt = x_opt(9 : 8+N_total_opt).';
% w_im_opt = x_opt(9+N_total_opt : 8+2*N_total_opt).';
% w_opt = w_re_opt + 1i*w_im_opt;
% 
% 
% %% 
% anglesRing1_opt = (0:N1_opt-1)*(360/N1_opt);  
% anglesRing2_opt = (0:N2_opt-1)*(360/N2_opt) + desfaseRel_opt;
% anglesRing3_opt = (0:N3_opt-1)*(360/N3_opt) + desfaseRel1_opt;
% 
% posRing1_opt = [r1_opt*cosd(anglesRing1_opt); r1_opt*sind(anglesRing1_opt)];
% posRing2_opt = [r2_opt*cosd(anglesRing2_opt); r2_opt*sind(anglesRing2_opt)];
% posRing3_opt = [r3_opt*cosd(anglesRing3_opt); r3_opt*sind(anglesRing3_opt)];
% elementPos_opt = [posRing1_opt, posRing2_opt, posRing3_opt];
% 
% 
% % matriz de apuntamiento
% stvmat_opt = zeros(N_total_opt, length(azimuth));
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
% %% Gráfico con leyenda
% figure; hold on;
% 
% hC1 = plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% hC2 = plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% hC3 = plot(r3_opt*cos(thetaC), r3_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% 
% % 2) Para cada elemento creamos su propio scatter, con DisplayName = "n: módulo ∠ fase"
% hElem = gobjects(N_total_opt,1);
% for k = 1:N_total_opt
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
% % 3) Ahora la leyenda la creamos sólo con los handles de los elementos
% legend(hElem, 'Location','eastoutside', 'FontSize',8);
% 
% 
% 
% %%
% function [state,options,optchanged] = costMonitor( ...
%         options, state, flag, ...
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
%         [~, cC, cE] = CostGA3AN( xBest, ...
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
%
%
% %% PARTE 2: SIMULACIÓN DE LA HUELLA (PFD) EN TIERRA
% %% Parámetros
% txPower = 10; % Watts
% 
% slant_range = 2121;
% off_nadir_lim = 54.33;
% 
% satAlt = 1000e3;
% satLat  = 40.4; % Latitud (aprox. Madrid)
% satLon  = -3.7; % Longitud (aprox. Madrid)
% 
% minElev = 20; % Elevación mínima de diseño (°)
% maxOffNadir = 90 - minElev; % Off-nadir máximo permitido (°)
% 
% azimuth = -90:90;                     
% pattern_interp = @(angle) interp1(azimuth, pattern_synth_opt, angle, 'linear', 'extrap');
% 
% %% Estaciones base
% stations = { 'Cadiz', 36.5, -6.3;
%              'Tarifa', 36, -5.6;
%              'Barcelona', 41.4, 2.2;
%              'La Coruna', 43.3, -8.4;
%              'Ponferrada', 42.9, -6.7;
%              'Madrid', 40.4, -3.7;
%              'Salamanca', 41, -5.96;
%              'Albacete', 39, -1.85;
%              'Paris', 48.8, 2.35;
%              'Londres', 51.7, 0};
% % NO ES PONFERRADA COMO TAL PERO ESTA CERCA (ES JUSTO EL LIMITE ENTRE
% % ASTURIAS Y LEON)
% 
% numStations = size(stations,1);
% PFD_stations = zeros(numStations,1);
% 
% fprintf('--- Estaciones Base ---\n');
% for i = 1:numStations
%     gsLat = stations{i,2};
%     gsLon = stations{i,3};
% 
%     dSurf = haversine(satLat, satLon, gsLat, gsLon);
%     R = sqrt(dSurf^2 + satAlt^2);
%     elev = atand(satAlt / dSurf);
%     offNadir = min(90 - elev, maxOffNadir);
% 
%     G = pattern_interp(offNadir);
%     EIRP = txPower * G;
%     PFD_stations(i) = EIRP / (4 * pi * R^2);
% 
%     PFD_dBW = 10 * log10(PFD_stations(i));
% 
%     fprintf('%s: off-nadir = %.2f° , PFD = %.2e W/m², = %.2f dBW, distancia = %.2f km\n', ...
%         stations{i,1}, offNadir, PFD_stations(i), PFD_dBW, R/1000);
% end
% 
% %% Calcular la PFD en una malla fina (región de interés)
% dLat = 0.1; % Paso en latitud
% dLon = 0.1; % Paso en longitud
% latVec = 25:dLat:57;  
% lonVec = -25:dLon:17;
% [lonGrid, latGrid] = meshgrid(lonVec, latVec);
% 
% PFD_map = zeros(size(latGrid));
% for idx = 1:numel(latGrid)
%     PFD_map(idx) = computePFD(latGrid(idx), lonGrid(idx), ...
%                               satLat, satLon, satAlt, ...
%                               txPower, pattern_interp, maxOffNadir);
% end
% 
% PFD_map_dBW = 10 * log10(PFD_map);
% 
% % puntos dentro del área de cobertura
% groundRange = slant_range * sind(off_nadir_lim);  % km
% 
% distances = arrayfun(@(lat, lon) haversine(satLat, satLon, lat, lon), latGrid, lonGrid);
% mask = distances <= groundRange * 1000;  % Solo puntos dentro del ground range
% 
% latInside = latGrid(mask);
% lonInside = lonGrid(mask);
% PFD_map_dBW_inside = PFD_map_dBW(mask);
% 
% %% Visualización
% figure
% ax = geoaxes;
% geobasemap(ax, 'satellite');
% geolimits(ax, [25 60], [-25 20]);
% hold(ax, 'on');
% 
% R_earth_km = 6378;
% angularRadius = rad2deg(groundRange / R_earth_km);  % grados
% [latCircle, lonCircle] = scircle1(satLat, satLon, angularRadius, R_earth_km, 'degrees');
% geoplot(ax, latCircle, lonCircle, 'r', 'LineWidth', 2);
% 
% % "Mapa de calor" de la PFD en dBW
% geoscatter(ax, latInside, lonInside, 100, PFD_map_dBW_inside, ...
%            'filled', 's', 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.2);
% colormap(ax, jet);
% cb = colorbar(ax);
% cb.Label.String = 'PFD (dBW)';
% 
% %estaciones base en el mapa
% for i = 1:numStations
%     gsLat = stations{i,2};
%     gsLon = stations{i,3};
%     geoplot(ax, gsLat, gsLon, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'red');
%     text(gsLat + 0.3, gsLon, stations{i,1}, 'Color', 'yellow', 'FontWeight', 'bold');
% end
% 
% title(ax, 'PFD en Tierra (dBW), Estaciones Base y Límite de Cobertura');
% hold(ax, 'off');
% 
% %% Funciones locales
% 
% function val = computePFD(lat, lon, satLat, satLon, satAlt, txPower, pattern_interp, maxOffNadir)
%     d = haversine(satLat, satLon, lat, lon);
%     R = sqrt(d^2 + satAlt^2);
%     elev = atand(satAlt / d);
%     offNadir = min(90 - elev, maxOffNadir);
%     G = pattern_interp(offNadir);
%     EIRP = txPower * G;
%     val = EIRP / (4 * pi * R^2);
% end
% 
% function d = haversine(lat1, lon1, lat2, lon2)
%     R = 6378e3;
%     dLat = deg2rad(lat2 - lat1);
%     dLon = deg2rad(lon2 - lon1);
%     a = sin(dLat/2)^2 + cosd(lat1) * cosd(lat2) * sin(dLon/2)^2;
%     c = 2 * atan2(sqrt(a), sqrt(1 - a));
%     d = R * c;
% end


%%
% SIN ELEMENTO CENTRAL
% CUANTIZACIÓN FASES
%%

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
% fprintf('λ = %.4f m\n', lmbda_eff);
% 
% slant_range = sqrt(r^2 + (r + h)^2 - 2*r*(r + h)*sin(elevacion + asin((r/(r+h))*cos(elevacion))));
% fprintf('slant range = %.4f m\n', slant_range);
% 
% %%
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
%     % --- valle isoflux geométrico ---
%     % interpolas la curva “exacta” y la elevas al nivel del pico
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
% N1 = 4;% Elementos en anillo 1
% N2 = 6;% Elementos en anillo 2
% N3 = 8;% Elementos en anillo 3
% N_total = N1 + N2 + N3;
% 
% %%
% objfun = @(x) CostGA3AN(x, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle);
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
% lb(6+N_total : 5+2*N_total) = 0;
% ub(6+N_total : 5+2*N_total) = 63;
% 
% % Configuro GA para variables enteras en las posiciones de índice de fase:
% IntCon = (6+N_total : 5+2*N_total);
% 
% %%
% 
% 
% % Arrays donde se irán guardando ambos costes
% centralHistory = [];
% externalHistory = [];
% 
% 
% ofcn = @(options,state,flag) costMonitor(...
%     options, state, flag, ...
%     N1, N2, N3, ...
%     azimuth, Beam_d, ...
%     lmbda_eff, intermediateAngle);
% 
% options = optimoptions('ga', ...
%     'Display','iter', ...
%     'FunctionTolerance',1e-8, ...
%     'ConstraintTolerance',1e-3, ...
%     'CrossoverFraction',0.7, ...
%     'PopulationSize',300, ...
%     'MaxGenerations',800, ...
%     'PlotFcn','gaplotbestf', ...
%     'OutputFcn',ofcn);
% 
% 
% [x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [], IntCon, options);
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
% % 1) Extrae magnitudes e índices de fase de la solución
% mags_opt   = x_opt(6 : 5+N_total);
% idx_opt    = x_opt(6+N_total : 5+2*N_total);  % enteros 0..63
% 
% % Asegurarnos de que son columnas
% mags_opt = mags_opt(:);
% idx_opt  = idx_opt(:);
% 
% % 2) Reconstruye el vector de pesos complejos
% phi_opt = (2*pi) * (idx_opt/64);                % fase en radianes
% w_opt   = mags_opt .* exp(1i * phi_opt);      % pesos finales
% 
% % 4) ¡FORZAMOS w_opt a columna de nuevo por si acaso!
% w_opt = w_opt(:);  % [N_total×1]
% 
% % calcula y dibuja las fases en grados
% phases_deg = idx_opt * (360/64);               % fase en grados [0,360)
% figure;
% stem(1:N_total, phases_deg, 'filled', 'LineWidth',1.5);
% xlabel('Element'); ylabel('Phase (°)');
% title('Cuantized phases (6 bits)');
% grid on; ylim([0 360]);
% 
% %% 
% % centro:
% anglesRing1_opt = (0:N1-1)*(360/N1);  
% anglesRing2_opt = (0:N2-1)*(360/N2) + desfaseRel_opt;
% anglesRing3_opt = (0:N3-1)*(360/N3) + desfaseRel1_opt;
% 
% posRing1_opt = [r1_opt*cosd(anglesRing1_opt); r1_opt*sind(anglesRing1_opt)];
% posRing2_opt = [r2_opt*cosd(anglesRing2_opt); r2_opt*sind(anglesRing2_opt)];
% posRing3_opt = [r3_opt*cosd(anglesRing3_opt); r3_opt*sind(anglesRing3_opt)];
% elementPos_opt = [posRing1_opt, posRing2_opt, posRing3_opt];
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
% %% Gráfico con leyenda
% figure; hold on;
% 
% % 1) Dibujamos los círculos de anillo SÓLO para referencia, SIN que salgan en la leyenda
% hC1 = plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% hC2 = plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% hC3 = plot(r3_opt*cos(thetaC), r3_opt*sin(thetaC), 'k--', 'HandleVisibility','off');
% 
% % 2) Para cada elemento creamos su propio scatter, con DisplayName = "n: módulo ∠ fase"
% hElem = gobjects(N_total,1);
% for k = 1:N_total
%     xk = elementPos_opt(1,k);
%     yk = elementPos_opt(2,k);
%     % Convertimos índice a grados
%     phase_deg = idx_opt(k) * (360/64);
%     % Preparamos la etiqueta
%     lbl = sprintf('%d: %.2f ∠ %.3f°', ...
%                   k, mags_opt(k), phase_deg);
% 
%     % scatter individual con DisplayName
%     hElem(k) = scatter(xk, yk, 80, 'b', 'filled', ...
%                        'DisplayName', lbl);
% 
%     % opcional: número en rojo encima
%     text(xk, yk, sprintf('%d',k), ...
%          'FontSize',8, 'FontWeight','bold', 'Color','r', ...
%          'HorizontalAlignment','center','VerticalAlignment','middle');
% end
% 
% axis equal; grid on;
% xlabel('X (m)'); ylabel('Y (m)');
% title('Array Elements Numbered with Magnitudes and Phase Indices');
% 
% legend(hElem, 'Location','eastoutside', 'FontSize',8);
% 
% 
% %%
% function [state,options,optchanged] = costMonitor( ...
%         options, state, flag, ...
%         N1, N2, N3, ...
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
%         [~, cC, cE] = CostGA3AN( xBest, ...
%                           N1, N2, N3, ...
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

%%
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------

%%
% CON ELEMENTO CENTRAL
% CUANTIZACIÓN FASES
%%

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
%     % interpolas la curva “exacta” y la elevas al nivel del pico
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
% Beam_sg = 10.^( Beam_dB_sg/20 );
% 
% Beam_d = Beam_sg;
% % ------------------------------------------------------------------------
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
% objfun = @(x) CostGA3AN(x, N0, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle);
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
% lb(6+N_total : 5+2*N_total) = 0;
% ub(6+N_total : 5+2*N_total) = 63;
% 
% % Configuro GA para variables enteras en las posiciones de índice de fase:
% IntCon = (6+N_total : 5+2*N_total);
% 
% %%
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
% [x_opt, fval, exitflag, output] = ga(objfun, dim, [], [], [], [], lb, ub, [],IntCon, options);
% 
% figure; hold on;
% gens = 0:length(centralHistory)-1;
% plot(gens, centralHistory, '-b','LineWidth',2);
% plot(gens, externalHistory, '-r','LineWidth',2);
% legend('Coste Central','Coste Externo','Location','best');
% xlabel('Generación'); ylabel('Coste');
% title('Evolución Costes Central vs Externo');
% grid on;
% 
% 
% r1_opt         = x_opt(1);
% r2_opt         = x_opt(2);
% r3_opt         = x_opt(3);
% desfaseRel_opt = x_opt(4);
% desfaseRel1_opt = x_opt(5);
% 
% 
% % 1) Extrae magnitudes e índices de fase de la solución
% mags_opt   = x_opt(6 : 5+N_total);
% idx_opt    = x_opt(6+N_total : 5+2*N_total);  % enteros 0..63
% 
% % Asegurarnos de que son columnas
% mags_opt = mags_opt(:);
% idx_opt  = idx_opt(:);
% 
% % 2) Reconstruye el vector de pesos complejos
% phi_opt = (2*pi) * (idx_opt/64);                % fase en radianes
% w_opt   = mags_opt .* exp(1i * phi_opt);      % pesos finales
% 
% % 4) ¡FORZAMOS w_opt a columna de nuevo por si acaso!
% w_opt = w_opt(:);  % [N_total×1]
% 
% % 3) (Opcional) calcula y dibuja las fases en grados
% phases_deg = idx_opt * (360/64);               % fase en grados [0,360)
% figure;
% stem(1:N_total, phases_deg, 'filled', 'LineWidth',1.5);
% xlabel('Element'); ylabel('Phase (°)');
% title('Cuantized phases (6 bits)');
% grid on; ylim([0 360]);
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
% %% Gráfico con leyenda
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
%     % Convertimos índice a grados
%     phase_deg = idx_opt(k) * (360/64);
%     % Preparamos la etiqueta
%     lbl = sprintf('%d: %.2f ∠ %.0f°', ...
%                   k, mags_opt(k), phase_deg);
% 
%     % scatter individual con DisplayName
%     hElem(k) = scatter(xk, yk, 80, 'b', 'filled', ...
%                        'DisplayName', lbl);
% 
%     % opcional: número en rojo encima
%     text(xk, yk, sprintf('%d',k), ...
%          'FontSize',8, 'FontWeight','bold', 'Color','r', ...
%          'HorizontalAlignment','center','VerticalAlignment','middle');
% end
% 
% axis equal; grid on;
% xlabel('X (m)'); ylabel('Y (m)');
% title('Array Elements Numbered with Magnitudes and Phase Indices');
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
%         [~, cC, cE] = CostGA3AN( xBest, ...
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