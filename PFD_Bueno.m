clear; clc; close all;


spacing   = 1;                          
lats      = -90:spacing:90;            
lons      = -180:spacing:180;          
[lonG, latG] = meshgrid(lons, lats);    % lonG=columnas, latG=filas
%% Parámetros básicos
Re       = 6378;      % km
h_km     = 1000;      % km
Emin_deg = 20;        % °
P_t = 1;

% Elevación mínima en radianes
Emin_rad = deg2rad(Emin_deg);

% Slant_ref (para la elevación mínima de 20°)
slant_ref = -Re*sin(Emin_rad) + sqrt((Re+h_km)^2 - (Re*cos(Emin_rad))^2);
% → ≃ 2121.45 km

%% Calculo φ_max y el off-axis real en ese punto (peakAngle)
phi_max = acos( ((Re+h_km)^2 + Re^2 - slant_ref^2) / (2*Re*(Re+h_km)) );
peakAngle_deg = rad2deg( atan2( Re*sin(phi_max), (Re+h_km)-Re*cos(phi_max) ) );
% → ≃ 54.3241°

%% Construyo mi vector de d_phi y theta_iso
N     = 2001;
phi   = linspace(-phi_max, +phi_max, N);            % rad
d_phi = sqrt((Re+h_km)^2 + Re^2 - 2*Re*(Re+h_km).*cos(phi));  % km
theta_iso_deg = rad2deg( atan2( Re*sin(phi), (Re+h_km)-Re*cos(phi) ) );  % °

%% Puntos donde quiero estaciones
anglesOff = [0, 10, 20, 30, 40, peakAngle_deg];   % incluyo el verdadero peakAngle

%% Interpolo el slant-range “exacto” para cada θ_off solicitado
slant_km = interp1(theta_iso_deg, d_phi, anglesOff, 'linear', NaN);
slant_m  = slant_km * 1e3;
range_m = d_phi * 1e3;  % pasa km→m para todo φ

%% Interpolo la ganancia en esas estaciones (en dBi → lineal)
D        = load('resultado_1.txt');
thetaVec = D(:,1);
GdBiVec  = D(:,2);
GlinGS   = 10.^( interp1(thetaVec, GdBiVec, anglesOff, 'linear', -Inf) / 10 );
GlinVec = 10.^(GdBiVec/10);

%% Calculo la PFD en W/m² y dBW/m²
pfd_W    = P_t * GlinGS ./ (4*pi*slant_m.^2);
pfd_dBWm = 10*log10(pfd_W);

%% Muestro resultados
fprintf('\n  off-axis | Slant [km] |    PFD [W/m²]   | PFD [dBW/m²]\n');
fprintf('-----------------------------------------------------------\n');
for k = 1:numel(anglesOff)
    fprintf('     %6.2f°   |   %8.2f  |  %9.3e  |   %6.2f\n', ...
        anglesOff(k), slant_km(k), pfd_W(k), pfd_dBWm(k));
end

%% 8) Plot de PFD vs Off–axis

figure; 
% eje izquierdo: PFD en W/m²
yyaxis left
plot(anglesOff, pfd_W, 'o-', 'LineWidth', 1.5, 'MarkerSize',8);
ylabel('PFD (W/m²)');
ylim([0, max(pfd_W)*1.1]);  % un poco de margen

% eje derecho: PFD en dBW/m²
yyaxis right
plot(anglesOff, pfd_dBWm, 's--', 'LineWidth', 1.5, 'MarkerSize',8);
ylabel('PFD (dBW/m²)');
ylim([min(pfd_dBWm)*1.1, 0]);  % ajusta para que se vea completo

% configuración común
xlabel('Off–axis angle rho (°)');
grid on;
title('PFD vs Off–axis angle for ground stations');
legend('PFD [W/m²]', 'PFD [dBW/m²]', 'Location','best');

%%
% máscara elevación mínima
elevMap = asind( ((Re+h_km).*cos(phi) - Re) ./ d_phi );
maskElev = elevMap < Emin_deg;


%% 4) Calculo de fluxMap (W/m²) sobre toda la grilla

% 4.1) Defino el subsatélite en (lat0, lon0). Aquí tomo (0,0) como ejemplo
lat0 = 0;
lon0 = 0;

% 4.2) Calculo el ángulo central φ para cada celda del mapa
cos_phi_map = sind(lat0).*sind(latG) + cosd(lat0).*cosd(latG).*cosd(lonG-lon0);
phi_map     = acos(cos_phi_map);    % radianes

% 4.3) Slant-range en km para cada φ
slant_map_km = sqrt((Re+h_km)^2 + Re^2 - 2*Re*(Re+h_km).*cos(phi_map));
range_map    = slant_map_km * 1e3;  % m

% 4.4) Calculo el off-axis angle θ_map en grados
theta_map = rad2deg( atan2( Re*sin(phi_map), (Re+h_km) - Re*cos(phi_map) ) );

Glin_map = interp1(thetaVec, GlinVec, theta_map, 'linear', 0);

% 4.6) máscara de elevación mínima (en grados)
elev_map = asind( ((Re+h_km).*cos(phi_map) - Re) ./ slant_map_km );
Glin_map(elev_map < Emin_deg) = NaN;

% 4.7)calculo la PFD en W/m² sobre toda la grilla
fluxMap = P_t .* Glin_map ./ (4*pi * range_map.^2);


%% 7) Mapa de huella en dBW/m²
figure;
ax = worldmap([-30 30],[-30 30]);
setm(ax,'MapProjection','eqdcylin','Frame','on','Grid','on',...
    'MLineLocation',10,'PLineLocation',10,...
    'MeridianLabel','on','ParallelLabel','on','FontSize',12);

% dibujo tierra
land = shaperead('landareas','UseGeoCoords',true);
geoshow(ax, land,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.5 0.5 0.5]);

% huella dBW/m² encima
h = surfacem(latG, lonG, fluxMap);
colormap(flipud(parula(256)));
c = colorbar('eastoutside');
c.Label.String = 'PFD (W/m^2)';
c.FontSize      = 12;

title('Isoflux pattern (W/m^2)', 'FontSize',14);


%% 7) Mapa de huella en dBW/m²
fluxMap_dBW = 10*log10(fluxMap);

figure;
ax = worldmap([-30 30],[-30 30]);
setm(ax,'MapProjection','eqdcylin','Frame','on','Grid','on',...
    'MLineLocation',10,'PLineLocation',10,...
    'MeridianLabel','on','ParallelLabel','on','FontSize',12);

% dibujo tierra
land = shaperead('landareas','UseGeoCoords',true);
geoshow(ax, land,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.5 0.5 0.5]);

% huella dBW/m² encima
h = surfacem(latG, lonG, fluxMap_dBW);
colormap(flipud(parula(256)));
c = colorbar('eastoutside');
c.Label.String = 'PFD (dBW/m^2)';
c.FontSize      = 12;

title('Isoflux pattern (dBW/m^2)', 'FontSize',14);