
%%
% RUNEAR SIEMPRE PRIMERO CON LA REFERENCIA 20º, 1000KM Y 15DB
% Si quiero ver como cambian las alturas y ganancias con otra elevación,
% tengo que ir al otro programa y runearlo y fijarme que ganancia hay a la
% elevación que yo quiero, y actualizar aquí h_ref, elevacion y G_ref_dB

R_E = 6378; % Radio de la Tierra (km)
elevacion = 20; % Elevación mínima en grados (puedes cambiarlo)
elev_rad = deg2rad(elevacion); % Radianes

% Ganancia y slant range de referencia (a 1000 km, 2121 km slant range y 15 dB ganancia)
G_ref_dB = 15; % dB
h_ref = 1000; % km
r_ref = sqrt(R_E^2 + (R_E + h_ref)^2 - 2*R_E*(R_E + h_ref)*sin(elev_rad + asin((R_E/(R_E+h_ref))*cos(elev_rad)))); % km (slant range de referencia)

% Vector de alturas
alturas = 500:100:1500; % km

% Inicialización
ganancias_dB = zeros(size(alturas));
slant_ranges = zeros(size(alturas));

for i = 1:length(alturas)
    h = alturas(i);
    
    % Cálculo del slant range
    % alpha = pi/2 - elev_rad; % ángulo complementario
    slant_ranges(i) = sqrt(R_E^2 + (R_E + h)^2 - 2*R_E*(R_E + h)*sin(elev_rad + asin((R_E/(R_E+h))*cos(elev_rad))));
    
    % Ganancia necesaria
    ganancias_dB(i) = G_ref_dB + 20 * log10(slant_ranges(i) / r_ref);
end

% Gráfica
figure;
plot(alturas, ganancias_dB, 'o-', 'LineWidth', 2);
xlabel('Satellite altitud (km)');
ylabel('Gain (dB)');
title(['Set for minimum elevation of ', num2str(elevacion), '°']);
grid on;

% Mostrar tabla
disp(table(alturas', slant_ranges', ganancias_dB', ...
    'VariableNames', {'Altitud_km', 'SlantRange_km', 'Gain_dB'}));