
%%
% RUNEAR SIEMPRE PRIMERO CON LA REFERENCIA 20º, 1000KM Y 15DB
% Si quiero ver como cambian las elevaciones y ganancias con otra altura,
% tengo que ir al otro programa y runearlo y fijarme que ganancia hay a la
% altura que yo quiero, y actualizar aquí h y G_ref_dB

R_E       = 6378;        % km, radio de la Tierra
h         = 1200;        % km, altura del satélite
elevacion = deg2rad(20); % en grados
r_ref = sqrt(R_E^2 + (R_E + h)^2 - 2*R_E*(R_E + h)*sin(elevacion + asin((R_E/(R_E+h))*cos(elevacion))));        % km

G_ref_dB = 16.272; % Ganancia de referencia en dB

elevaciones = 5:5:90; % Elevaciones mínimas en grados

% Inicializamos arrays
ganancias_dB = zeros(size(elevaciones));
slant_range = zeros(size(elevaciones));

% Cálculo para cada elevación
for i = 1:length(elevaciones)
    elev = elevaciones(i);
    elev_rad = deg2rad(elev);
    
    % alpha = pi/2 - elev_rad; % ángulo complementario
    slant_range(i) = sqrt(R_E^2 + (R_E + h)^2 - 2*R_E*(R_E + h)*sin(elev_rad + asin((R_E/(R_E+h))*cos(elev_rad))));
    
    % Cálculo de la ganancia necesaria
    ganancias_dB(i) = G_ref_dB + 20 * log10(slant_range(i) / r_ref);
end

% Gráfica
figure;
plot(elevaciones, ganancias_dB, 'o-', 'LineWidth', 2);
xlabel('Minimum elevation (°)');
ylabel('Gain (dB)');
title(['With satellite at ', num2str(h), 'Km']);
grid on;

% Mostrar tabla
disp(table(elevaciones', slant_range', ganancias_dB', ...
    'VariableNames', {'Elevation_degrees', 'SlantRange_km', 'Gain_dB'}));