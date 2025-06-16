%% RESTRICCIONES 1/2/3/4/5/6/7


function cost = macroArrayCostGAPRUEBA(x_wb, z, n_elem, elementPos, ringIndex, azimuth, Beam_d, lmbda_eff, intermediateAngle)
    % Función de coste para macro array con geometría fija,
    % optimizando solo pesos complejos y máscara on/off.

    N_total = length(ringIndex);
    % Extraer pesos y máscara
    w_re = x_wb(1:N_total);
    w_im = x_wb(N_total+1:2*N_total);
    b_vec = x_wb(2*N_total+1:3*N_total);

    % Construir peso complejo bruto
    w_raw = w_re + 1i*w_im;

    nBits      = 6;
    Nlevels    = 2^nBits;
    phaseStep  = 360 / Nlevels;

    % Calcular desfase continuo por elemento en radianes
    % fases en grados: -180 + n_elem*5.625
    phase_deg = 0 + n_elem * phaseStep;
    phase_rad = deg2rad(phase_deg);

    % Aplicar máscara, desfase y peso
    w_eff = (b_vec .* w_raw) .* exp(1i * phase_rad);  % 1×N_total

    % Cálculo de la matriz de apuntamiento (steering vector)
    numAngles = length(azimuth);
    stvmat = zeros(N_total, numAngles);
    for i = 1:numAngles
        phi = azimuth(i);
        stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)' * cosd(phi) + elementPos(2,:)' * sind(phi) ) / lmbda_eff );
    end

    % Patrón sintetizado: usar w_eff.' para asegurar 1×N_total
    pattern_synth = abs(w_eff * stvmat);   % 1×numAngles
    pattern_synth_dB = 20 * log10(pattern_synth);


    % Definir la región central con ángulo intermedio (40°)
    % intermediateAngle = 70; % grados
    central_idx = find(abs(azimuth) <= intermediateAngle);
    external_idx = find(abs(azimuth) > intermediateAngle);

    % Tolerancia para ubicación de puntos críticos
    tol = 1;
    W = ones(size(azimuth));

    % Pesos en puntos específicos
    % w_valley = 20; w_peak = 20;
    % w_a = 13; w_b = 13; w_c = 13;
    % w_d = 10; w_e = 10; w_f = 10; w_g = 10;
    % 
    % % Asignación de pesos
    % W(abs(azimuth) < tol) = w_valley;
    % W(abs(azimuth - 54) < tol | abs(azimuth + 54) < tol) = w_peak;
    % W(abs(azimuth - 5)  < tol | abs(azimuth + 5)  < tol) = w_a;
    % W(abs(azimuth - 10) < tol | abs(azimuth + 10) < tol) = w_b;
    % W(abs(azimuth - 15) < tol | abs(azimuth + 15) < tol) = w_c;
    % W(abs(azimuth - 25) < tol | abs(azimuth + 25) < tol) = w_d;
    % W(abs(azimuth - 30) < tol | abs(azimuth + 30) < tol) = w_e;
    % W(abs(azimuth - 35) < tol | abs(azimuth + 35) < tol) = w_f;
    % W(abs(azimuth - 40) < tol | abs(azimuth + 40) < tol) = w_g;
    w_valley = 1;
    w_peak   = 1;
    w_a = 1;
    w_b = 1;
    w_c = 1;
    w_d = 1;
    w_e = 1;
    w_f = 1;
    w_g = 1;
    w_h = 1;
    w_i = 1;
    w_j = 1;

    % Para el valle (ángulo 0°)
    W(abs(azimuth - 0) < tol) = w_valley;
    W(abs(azimuth + 0) < tol) = w_valley;
    % W(abs(azimuth) < tol) = w_valley;
    % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
    W(abs(azimuth - 54) < tol) = w_peak;
    W(abs(azimuth + 54) < tol) = w_peak;

    W(abs(azimuth - 5) < tol) = w_a;
    W(abs(azimuth + 5) < tol) = w_a;

    W(abs(azimuth - 3) < tol) = 1;
    W(abs(azimuth + 3) < tol) = 1;

    W(abs(azimuth - 10) < tol) = w_b;
    W(abs(azimuth + 10) < tol) = w_b;

    W(abs(azimuth - 12) < tol) = 1;
    W(abs(azimuth + 12) < tol) = 1;

    W(abs(azimuth - 15) < tol) = w_c;
    W(abs(azimuth + 15) < tol) = w_c;

    W(abs(azimuth - 17) < tol) = 1;
    W(abs(azimuth + 17) < tol) = 1;

    W(abs(azimuth - 20) < tol) = w_h;
    W(abs(azimuth + 20) < tol) = w_h;

    W(abs(azimuth - 25) < tol) = w_d;
    W(abs(azimuth + 25) < tol) = w_d;

    W(abs(azimuth - 30) < tol) = w_e;
    W(abs(azimuth + 30) < tol) = w_e;

    W(abs(azimuth - 35) < tol) = w_f;
    W(abs(azimuth + 35) < tol) = w_f;

    W(abs(azimuth - 40) < tol) = w_g;
    W(abs(azimuth + 40) < tol) = w_g;

    W(abs(azimuth - 45) < tol) = 1;
    W(abs(azimuth + 45) < tol) = 1;

    W(abs(azimuth - 50) < tol) = 1;
    W(abs(azimuth + 50) < tol) = 1;

    W(abs(azimuth - 60) < tol) = w_i;
    W(abs(azimuth + 60) < tol) = w_i;

    W(abs(azimuth - 65) < tol) = w_j;
    W(abs(azimuth + 65) < tol) = w_j;

    % Ponderaciones por región
    W_central = W(central_idx);
    W_external = W(external_idx);

    % Cálculo de error en región central y externa
    error_central = pattern_synth(central_idx) - Beam_d(central_idx);
    cost_central = sqrt(sum((W_central .* error_central).^2));

    error_external = pattern_synth(external_idx) - Beam_d(external_idx);
    cost_external = sqrt(sum((W_external .* error_external).^2));

    % Costo total con ponderación
    cost = cost_central + cost_external;
end




%%
%%
%%
%%
% RESTRICCIONES 1/2/3/4/5/7

% function cost = macroArrayCostGAPRUEBA(x_wb, z, elementPos, ringIndex, azimuth, Beam_d, lmbda_eff, intermediateAngle)
%     % Función de coste para macro array con geometría fija,
%     % optimizando solo pesos complejos y máscara on/off.
% 
%     N_total = length(ringIndex);
%     % Extraer pesos y máscara
%     w_re = x_wb(1:N_total);
%     w_im = x_wb(N_total+1:2*N_total);
%     b_vec = x_wb(2*N_total+1:3*N_total);
% 
%     % Construir peso complejo bruto
%     w_raw = w_re + 1i*w_im;
% 
%     % Aplicar máscara, desfase y peso
%     w_eff = (b_vec .* w_raw);  % 1×N_total
% 
%     % Cálculo de la matriz de apuntamiento (steering vector)
%     % numAngles = length(azimuth);
%     % stvmat = zeros(N_total, numAngles);
%     % for i = 1:numAngles
%     %     phi = azimuth(i);
%     %     stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)' * cosd(phi) + elementPos(2,:)' * sind(phi) ) / lmbda_eff );
%     % end
% 
%     cos_phi = cosd(azimuth);
%     sin_phi = sind(azimuth);
%     % xy = elementPos' * [cos_phi; sin_phi];
%     % stvmat = exp(1i * 2*pi * xy / lmbda_eff);
% 
%     phase = elementPos(1,:).' * cos_phi + elementPos(2,:).' * sin_phi; % [N_total x numAngles]
%     stvmat = exp(1i * 2*pi * phase / lmbda_eff);   % [N_total x numAngles]
% 
%     % Patrón sintetizado: usar w_eff.' para asegurar 1×N_total
%     pattern_synth = abs(w_eff * stvmat);   % 1×numAngles
%     pattern_synth_dB = 20 * log10(pattern_synth);
% 
% 
%     % Definir la región central con ángulo intermedio (40°)
%     % intermediateAngle = 70; % grados
%     central_idx = find(abs(azimuth) <= intermediateAngle);
%     external_idx = find(abs(azimuth) > intermediateAngle);
% 
%     % Tolerancia para ubicación de puntos críticos
%     % tol = 0.1;
%     tol = 1;
%     W = ones(size(azimuth));
% 
%     % Pesos en puntos específicos
%     % w_valley = 20; % por ejemplo, darle tres veces más importancia al valle en 0°
%     % w_peak = 20;   % y tres veces más importancia a los picos en ±peakAngle
%     % 
%     % w_a = 13;
%     % w_b = 13;
%     % w_c = 13;
%     % w_d = 10;
%     % w_e = 10;
%     % w_f = 10;
%     % w_g = 10;
%     % w_h = 10;
%     % w_i = 13;
%     % 
%     % % Para el valle (ángulo 0°)
%     % W(abs(azimuth) < tol) = w_valley;
%     % % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
%     % W(abs(azimuth - 54) < tol) = w_peak;
%     % W(abs(azimuth + 54) < tol) = w_peak;
%     % 
%     % W(abs(azimuth - 5) < tol) = w_a;
%     % W(abs(azimuth + 5) < tol) = w_a;
%     % 
%     % W(abs(azimuth - 10) < tol) = w_b;
%     % W(abs(azimuth + 10) < tol) = w_b;
%     % 
%     % W(abs(azimuth - 15) < tol) = w_c;
%     % W(abs(azimuth + 15) < tol) = w_c;
%     % 
%     % W(abs(azimuth - 20) < tol) = w_h;
%     % W(abs(azimuth + 20) < tol) = w_h;
%     % 
%     % W(abs(azimuth - 25) < tol) = w_d;
%     % W(abs(azimuth + 25) < tol) = w_d;
%     % 
%     % W(abs(azimuth - 30) < tol) = w_e;
%     % W(abs(azimuth + 30) < tol) = w_e;
%     % 
%     % W(abs(azimuth - 35) < tol) = w_f;
%     % W(abs(azimuth + 35) < tol) = w_f;
%     % 
%     % W(abs(azimuth - 40) < tol) = w_g;
%     % W(abs(azimuth + 40) < tol) = w_g;
%     % 
%     % W(abs(azimuth - 60) < tol) = w_i;
%     % W(abs(azimuth + 60) < tol) = w_i;
% 
% 
%     w_valley = 1;
%     w_peak   = 1;
%     w_a = 1;
%     w_b = 1;
%     w_c = 1;
%     w_d = 1;
%     w_e = 1;
%     w_f = 1;
%     w_g = 1;
%     w_h = 1;
%     w_i = 1;
%     w_j = 1;
% 
%     % Para el valle (ángulo 0°)
%     W(abs(azimuth - 0) < tol) = w_valley;
%     W(abs(azimuth + 0) < tol) = w_valley;
%     % W(abs(azimuth) < tol) = w_valley;
%     % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
%     W(abs(azimuth - 54) < tol) = w_peak;
%     W(abs(azimuth + 54) < tol) = w_peak;
% 
%     W(abs(azimuth - 5) < tol) = w_a;
%     W(abs(azimuth + 5) < tol) = w_a;
% 
%     W(abs(azimuth - 3) < tol) = 1;
%     W(abs(azimuth + 3) < tol) = 1;
% 
%     W(abs(azimuth - 10) < tol) = w_b;
%     W(abs(azimuth + 10) < tol) = w_b;
% 
%     W(abs(azimuth - 12) < tol) = 1;
%     W(abs(azimuth + 12) < tol) = 1;
% 
%     W(abs(azimuth - 15) < tol) = w_c;
%     W(abs(azimuth + 15) < tol) = w_c;
% 
%     W(abs(azimuth - 17) < tol) = 1;
%     W(abs(azimuth + 17) < tol) = 1;
% 
%     W(abs(azimuth - 20) < tol) = w_h;
%     W(abs(azimuth + 20) < tol) = w_h;
% 
%     W(abs(azimuth - 25) < tol) = w_d;
%     W(abs(azimuth + 25) < tol) = w_d;
% 
%     W(abs(azimuth - 30) < tol) = w_e;
%     W(abs(azimuth + 30) < tol) = w_e;
% 
%     W(abs(azimuth - 35) < tol) = w_f;
%     W(abs(azimuth + 35) < tol) = w_f;
% 
%     W(abs(azimuth - 40) < tol) = w_g;
%     W(abs(azimuth + 40) < tol) = w_g;
% 
%     W(abs(azimuth - 45) < tol) = 1;
%     W(abs(azimuth + 45) < tol) = 1;
% 
%     W(abs(azimuth - 50) < tol) = 1;
%     W(abs(azimuth + 50) < tol) = 1;
% 
%     W(abs(azimuth - 60) < tol) = w_i;
%     W(abs(azimuth + 60) < tol) = w_i;
% 
%     W(abs(azimuth - 65) < tol) = w_j;
%     W(abs(azimuth + 65) < tol) = w_j;
% 
% 
%     % Ponderaciones por región
%     W_central = W(central_idx);
%     W_external = W(external_idx);
% 
%     % Cálculo de error en región central y externa
%     error_central = pattern_synth(central_idx) - Beam_d(central_idx);
%     cost_central = sqrt(sum((W_central .* error_central).^2));
% 
%     % error_external = pattern_synth(external_idx) - Beam_d(external_idx);
%     % cost_external = sqrt(sum((W_external .* error_external).^2));
% 
%     % Costo total con ponderación
%     cost = cost_central;
% end