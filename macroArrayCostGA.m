function cost = macroArrayCostGA(x, elementPos, ringIndex, azimuth, Beam_d, lmbda_eff)
    % Función de coste para macro array con 20 anillos,
    % optimizando desfases por anillo, pesos complejos y máscara on/off.
    
    numRings = max(ringIndex);
    N_total = length(ringIndex);
    
    % Extraer los desfases para cada anillo (en grados)
    phase_offsets = x(1:numRings);
    
    % Extraer los pesos complejos: partes real e imaginaria
    w_re = x(numRings+1 : numRings+N_total).';
    w_im = x(numRings+N_total+1 : numRings+2*N_total).';
    w_raw = w_re + 1i * w_im;
    

    b_start = numRings + 2*N_total + 1;
    b = x(b_start : b_start + N_total - 1).';    % 1×N_total máscara on/off
    
    % Calcular los pesos efectivos: se aplica el desfase de cada anillo y la máscara
    w_eff = zeros(1, N_total);
    for idx = 1:N_total
        ring_num = ringIndex(idx);
        w_eff(idx) = b(idx) * w_raw(idx) * exp(1i * deg2rad(phase_offsets(ring_num)));
    end
    
    % Cálculo de la matriz de apuntamiento (steering vector)
    numAngles = length(azimuth);
    stvmat = zeros(N_total, numAngles);
    for i = 1:numAngles
        phi = azimuth(i);
        stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)' * cosd(phi) + elementPos(2,:)' * sind(phi) ) / lmbda_eff );
    end
    
    % Patrón sintetizado
    pattern_synth = abs(w_eff * stvmat);
    pattern_synth_dB = 20 * log10(pattern_synth);
    
    % Definir la región central usando un ángulo intermedio (40°)
    intermediateAngle = 70; % grados
    central_idx = find(abs(azimuth) <= intermediateAngle);
    external_idx = find(abs(azimuth) > intermediateAngle);
    
    tol = 0.1;  % Tolerancia en grados para asignar pesos
    % Vector de ponderación para puntos críticos del patrón
    W = ones(size(azimuth));
    
    % Valores de peso en puntos críticos
    w_valley = 10;
    w_peak   = 10;
    w_a = 6;
    w_b = 6;
    w_c = 6;
    w_d = 6;
    w_e = 6;
    w_f = 6;
    w_g = 6;
    w_h = 6;
    w_i = 6;
    w_j = 6;

    % Para el valle (ángulo 0°)
    W(abs(azimuth - 0) < tol) = w_valley;
    W(abs(azimuth + 0) < tol) = w_valley;
    % W(abs(azimuth) < tol) = w_valley;
    % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
    W(abs(azimuth - 54) < tol) = w_peak;
    W(abs(azimuth + 54) < tol) = w_peak;

    W(abs(azimuth - 5) < tol) = w_a;
    W(abs(azimuth + 5) < tol) = w_a;

    W(abs(azimuth - 3) < tol) = 7;
    W(abs(azimuth + 3) < tol) = 7;

    W(abs(azimuth - 10) < tol) = w_b;
    W(abs(azimuth + 10) < tol) = w_b;

    W(abs(azimuth - 12) < tol) = 7;
    W(abs(azimuth + 12) < tol) = 7;

    W(abs(azimuth - 15) < tol) = w_c;
    W(abs(azimuth + 15) < tol) = w_c;

    W(abs(azimuth - 17) < tol) = 7;
    W(abs(azimuth + 17) < tol) = 7;

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

    W(abs(azimuth - 45) < tol) = 7;
    W(abs(azimuth + 45) < tol) = 7;

    W(abs(azimuth - 50) < tol) = 7;
    W(abs(azimuth + 50) < tol) = 7;

    W(abs(azimuth - 60) < tol) = w_i;
    W(abs(azimuth + 60) < tol) = w_i;

    W(abs(azimuth - 65) < tol) = w_j;
    W(abs(azimuth + 65) < tol) = w_j;
    
    % Separar la ponderación para las regiones central y externa
    W_central = W(central_idx);
    W_external = W(external_idx);
    
    % Error en la región central
    error_central = pattern_synth(central_idx) - Beam_d(central_idx);
    cost_central = sqrt(sum((W_central .* error_central).^2));
    
    mask_threshold_dB = -5;
    penalty_factor = 0.1;  % Factor de penalización, ajustable
    % penalty = 0;
    cost_external_1 = 0;
    for i = 1:length(external_idx)
        idx = external_idx(i);
        if pattern_synth_dB(idx) > mask_threshold_dB
            overshoot = pattern_synth_dB(idx) - mask_threshold_dB;
            cost_external_1 = cost_external_1 + penalty_factor * overshoot^2;
        end
    end
    
    cost = 2.5 * cost_central + cost_external_1;
end
