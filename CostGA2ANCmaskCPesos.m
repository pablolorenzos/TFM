function cost = CostGA2ANCmaskCPesos(x, N1, N2, azimuth, Beam_d, lmbda, intermediateAngle)

    x = x(:);
    N_total = N1 + N2;
    
    r1 = x(1);
    r2 = x(2);
    desfaseRel = x(3);
    
    w_re = x(4 : 3+N_total).';
    w_im = x(4+N_total : 3+2*N_total).';
    w = w_re + 1i*w_im;
    w = w(:);  % [N_total x 1]
    
    anglesRing1 = (0:N1-1) * (360/N1);
    anglesRing2 = (0:N2-1) * (360/N2) + desfaseRel;
    
    posRing1 = [r1 * cosd(anglesRing1); r1 * sind(anglesRing1)];  % 2 x N1
    posRing2 = [r2 * cosd(anglesRing2); r2 * sind(anglesRing2)];  % 2 x N2
    elementPos = [posRing1, posRing2];  % 2 x N_total
    
    % stvmat [N_total x numAngles]
    numAngles = length(azimuth);
    stvmat = zeros(N_total, numAngles);
    for i = 1:numAngles
        phi = azimuth(i);
        stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda );
    end
    
    % w' es de tamaño [1 x N_total] y stvmat es [N_total x numAngles], dando un resultado de [1 x numAngles]
    pattern_synth = abs(w' * stvmat);  % Vector fila de tamaño [1 x numAngles]


    W = ones(size(azimuth)); % peso base 1 para todos los ángulos
    tol = 0.5;  % Tolerancia en grados para identificar el "punto exacto"
    % valores de peso en los puntos críticos:
    w_valley = 20;
    w_peak = 20;

    w_a = 13;
    w_b = 13;
    w_c = 13;
    w_d = 10;
    w_e = 10;
    w_f = 10;
    w_g = 10;
    w_h = 10;
    w_i = 13;

    % Para el valle (ángulo 0°)
    W(abs(azimuth) < tol) = w_valley;
    % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
    W(abs(azimuth - 54) < tol) = w_peak;
    W(abs(azimuth + 54) < tol) = w_peak;

    W(abs(azimuth - 5) < tol) = w_a;
    W(abs(azimuth + 5) < tol) = w_a;

    W(abs(azimuth - 7) < tol) = 10;
    W(abs(azimuth + 7) < tol) = 10;

    W(abs(azimuth - 10) < tol) = w_b;
    W(abs(azimuth + 10) < tol) = w_b;

    W(abs(azimuth - 15) < tol) = w_c;
    W(abs(azimuth + 15) < tol) = w_c;

    W(abs(azimuth - 17) < tol) = 10;
    W(abs(azimuth + 17) < tol) = 10;

    W(abs(azimuth - 20) < tol) = w_h;
    W(abs(azimuth + 20) < tol) = w_h;

    W(abs(azimuth - 25) < tol) = w_d;
    W(abs(azimuth + 25) < tol) = w_d;

    W(abs(azimuth - 30) < tol) = w_e;
    W(abs(azimuth + 30) < tol) = w_e;

    W(abs(azimuth - 32) < tol) = 10;
    W(abs(azimuth + 32) < tol) = 10;

    W(abs(azimuth - 35) < tol) = w_f;
    W(abs(azimuth + 35) < tol) = w_f;

    W(abs(azimuth - 40) < tol) = w_g;
    W(abs(azimuth + 40) < tol) = w_g;

    W(abs(azimuth - 45) < tol) = 10;
    W(abs(azimuth + 45) < tol) = 10;

    W(abs(azimuth - 50) < tol) = 10;
    W(abs(azimuth + 50) < tol) = 10;

    W(abs(azimuth - 60) < tol) = w_i;
    W(abs(azimuth + 60) < tol) = w_i;

    W(abs(azimuth - 65) < tol) = 10;
    W(abs(azimuth + 65) < tol) = 10;
    

    errorVec = pattern_synth - Beam_d; % error en cada ángulo entre el patrón sintetizado y el deseado.
    % W .* errorVec ponderación que hace que ciertos ángulos cuenten más
    % Se suma todos los valores y se hace la raíz cuadrada (calcular la
    % norma)
    cost_main = sqrt(sum((W .* errorVec).^2));

    % Penalización: para |phi| >= intermediateAngle, forzamos que no supere -5 dB.
    mask_threshold_dB = -5; % umbral en dB
    penalty_factor = 0.01; % factor de penalización ---------------------------COMO SE CUAL ES EL ADECUADO
    penalty = 0;
    pattern_synth_dB = mag2db(pattern_synth);
    for i = 1:length(azimuth)
        phi = azimuth(i);
        if abs(phi) >= intermediateAngle
            if pattern_synth_dB(i) > mask_threshold_dB
                overshoot = pattern_synth_dB(i) - mask_threshold_dB;
                penalty = penalty + penalty_factor * overshoot^2;
            end
        end
    end

    cost = cost_main + penalty;
end
