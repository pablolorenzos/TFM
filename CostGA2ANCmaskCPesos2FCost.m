function cost = CostGA2ANCmaskCPesos2FCost(x, N1, N2, azimuth, Beam_d, lmbda_eff, intermediateAngle)
    % Forzar x a ser vector columna
    x = x(:);
    N_total = N1 + N2;

    r1 = x(1);
    r2 = x(2);
    desfaseRel = x(3);

    w_re = x(4 : 3+N_total).';
    w_im = x(4+N_total : 3+2*N_total).';
    w = w_re + 1i * w_im;
    w = w(:);


    anglesRing1 = (0:N1-1) * (360/N1);
    anglesRing2 = (0:N2-1) * (360/N2) + desfaseRel;

    posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
    posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
    elementPos = [posRing1, posRing2];


    numAngles = length(azimuth);
    stvmat = zeros(N_total, numAngles);
    for i = 1:numAngles
        phi = azimuth(i);
        stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)' * cosd(phi) + elementPos(2,:)' * sind(phi) ) / lmbda_eff );
    end


    pattern_synth = abs(w' * stvmat);
    pattern_synth_dB = 20 * log10(pattern_synth);

    central_idx = find(abs(azimuth) <= intermediateAngle);
    external_idx = find(abs(azimuth) > intermediateAngle);


    tol = 0.1;  % Tolerancia en grados
    % Inicialmente, asignamos un peso base 1
    W = ones(size(azimuth));

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



    % Extraer el vector de ponderación para la región central
    W_central = W(central_idx);
    W_external = W(external_idx);

    % Calcular el error en la región central:
    error_central = pattern_synth(central_idx) - Beam_d(central_idx);
    cost_central = sqrt(sum((W_central .* error_central).^2));

    error_external = pattern_synth(external_idx) - Beam_d(external_idx);
    cost_external = sqrt(sum((W_external .* error_external).^2));

    mask_threshold_dB = -5;
    penalty_factor = 0.01;  % Ajusta este valor según lo deseado
    % penalty = 0;
    cost_external=0;
    for i = 1:length(external_idx)
        idx = external_idx(i);
        if pattern_synth_dB(idx) > mask_threshold_dB
            overshoot = pattern_synth_dB(idx) - mask_threshold_dB;
            cost_external = cost_external + penalty_factor * overshoot^2;
        end
    end
    cost = 2.5*cost_central + cost_external;
    fprintf('El valor central de coste es: %.2f\n', cost_central);
    fprintf('El valor externo de coste es: %.2f\n', cost_external);
end        