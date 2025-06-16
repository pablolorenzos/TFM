function cost = CostGA2ANCMask(x, N1, N2, azimuth, Beam_d, lmbda, intermediateAngle)
    x = x(:); %-------------------------- FORZAMOS X A COLUMNA
    N_total = N1 + N2;

    r1 = x(1);
    r2 = x(2);
    desfaseRel = x(3);

    % ---------------------FORZAMOS A Q SEAN VECTORES COLUMNA
    w_re = x(4 : 3+N_total).';
    w_im = x(4+N_total : 3+2*N_total).';
    w = w_re + 1i*w_im;

    % ángulos de cada anillo:
    anglesRing1 = (0:N1-1)*(360/N1);
    anglesRing2 = (0:N2-1)*(360/N2) + desfaseRel;

    % Calcular posiciones (x,y)
    posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
    posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
    elementPos = [posRing1, posRing2];

    % matriz de apuntamiento
    stvmat = zeros(N_total, length(azimuth));
    for i = 1:length(azimuth)
        phi = azimuth(i);
        stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda );
    end

    % ------------------------- FORZAMOS A Q SEA COLUMNA
    w = w(:);
    pattern_synth = abs(w' * stvmat);


% ----------------------------- AQUÍ INNTRODUCIMOS PENNALIZACIÓN
    % función de coste
    cost_main = norm(pattern_synth - Beam_d, 2);

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