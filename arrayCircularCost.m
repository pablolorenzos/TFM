function cost = arrayCircularCost(x, N1, N2, azimuth, Beam_d, lmbda, intermediateAngle)

N_total = N1 + N2;

r1 = x(1);
r2 = x(2);
desfaseRel = x(3);  % desfase relativo entre anillo 2 y anillo 1

w_re = x(4 : 3+N_total);
w_im = x(4+N_total : 3+2*N_total);
w    = w_re + 1i*w_im;

anglesRing1 = (0:N1-1) * (360/N1);           
anglesRing2 = (0:N2-1) * (360/N2) + desfaseRel;

% Calcular las posiciones (x,y)
posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];  % 2xN1
posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];  % 2xN2
elementPos = [posRing1, posRing2];  % 2 x (N1+N2)

% matriz de apuntamiento
stvmat = zeros(N_total, length(azimuth));
for i = 1:length(azimuth)
    phi = azimuth(i);
    stvmat(:, i) = exp(1i * 2*pi * (elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi)) / lmbda );
end

% Patrón sintetizado
pattern_synth = abs(w'*stvmat);

%% PENALIZACION
cost_main = norm(pattern_synth - Beam_d, 2);

% pattern_synth a dB
pattern_synth_dB = mag2db(pattern_synth);

mask_threshold_dB = -5;
penalty_factor = 0.5;

penalty = 0;


for i = 1:length(azimuth)
    phi = azimuth(i);
    % Penalizamos solo si |phi| >= intermediateAngle
    if abs(phi) >= intermediateAngle
        % Si el patrón supera -5 dB, se añade penalización
        if pattern_synth_dB(i) > mask_threshold_dB
            overshoot = pattern_synth_dB(i) - mask_threshold_dB;
            penalty = penalty + penalty_factor * overshoot^2;
        end
    end
end

 % Costo total
 cost = cost_main + penalty;
end
