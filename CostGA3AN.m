%%
% AQUÍ VAMOS A INTRODUCIR UN TERCER ANILLO
% SIN ELEMENTO CENTRAL

function [cost, cost_central, cost_external_1] = CostGA3AN(x, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle)
    x = x(:); %-------------------------- FORZAMOS X A COLUMNA

    N_total = N1 + N2 + N3;

    r1 = x(1);
    r2 = x(2);
    r3 = x(3);
    desfaseRel = x(4);
    desfaseRel1 = x(5);

    % ---------------------FORZAMOS A Q SEAN VECTORES COLUMNA
    w_re = x(6 : 5+N_total).';
    w_im = x(6+N_total : 5+2*N_total).';
    w = w_re + 1i*w_im;

    %   Elemento central

    % ángulos de cada anillo:
    anglesRing1 = (0:N1-1)*(360/N1);
    anglesRing2 = (0:N2-1)*(360/N2) + desfaseRel;
    anglesRing3 = (0:N3-1)*(360/N3) + desfaseRel1;

    % Calcular posiciones (x,y)
    posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
    posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
    posRing3 = [r3*cosd(anglesRing3); r3*sind(anglesRing3)];
    elementPos = [posRing1, posRing2, posRing3];

    % matriz de apuntamiento
    stvmat = zeros(N_total, length(azimuth));
    for i = 1:length(azimuth)
        phi = azimuth(i);
        stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda_eff );
    end

    % ------------------------- FORZAMOS A Q SEA COLUMNA
    w = w(:);
    pattern_synth = abs(w' * stvmat);
    pattern_synth_dB = 20 * log10(pattern_synth);

    %%
    central_idx = find(abs(azimuth) <= intermediateAngle);
    external_idx = find(abs(azimuth) > intermediateAngle);

    % IMPORTANCIA A ANGULOS CENTRALES 
        % Aquí definimos un vector de pesos para dar mayor importancia a ciertos ángulos.
    % Por ejemplo, damos más peso a los puntos en 0° (valle) y en ±peakAngle (picos).
    W = ones(size(azimuth)); % peso base 1 para todos los ángulos
    tol = 1;  % Tolerancia en grados para identificar el "punto exacto"
    % Definir valores de peso mayores en los puntos críticos:
    % w_valley = 20; % por ejemplo, darle tres veces más importancia al valle en 0°
    % w_peak = 20;   % y tres veces más importancia a los picos en ±peakAngle

    w_a = 13;
    w_b = 13;
    w_c = 13;
    w_d = 10;
    w_e = 10;
    w_f = 10;
    w_g = 10;
    w_h = 10;
    w_i = 13;


    Para el valle (ángulo 0°)
    W(abs(azimuth) < tol) = w_valley;
    Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
    W(abs(azimuth - 54) < tol) = w_peak;
    W(abs(azimuth + 54) < tol) = w_peak;

    W(abs(azimuth - 5) < tol) = w_a;
    W(abs(azimuth + 5) < tol) = w_a;

    W(abs(azimuth - 10) < tol) = w_b;
    W(abs(azimuth + 10) < tol) = w_b;

    W(abs(azimuth - 15) < tol) = w_c;
    W(abs(azimuth + 15) < tol) = w_c;

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

    W(abs(azimuth - 60) < tol) = w_i;
    W(abs(azimuth + 60) < tol) = w_i;


    % Extraer el vector de ponderación para la región central
    W_central = W(central_idx);
    W_external = W(external_idx);

    % Calcular el error en la región central:
    error_central = pattern_synth(central_idx) - Beam_d(central_idx);
    cost_central = sqrt(sum((W_central .* error_central).^2));

    % Penalización: para |phi| >= intermediateAngle, forzamos que no supere -5 dB.
    mask_threshold_dB = 15 - 20; % umbral en dB
    penalty_factor = 0.2; % factor de penalización
    % penalty = 0;
    cost_external_1 = 0;

    for i = 1:length(external_idx)
        idx = external_idx(i);
        if pattern_synth_dB(idx) > mask_threshold_dB
            overshoot = pattern_synth_dB(idx) - mask_threshold_dB;
            cost_external_1 = cost_external_1 + penalty_factor * overshoot^2;
        end
    end

    cost = 2*cost_central + cost_external_1;

    cost_central  = cost_central;
    cost_external_1 = cost_external_1;
end

%%
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
% AQUÍ VAMOS A INTRODUCIR UN TERCER ANILLO
% CON ELEMENTO CENTRAL

% function [cost, cost_central, cost_external_1] = CostGA3AN(x, N0, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle)
%     x = x(:); % FORZAMOS X A COLUMNA
%     N_total = N0 + N1 + N2 + N3;
% 
%     r1 = x(1);
%     r2 = x(2);
%     r3 = x(3);
%     desfaseRel = x(4);
%     desfaseRel1 = x(5);
% 
%     % FORZAMOS A Q SEAN VECTORES COLUMNA
%     w_re = x(6 : 5+N_total).';
%     w_im = x(6+N_total : 5+2*N_total).';
%     w = w_re + 1i*w_im;
% 
%       % Elemento central
%     pos0 = [0;0];
% 
%     % ángulos de cada anillo:
%     anglesRing1 = (0:N1-1)*(360/N1);
%     anglesRing2 = (0:N2-1)*(360/N2) + desfaseRel;
%     anglesRing3 = (0:N3-1)*(360/N3) + desfaseRel1;
% 
%     % Calcular posiciones (x,y)
%     posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
%     posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
%     posRing3 = [r3*cosd(anglesRing3); r3*sind(anglesRing3)];
%     elementPos = [pos0, posRing1, posRing2, posRing3];
% 
%     % matriz de apuntamiento
%     stvmat = zeros(N_total, length(azimuth));
%     for i = 1:length(azimuth)
%         phi = azimuth(i);
%         stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda_eff );
%     end
% 
%     % FORZAMOS A Q SEA COLUMNA
%     w = w(:);
%     pattern_synth = abs(w' * stvmat);
%     pattern_synth_dB = 20 * log10(pattern_synth);
% 
% 
%     central_idx = find(abs(azimuth) <= intermediateAngle);
%     external_idx = find(abs(azimuth) > intermediateAngle);
% 
% 
%     W = ones(size(azimuth)); % peso base 1 para todos los ángulos
%     tol = 0.5;  % Tolerancia en grados para identificar el "punto exacto"
% 
%     w_valley = 20; % por ejemplo, darle tres veces más importancia al valle en 0°
%     w_peak = 20;   % y tres veces más importancia a los picos en ±peakAngle
% 
%     w_a = 13;
%     w_b = 13;
%     w_c = 13;
%     w_d = 10;
%     w_e = 10;
%     w_f = 10;
%     w_g = 10;
%     w_h = 10;
%     w_i = 13;
% 
%     % Para el valle (ángulo 0°)
%     W(abs(azimuth) < tol) = w_valley;
%     % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
%     W(abs(azimuth - 54) < tol) = w_peak;
%     W(abs(azimuth + 54) < tol) = w_peak;
% 
%     W(abs(azimuth - 5) < tol) = w_a;
%     W(abs(azimuth + 5) < tol) = w_a;
% 
%     W(abs(azimuth - 10) < tol) = w_b;
%     W(abs(azimuth + 10) < tol) = w_b;
% 
%     W(abs(azimuth - 15) < tol) = w_c;
%     W(abs(azimuth + 15) < tol) = w_c;
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
%     W(abs(azimuth - 60) < tol) = w_i;
%     W(abs(azimuth + 60) < tol) = w_i;
% 
% 
%     W_central = W(central_idx);
%     W_external = W(external_idx);
% 
%     error_central = pattern_synth(central_idx) - Beam_d(central_idx);
%     cost_central = sqrt(sum((W_central .* error_central).^2));
% 
%    % error_external = pattern_synth(external_idx) - Beam_d(external_idx);
%    % cost_external = sqrt(sum((W_external .* error_external).^2));
% 
%     % Penalización: para |phi| >= intermediateAngle, forzamos que no supere -5 dB.
%     mask_threshold_dB = -5; % umbral en dB
%     penalty_factor = 0.1; % factor de penalización
%     % penalty = 0;
%     cost_external_1 = 0;
% 
%     for i = 1:length(external_idx)
%         idx = external_idx(i);
%         if pattern_synth_dB(idx) > mask_threshold_dB
%             overshoot = pattern_synth_dB(idx) - mask_threshold_dB;
%             cost_external_1 = cost_external_1 + penalty_factor * overshoot^2;
%             % penalty = penalty + penalty_factor * overshoot^2;
%         end
%     end
%     % cost_external_1 = cost_external + sqrt(penalty);
% 
%     cost = 2.5*cost_central + cost_external_1;
%     cost_central  = cost_central;
%     cost_external_1 = cost_external_1;
% 
% end

%%
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------
% OPTIMIZAMOS NUMERO DE ELEMENTOS CON ELEMENTO CENTRAL


% function [cost, cost_central, cost_external_1] = CostGA3AN(x, azimuth, Beam_d, lmbda_eff, intermediateAngle)
%     x = x(:); % FORZAMOS X A COLUMNA
% 
% 
%     r1 = x(1);
%     r2= x(2);
%     r3 = x(3);
%     desfaseRel = x(4);
%     desfaseRel1 = x(5);
%     N0 = round(x(6));
%     N1 = round(x(7));
%     N2 = round(x(8));
%     N3 = round(x(9));
% 
%         % Ensure N1 and N2 are integers within bounds
%     N0 = max(0, min(1, N0)); % Limit N1 between 4 and 16
%     N1 = max(2, min(4, N1)); % Limit N1 between 4 and 16
%     N2 = max(4, min(8, N2)); % Limit N2 between 4 and 16
%     N3 = max(4, min(10, N3)); % Limit N2 between 4 and 16
% 
%     % Total number of elements
%     N_total = N0 + N1 + N2 + N3;
% 
%     % FORZAMOS A Q SEAN VECTORES COLUMNA
%     w_re = x(10 : 9+N_total).';
%     w_im = x(10+N_total : 9+2*N_total).';
%     w = w_re + 1i*w_im;
%     % Ensure weights are column vectors
%     w = w(:);
% 
%     % Central position (if N0 = 1)
%     if N0 > 0
%         pos0 = [0; 0];
%     else
%         pos0 = [];
%     end
% 
%     % ángulos de cada anillo:
%     anglesRing1 = (0:N1-1)*(360/N1);
%     anglesRing2 = (0:N2-1)*(360/N2) + desfaseRel;
%     anglesRing3 = (0:N3-1)*(360/N3) + desfaseRel1;
% 
%     % Calcular posiciones (x,y)
%     posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
%     posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
%     posRing3 = [r3*cosd(anglesRing3); r3*sind(anglesRing3)];
%     elementPos = [pos0, posRing1, posRing2, posRing3];
% 
%     % matriz de apuntamiento
%     stvmat = zeros(N_total, length(azimuth));
%     for i = 1:length(azimuth)
%         phi = azimuth(i);
%         stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda_eff );
%     end
% 
%     % FORZAMOS A Q SEA COLUMNA
%     w = w(:);
%     pattern_synth = abs(w' * stvmat);
%     pattern_synth_dB = 20 * log10(pattern_synth);
% 
% 
%     %%
%     central_idx = find(abs(azimuth) <= intermediateAngle);
%     external_idx = find(abs(azimuth) > intermediateAngle);
% 
%     W = ones(size(azimuth)); % peso base 1 para todos los ángulos
%     tol = 0.5;  % Tolerancia en grados para identificar el "punto exacto"
%     % Definir valores de peso mayores en los puntos críticos:
%     w_valley = 20; % por ejemplo, darle tres veces más importancia al valle en 0°
%     w_peak = 20;   % y tres veces más importancia a los picos en ±peakAngle
% 
%     w_a = 13;
%     w_b = 13;
%     w_c = 13;
%     w_d = 10;
%     w_e = 10;
%     w_f = 10;
%     w_g = 10;
%     w_h = 10;
%     w_i = 13;
% 
%     % Para el valle (ángulo 0°)
%     W(abs(azimuth) < tol) = w_valley;
%     % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
%     W(abs(azimuth - 54) < tol) = w_peak;
%     W(abs(azimuth + 54) < tol) = w_peak;
% 
%     W(abs(azimuth - 5) < tol) = w_a;
%     W(abs(azimuth + 5) < tol) = w_a;
% 
%     W(abs(azimuth - 10) < tol) = w_b;
%     W(abs(azimuth + 10) < tol) = w_b;
% 
%     W(abs(azimuth - 15) < tol) = w_c;
%     W(abs(azimuth + 15) < tol) = w_c;
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
%     W(abs(azimuth - 60) < tol) = w_i;
%     W(abs(azimuth + 60) < tol) = w_i;
% 
%     % Extraer el vector de ponderación para la región central
%     W_central = W(central_idx);
%     W_external = W(external_idx);
% 
%     % Calcular el error en la región central:
%     error_central = pattern_synth(central_idx) - Beam_d(central_idx);
%     cost_central = sqrt(sum((W_central .* error_central).^2));
% 
%     % error_external = pattern_synth(external_idx) - Beam_d(external_idx);
%     % cost_external = sqrt(sum((W_external .* error_external).^2));
% 
% 
%     % Penalización: para |phi| >= intermediateAngle, forzamos que no supere -5 dB.
%     mask_threshold_dB = -5; % umbral en dB
%     penalty_factor = 0.1; % factor de penalización
%     % penalty = 0;
%     cost_external_1 = 0;
% 
%     for i = 1:length(external_idx)
%         idx = external_idx(i);
%         if pattern_synth_dB(idx) > mask_threshold_dB
%             overshoot = pattern_synth_dB(idx) - mask_threshold_dB;
%             cost_external_1 = cost_external_1 + penalty_factor * overshoot^2;
%             % penalty = penalty + penalty_factor * overshoot^2;
%         end
%     end
%     % cost_external_1 = cost_external + sqrt(penalty);
% 
%     cost = 2*cost_central + cost_external_1;
%     cost_central  = cost_central;
%     cost_external_1 = cost_external_1;
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

% function [cost, cost_central, cost_external_1] = CostGA3AN(x, azimuth, Beam_d, lmbda_eff, intermediateAngle)
%     x = x(:); % FORZAMOS X A COLUMNA
% 
% 
%     r1 = x(1);
%     r2= x(2);
%     r3 = x(3);
%     desfaseRel = x(4);
%     desfaseRel1 = x(5);
%     N1 = round(x(6));
%     N2 = round(x(7));
%     N3 = round(x(8));
% 
%         % Ensure N1 and N2 are integers within bounds
%     N1 = max(2, min(4, N1)); % Limit N1 between 4 and 16
%     N2 = max(4, min(8, N2)); % Limit N2 between 4 and 16
%     N3 = max(4, min(10, N3)); % Limit N2 between 4 and 16
% 
%     % Total number of elements
%     N_total = N1 + N2 + N3;
% 
%     % FORZAMOS A Q SEAN VECTORES COLUMNA
%     w_re = x(9 : 8+N_total).';
%     w_im = x(9+N_total : 8+2*N_total).';
%     w = w_re + 1i*w_im;
%     % Ensure weights are column vectors
%     w = w(:);
% 
%     % ángulos de cada anillo:
%     anglesRing1 = (0:N1-1)*(360/N1);
%     anglesRing2 = (0:N2-1)*(360/N2) + desfaseRel;
%     anglesRing3 = (0:N3-1)*(360/N3) + desfaseRel1;
% 
%     % Calcular posiciones (x,y)
%     posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
%     posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
%     posRing3 = [r3*cosd(anglesRing3); r3*sind(anglesRing3)];
%     elementPos = [posRing1, posRing2, posRing3];
% 
%     % matriz de apuntamiento
%     stvmat = zeros(N_total, length(azimuth));
%     for i = 1:length(azimuth)
%         phi = azimuth(i);
%         stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda_eff );
%     end
% 
%     % FORZAMOS A Q SEA COLUMNA
%     w = w(:);
%     pattern_synth = abs(w' * stvmat);
%     pattern_synth_dB = 20 * log10(pattern_synth);
% 
% 
%     %%
%     central_idx = find(abs(azimuth) <= intermediateAngle);
%     external_idx = find(abs(azimuth) > intermediateAngle);
% 
%     W = ones(size(azimuth)); % peso base 1 para todos los ángulos
%     tol = 0.5;  % Tolerancia en grados para identificar el "punto exacto"
%     % Definir valores de peso mayores en los puntos críticos:
%     w_valley = 20; % por ejemplo, darle tres veces más importancia al valle en 0°
%     w_peak = 20;   % y tres veces más importancia a los picos en ±peakAngle
% 
%     w_a = 13;
%     w_b = 13;
%     w_c = 13;
%     w_d = 10;
%     w_e = 10;
%     w_f = 10;
%     w_g = 10;
%     w_h = 10;
%     w_i = 13;
% 
%     % Para el valle (ángulo 0°)
%     W(abs(azimuth) < tol) = w_valley;
%     % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
%     W(abs(azimuth - 54) < tol) = w_peak;
%     W(abs(azimuth + 54) < tol) = w_peak;
% 
%     W(abs(azimuth - 5) < tol) = w_a;
%     W(abs(azimuth + 5) < tol) = w_a;
% 
%     W(abs(azimuth - 10) < tol) = w_b;
%     W(abs(azimuth + 10) < tol) = w_b;
% 
%     W(abs(azimuth - 15) < tol) = w_c;
%     W(abs(azimuth + 15) < tol) = w_c;
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
%     W(abs(azimuth - 60) < tol) = w_i;
%     W(abs(azimuth + 60) < tol) = w_i;
% 
%     % Extraer el vector de ponderación para la región central
%     W_central = W(central_idx);
%     W_external = W(external_idx);
% 
%     % Calcular el error en la región central:
%     error_central = pattern_synth(central_idx) - Beam_d(central_idx);
%     cost_central = sqrt(sum((W_central .* error_central).^2));
% 
%     % error_external = pattern_synth(external_idx) - Beam_d(external_idx);
%     % cost_external = sqrt(sum((W_external .* error_external).^2));
% 
%     % Penalización: para |phi| >= intermediateAngle, forzamos que no supere -5 dB.
%     mask_threshold_dB = -5; % umbral en dB
%     penalty_factor = 0.1; % factor de penalización
%     % penalty = 0;
%     cost_external_1 = 0;
% 
%     for i = 1:length(external_idx)
%         idx = external_idx(i);
%         if pattern_synth_dB(idx) > mask_threshold_dB
%             overshoot = pattern_synth_dB(idx) - mask_threshold_dB;
%             cost_external_1 = cost_external_1 + penalty_factor * overshoot^2;
%             % penalty = penalty + penalty_factor * overshoot^2;
%         end
%     end
%     % cost_external_1 = cost_external + sqrt(penalty);
% 
%     cost = 2*cost_central + cost_external_1;
% 
%     cost_central  = cost_central;
%     cost_external_1 = cost_external_1;
% end



%%
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
%
% SIN ELEMENTO CENTRAL
% CUANTIZACIÓN FASES

% function [cost, cost_central, cost_external_1] = CostGA3AN(x, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle)
%     x = x(:); % FORZAMOS X A COLUMNA
% 
%     N_total = N1 + N2 + N3;
% 
%     r1 = x(1);
%     r2 = x(2);
%     r3 = x(3);
%     desfaseRel = x(4);
%     desfaseRel1 = x(5);
% 
%         % Magnitudes
%     mags = x(6 : 5+N_total);
%     idxPhase = x(6+N_total : 5+2*N_total);      % índices enteros 0…63
%     phi = 2*pi * (idxPhase/64);       % radianes
%     w   = mags .* exp(1i*phi);        % vector de pesos
% 
%     % ángulos de cada anillo:
%     anglesRing1 = (0:N1-1)*(360/N1);
%     anglesRing2 = (0:N2-1)*(360/N2) + desfaseRel;
%     anglesRing3 = (0:N3-1)*(360/N3) + desfaseRel1;
% 
%     % Calcular posiciones (x,y)
%     posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
%     posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
%     posRing3 = [r3*cosd(anglesRing3); r3*sind(anglesRing3)];
%     elementPos = [posRing1, posRing2, posRing3];
% 
%     % matriz de apuntamiento
%     stvmat = zeros(N_total, length(azimuth));
%     for i = 1:length(azimuth)
%         phi = azimuth(i);
%         stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda_eff );
%     end
% 
%     % FORZAMOS A Q SEA COLUMNA
%     w = w(:);
%     pattern_synth = abs(w' * stvmat);
%     pattern_synth_dB = 20 * log10(pattern_synth);
% 
% 
%     %%
%     central_idx = find(abs(azimuth) <= intermediateAngle);
%     external_idx = find(abs(azimuth) > intermediateAngle);
% 
%     % IMPORTANCIA A ANGULOS CENTRALES // LO Q ESTA ENTRE %% ES QUITAR Y
%     % PONER
%         % Aquí definimos un vector de pesos para dar mayor importancia a ciertos ángulos.
%     % Por ejemplo, damos más peso a los puntos en 0° (valle) y en ±peakAngle (picos).
%     W = ones(size(azimuth)); % peso base 1 para todos los ángulos
%     tol = 0.5;  % Tolerancia en grados para identificar el "punto exacto"
%     % Definir valores de peso mayores en los puntos críticos:
%     w_valley = 20; % por ejemplo, darle tres veces más importancia al valle en 0°
%     w_peak = 20;   % y tres veces más importancia a los picos en ±peakAngle
% 
%     w_a = 13;
%     w_b = 13;
%     w_c = 13;
%     w_d = 10;
%     w_e = 10;
%     w_f = 10;
%     w_g = 10;
%     w_h = 10;
%     w_i = 13;
% 
%     % Para el valle (ángulo 0°)
%     W(abs(azimuth) < tol) = w_valley;
%     % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
%     W(abs(azimuth - 54) < tol) = w_peak;
%     W(abs(azimuth + 54) < tol) = w_peak;
% 
%     W(abs(azimuth - 5) < tol) = w_a;
%     W(abs(azimuth + 5) < tol) = w_a;
% 
%     W(abs(azimuth - 10) < tol) = w_b;
%     W(abs(azimuth + 10) < tol) = w_b;
% 
%     W(abs(azimuth - 15) < tol) = w_c;
%     W(abs(azimuth + 15) < tol) = w_c;
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
%     W(abs(azimuth - 60) < tol) = w_i;
%     W(abs(azimuth + 60) < tol) = w_i;
% 
%     % Extraer el vector de ponderación para la región central
%     W_central = W(central_idx);
%     W_external = W(external_idx);
% 
%     % Calcular el error en la región central:
%     error_central = pattern_synth(central_idx) - Beam_d(central_idx);
%     cost_central = sqrt(sum((W_central .* error_central).^2));
% 
% 
%     % Penalización: para |phi| >= intermediateAngle, forzamos que no supere -5 dB.
%     mask_threshold_dB = -5; % umbral en dB
%     penalty_factor = 0.1; % factor de penalización
%     cost_external_1 = 0;
% 
%     for i = 1:length(external_idx)
%         idx = external_idx(i);
%         if pattern_synth_dB(idx) > mask_threshold_dB
%             overshoot = pattern_synth_dB(idx) - mask_threshold_dB;
%             cost_external_1 = cost_external_1 + penalty_factor * overshoot^2;
%         end
%     end
%     cost = 2.5*cost_central + cost_external_1;
%     cost_central  = cost_central;
%     cost_external_1 = cost_external_1;
% 
% end


%%
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
% 

% CON ELEMENTO CENTRAL
% CUANTIZACIÓN FASES

% function [cost, cost_central, cost_external_1] = CostGA3AN(x, N0, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle)
%     x = x(:); % FORZAMOS X A COLUMNA
%     N_total = N0 + N1 + N2 + N3;
% 
%     r1 = x(1);
%     r2 = x(2);
%     r3 = x(3);
%     desfaseRel = x(4);
%     desfaseRel1 = x(5);
% 
%         % Magnitudes
%     mags = x(6 : 5+N_total);
%     idxPhase = x(6+N_total : 5+2*N_total);      % índices enteros 0…63
%     phi = 2*pi * (idxPhase/64);       % radianes
%     w   = mags .* exp(1i*phi);        % vector de pesos
% 
%       % Elemento central
%     pos0 = [0;0];
% 
%     % ángulos de cada anillo:
%     anglesRing1 = (0:N1-1)*(360/N1);
%     anglesRing2 = (0:N2-1)*(360/N2) + desfaseRel;
%     anglesRing3 = (0:N3-1)*(360/N3) + desfaseRel1;
% 
%     % Calcular posiciones (x,y)
%     posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
%     posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
%     posRing3 = [r3*cosd(anglesRing3); r3*sind(anglesRing3)];
%     elementPos = [pos0, posRing1, posRing2, posRing3];
% 
%     % matriz de apuntamiento
%     stvmat = zeros(N_total, length(azimuth));
%     for i = 1:length(azimuth)
%         phi = azimuth(i);
%         stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda_eff );
%     end
% 
%     % FORZAMOS A Q SEA COLUMNA
%     w = w(:);
%     pattern_synth = abs(w' * stvmat);
%     pattern_synth_dB = 20 * log10(pattern_synth);
% 
% 
%     central_idx = find(abs(azimuth) <= intermediateAngle);
%     external_idx = find(abs(azimuth) > intermediateAngle);
% 
% 
%     W = ones(size(azimuth)); % peso base 1 para todos los ángulos
%     tol = 0.5;  % Tolerancia en grados para identificar el "punto exacto"
% 
%     w_valley = 20; % por ejemplo, darle tres veces más importancia al valle en 0°
%     w_peak = 20;   % y tres veces más importancia a los picos en ±peakAngle
% 
%     w_a = 13;
%     w_b = 13;
%     w_c = 13;
%     w_d = 10;
%     w_e = 10;
%     w_f = 10;
%     w_g = 10;
%     w_h = 10;
%     w_i = 13;
% 
%     % Para el valle (ángulo 0°)
%     W(abs(azimuth) < tol) = w_valley;
%     % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
%     W(abs(azimuth - 54) < tol) = w_peak;
%     W(abs(azimuth + 54) < tol) = w_peak;
% 
%     W(abs(azimuth - 5) < tol) = w_a;
%     W(abs(azimuth + 5) < tol) = w_a;
% 
%     W(abs(azimuth - 10) < tol) = w_b;
%     W(abs(azimuth + 10) < tol) = w_b;
% 
%     W(abs(azimuth - 15) < tol) = w_c;
%     W(abs(azimuth + 15) < tol) = w_c;
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
%     W(abs(azimuth - 60) < tol) = w_i;
%     W(abs(azimuth + 60) < tol) = w_i;
% 
% 
%     W_central = W(central_idx);
%     W_external = W(external_idx);
% 
%     error_central = pattern_synth(central_idx) - Beam_d(central_idx);
%     cost_central = sqrt(sum((W_central .* error_central).^2));
% 
%    % error_external = pattern_synth(external_idx) - Beam_d(external_idx);
%    % cost_external = sqrt(sum((W_external .* error_external).^2));
% 
%     % Penalización: para |phi| >= intermediateAngle, forzamos que no supere -5 dB.
%     mask_threshold_dB = -5; % umbral en dB
%     penalty_factor = 0.1; % factor de penalización
%     % penalty = 0;
%     cost_external_1 = 0;
% 
%     for i = 1:length(external_idx)
%         idx = external_idx(i);
%         if pattern_synth_dB(idx) > mask_threshold_dB
%             overshoot = pattern_synth_dB(idx) - mask_threshold_dB;
%             cost_external_1 = cost_external_1 + penalty_factor * overshoot^2;
%         end
%     end
% 
%     cost = 2.5*cost_central + cost_external_1;
%     cost_central  = cost_central;
%     cost_external_1 = cost_external_1;
% 
% end
