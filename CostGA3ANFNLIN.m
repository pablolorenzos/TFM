% AQUÍ VAMOS A INTRODUCIR UN TERCER ANILLO
% CON FUNCIÓN NO LINEAL
% SIN ELEMENTO CENTRAL

function [cost, cost_central, cost_external_1] = CostGA3ANFNLIN(x, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle)
    x = x(:); % FORZAMOS X A COLUMNA

    N_total = N1 + N2 + N3;

    r1 = x(1);
    r2 = x(2);
    r3 = x(3);
    desfaseRel = x(4);
    desfaseRel1 = x(5);

    % FORZAMOS A Q SEAN VECTORES COLUMNA
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

    % FORZAMOS A Q SEA COLUMNA
    w = w(:);
    pattern_synth = abs(w' * stvmat);
    pattern_synth_dB = 20 * log10(pattern_synth);


    %%
    central_idx = find(abs(azimuth) <= intermediateAngle);
    external_idx = find(abs(azimuth) > intermediateAngle);

    % IMPORTANCIA A ANGULOS CENTRALES // LO Q ESTA ENTRE %% ES QUITAR Y
    % PONER
        % Aquí definimos un vector de pesos para dar mayor importancia a ciertos ángulos.
    % Por ejemplo, damos más peso a los puntos en 0° (valle) y en ±peakAngle (picos).
    W = ones(size(azimuth)); % peso base 1 para todos los ángulos
    tol = 0.5;  % Tolerancia en grados para identificar el "punto exacto"
    % Definir valores de peso mayores en los puntos críticos:
    % w_valley = 20; % por ejemplo, darle tres veces más importancia al valle en 0°
    % w_peak = 20;   % y tres veces más importancia a los picos en ±peakAngle
    % 
    % w_a = 13;
    % w_b = 13;
    % w_c = 13;
    % w_d = 10;
    % w_e = 10;
    % w_f = 10;
    % w_g = 10;
    % w_h = 10;
    % w_i = 13;

    w_valley = 1;
    w_peak = 1;

    w_a = 1;
    w_b = 1;
    w_c = 1;
    w_d = 1;
    w_e = 1;
    w_f = 1;
    w_g = 1;
    w_h = 1;
    w_i = 1;


    % Para el valle (ángulo 0°)
    W(abs(azimuth) < tol) = w_valley;
    % Para los picos (ángulos cercanos a +peakAngle y -peakAngle)
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

    error_external = pattern_synth(external_idx) - Beam_d(external_idx);
    cost_external_1 = sqrt(sum((W_external .* error_external).^2));

    cost = cost_central + 2*cost_external_1;
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
% CON FUNCIÓN NO LINEAL
% CON ELEMENTO CENTRAL

% function [cost, cost_central, cost_external_1] = CostGA3ANFNLIN(x, N0, N1, N2, N3, azimuth, Beam_d, lmbda_eff, intermediateAngle)
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
%     %  FORZAMOS A Q SEA COLUMNA
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
% 
%     w_valley = 1;
%     w_peak = 1;
% 
%     w_a = 1;
%     w_b = 1;
%     w_c = 1;
%     w_d = 1;
%     w_e = 1;
%     w_f = 1;
%     w_g = 1;
%     w_h = 1;
%     w_i = 1;
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
%     error_external = pattern_synth(external_idx) - Beam_d(external_idx);
%     cost_external_1 = sqrt(sum((W_external .* error_external).^2));
% 
%     cost = cost_central + 2*cost_external_1;
%     cost_central  = cost_central;
%     cost_external_1 = cost_external_1;
% 
% end