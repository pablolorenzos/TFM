% AQUÍ VAMOS A INTRODUCIR UN TERCER ANILLO
% CON FUNCIÓN NO LINEAL
% SIN ELEMENTO CENTRAL

function [c, ceq] = my3ANFNLin(x, N1, N2, N3, azimuth, lmbda, intermediateAngle, mask_dB)

    N_total = N1 + N2 + N3;

    x = x(:);

    r1 = x(1);
    r2 = x(2);
    r3 = x(3);
    desfaseRel = x(4);
    desfaseRel1 = x(5);

    % FORZAMOS A Q SEAN VECTORES COLUMNA
    w_re = x(6 : 5+N_total).';
    w_im = x(6+N_total : 5+2*N_total).';
    w = w_re + 1i*w_im;

    % ángulos de cada anillo:
    anglesRing1 = (0:N1-1)*(360/N1);
    anglesRing2 = (0:N2-1)*(360/N2) + desfaseRel;
    anglesRing3 = (0:N3-1)*(360/N3) + desfaseRel1;

    % Calcular posiciones (x,y)
    posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];
    posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];
    posRing3 = [r3*cosd(anglesRing3); r3*sind(anglesRing3)];
    elementPos = [posRing1, posRing2, posRing3];


    stvmat = zeros(N_total, length(azimuth));
    for i = 1:length(azimuth)
        phi = azimuth(i);
        stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda );
    end


    w = w(:);  % vector columna
    pattern = abs(w' * stvmat);
    pattern_dB = 20*log10(pattern);

    %  para cada ángulo con |phi| >= intermediateAngle,
    % se requiere que pattern_dB(i) <= mask_dB. que es lo mismo que poner
    % c_i(x) = pattern_dB(i) - mask_dB <= 0.
    c = [];
    for i = 1:length(azimuth)
        phi = azimuth(i);
        if abs(phi) >= intermediateAngle
            c(end+1,1) = pattern_dB(i) - mask_dB;
        end
    end

    ceq = [];  % Se queda vacio porque no tenemos restricciones de igualdad
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

% function [c, ceq] = my3ANFNLin(x, N0, N1, N2, N3, azimuth, lmbda, intermediateAngle, mask_dB)
% 
%     N_total = N0 + N1 + N2 + N3;
% 
%     x = x(:);
% 
%     r1 = x(1);
%     r2 = x(2);
%     r3 = x(3);
%     desfaseRel = x(4);
%     desfaseRel1 = x(5);
% 
%       % Elemento central
%     pos0 = [0;0];
% 
%     % FORZAMOS A Q SEAN VECTORES COLUMNA
%     w_re = x(6 : 5+N_total).';
%     w_im = x(6+N_total : 5+2*N_total).';
%     w = w_re + 1i*w_im;
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
% 
%     stvmat = zeros(N_total, length(azimuth));
%     for i = 1:length(azimuth)
%         phi = azimuth(i);
%         stvmat(:, i) = exp(1i * 2*pi * ( elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi) ) / lmbda );
%     end
% 
% 
%     w = w(:);  % vector columna
%     pattern = abs(w' * stvmat);
%     pattern_dB = 20*log10(pattern);
% 
%     % para cada ángulo con |phi| >= intermediateAngle,
%     % se requiere que pattern_dB(i) <= mask_dB. que es lo mismo que poner
%     % c_i(x) = pattern_dB(i) - mask_dB <= 0.
%     c = [];
%     for i = 1:length(azimuth)
%         phi = azimuth(i);
%         if abs(phi) >= intermediateAngle
%             c(end+1,1) = pattern_dB(i) - mask_dB;
%         end
%     end
% 
%     ceq = [];  % Se queda vacio porque no tenemos restricciones de igualdad
% end