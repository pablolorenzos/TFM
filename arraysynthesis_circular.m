clear; clc; close all;

%% 0. Frecuencia y longitud de onda
f0 = 1.3e9;                             % Frecuencia de operación: 1,3 GHz
lmbda = physconst('lightspeed')/f0;     % Longitud de onda
fprintf('λ = %.4f m\n', lmbda);

%% 1. Definir el patrón deseado (en azimut)
azimuth = -90:90;

% Niveles en dB para el patrón
peak_dB         = 15;  
valley_dB       = peak_dB - 20*log10(2121/1000);  
intermediate_dB = -5; 
outer_dB        = -5;  

% Ángulos clave
peakAngle        = 20;  
intermediateAngle= 40;  
outerAngle       = 80;  

Beam_d_lin = zeros(size(azimuth));
for k = 1:length(azimuth)
    th = azimuth(k);
    x  = abs(th);
    if x <= peakAngle
        xNorm = x / peakAngle;
        dB_val = valley_dB + (peak_dB - valley_dB)*(0.5 - 0.5*cos(pi*xNorm));
    elseif x <= intermediateAngle
        xNorm = (x - peakAngle)/(intermediateAngle - peakAngle);
        dB_val = peak_dB + (intermediate_dB - peak_dB)*(0.5 - 0.5*cos(pi*xNorm));
    elseif x <= outerAngle
        xNorm = (x - intermediateAngle)/(outerAngle - intermediateAngle);
        dB_val = intermediate_dB + (outer_dB - intermediate_dB)*(0.5 - 0.5*cos(pi*xNorm));
    else
        dB_val = outer_dB;
    end
    Beam_d_lin(k) = 10^(dB_val/20);
end
Beam_d = Beam_d_lin;

% Graficar patrón deseado
figure;
plot(azimuth, mag2db(Beam_d), 'b','LineWidth',2);
xlabel('\theta (°)'); ylabel('Ganancia (dB)');
title('Patrón deseado con transiciones half-cosine');
grid on; ylim([-10 20]);

%% 2. Array circular de 16 elementos en dos anillos (8 + 8)
% Radios de los dos anillos
r1 = 0.5*lmbda;  % anillo interior
r2 = 1*lmbda;  % anillo exterior

N1 = 8; % elementos anillo 1
N2 = 8; % elementos anillo 2
desfase1 = 0; % desfase anillo 1
desfase2 = 30; % desfase anillo 2

% Ángulos de los 8 elementos en cada anillo
anglesRing1 = (0:N1-1)*fix(360/N1) + desfase1;  % 0°, 45°, 90°, ..., 315°
anglesRing2 = (0:N2-1)*fix(360/N2) + desfase2;

% Posiciones (x,y) de cada anillo
posRing1 = [r1*cosd(anglesRing1); r1*sind(anglesRing1)];  % 2x8
posRing2 = [r2*cosd(anglesRing2); r2*sind(anglesRing2)];  % 2x8

% Concatenamos en una sola matriz 2x16
elementPos = [posRing1, posRing2];

%% 3. Construir la matriz de apuntamiento (steering vector)
N_total = size(elementPos,2);   % 16 elementos
numAngles = length(azimuth);
stvmat = zeros(N_total, numAngles);

for i = 1:numAngles
    phi = azimuth(i);
    % Elevación=0 => [cosd(phi); sind(phi)]
    % Fase = 2π/λ*(x*cosφ + y*sinφ)
    stvmat(:, i) = exp(1i * 2*pi * ...
        (elementPos(1,:)'*cosd(phi) + elementPos(2,:)'*sind(phi))/lmbda );
end

%% 4. Optimización de pesos complejos
w_i_re = ones(N_total,1);
w_i_im = zeros(N_total,1);
x_ini = [w_i_re; w_i_im];

lb = -ones(2*N_total,1);
ub =  ones(2*N_total,1);

objfun = @(x) norm( abs((x(1:N_total) + 1i*x(N_total+1:end)).'*stvmat) - Beam_d );

options = optimoptions('fmincon','Display','iter','MaxIterations',600,MaxFunctionEvaluations=80000,OptimalityTolerance=1e-6,ConstraintTolerance=1e-6);
x_opt = fmincon(objfun, x_ini, [], [], [], [], lb, ub, [], options);

% Pesos optimizados
w_opt = x_opt(1:N_total) + 1i*x_opt(N_total+1:end);
pattern_synth = abs(w_opt'*stvmat);

%% 5. Graficar comparación (dB)
figure;
plot(azimuth, mag2db([Beam_d; pattern_synth])','LineWidth',2);
legend('Desired','Synthesized','Location','Best');
xlabel('\theta (°)'); ylabel('Ganancia (dB)');
title('Síntesis con array circular');
grid on; ylim([-10 20]);

%% 6. Gráfico con las posiciones de los elementos (estilo constelación)
figure;
scatter(elementPos(1,:), elementPos(2,:), 60, 'filled', 'b');
hold on;
% Opcional: dibujar los círculos correspondientes a r1 y r2
thetaC = linspace(0,2*pi,200);
plot(r1*cos(thetaC), r1*sin(thetaC), 'k--');  % anillo interior
plot(r2*cos(thetaC), r2*sin(thetaC), 'k--');  % anillo exterior

axis equal; grid on;
xlabel('Eje X (m)');
ylabel('Eje Y (m)');
title('Posiciones de los elementos en el array circular');
