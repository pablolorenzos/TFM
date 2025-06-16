%%
clear; clc; close all;

%% 0. Frecuencia y longitud de onda
f0 = 1.3e9;
lmbda = physconst('lightspeed')/f0;
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
Beam_d = Beam_d_lin;  % Patrón deseado

figure;
plot(azimuth, mag2db(Beam_d), 'b','LineWidth',2);
xlabel('\theta (°)'); ylabel('Ganancia (dB)');
title('Patrón deseado (half-cosine)');
grid on; ylim([-10 20]);

%% 2. Array circular de 16 elementos en dos anillos (8 + 8)
N1 = 8;  % Elementos en anillo 1
N2 = 8;  % Elementos en anillo 2
N_total = N1 + N2;

%% 3. Función objetivo
objfun = @(x) arrayCircularCost(x, N1, N2, azimuth, Beam_d, lmbda, intermediateAngle);

%% 4. Vector inicial y límites
x0 = zeros(3 + 2*N_total,1);
x0(1) = 0.5*lmbda;
x0(2) = 1.0*lmbda;
x0(3) = 0;
x0(4 : 3+N_total) = 1;
x0(4+N_total : 3+2*N_total) = 0;

% Límites
lb = zeros(size(x0));
ub = zeros(size(x0));

% r1
lb(1) = 0.3*lmbda;
ub(1) = 0.7*lmbda;
% r2
lb(2) = 0.8*lmbda;
ub(2) = 1.5*lmbda;
% desfase relativo
lb(3) = 0;
ub(3) = 360;
% pesos reales
lb(4 : 3+N_total) = -10;
ub(4 : 3+N_total) = 10;
% pesos imaginarios
lb(4+N_total : 3+2*N_total) = -10;  
ub(4+N_total : 3+2*N_total) = 10;

%% 5. Ejecutar optimización con fmincon (o patternsearch, según prefieras)
options = optimoptions('fmincon','Display','iter','MaxIterations',2000,MaxFunctionEvaluations=100000,OptimalityTolerance=1e-6,ConstraintTolerance=1e-6);
x_opt = fmincon(objfun, x0, [], [], [], [], lb, ub, [], options);

% Extraer variables optimizadas
r1_opt       = x_opt(1);
r2_opt       = x_opt(2);
desfaseRel_opt = x_opt(3);  % Este es el desfase relativo anillo 2
w_re_opt     = x_opt(4 : 3+N_total);
w_im_opt     = x_opt(4+N_total : 3+2*N_total);
w_opt        = w_re_opt + 1i*w_im_opt;

% en módulo y fase
mod_w = abs(w_opt);
phase_w = rad2deg(angle(w_opt));

%% 6. Reconstruir la geometría optimizada y calcular el patrón final
anglesRing1_opt = (0:N1-1)*(360/N1);
anglesRing2_opt = (0:N2-1)*(360/N2) + desfaseRel_opt;

posRing1_opt = [r1_opt*cosd(anglesRing1_opt); r1_opt*sind(anglesRing1_opt)];
posRing2_opt = [r2_opt*cosd(anglesRing2_opt); r2_opt*sind(anglesRing2_opt)];
elementPos_opt = [posRing1_opt, posRing2_opt];

stvmat_opt = zeros(N_total, length(azimuth));
for i = 1:length(azimuth)
    phi = azimuth(i);
    stvmat_opt(:, i) = exp(1i * 2*pi * (elementPos_opt(1,:)'*cosd(phi) + elementPos_opt(2,:)'*sind(phi)) / lmbda );
end
pattern_synth_opt = abs(w_opt'*stvmat_opt);

% Resultado final
figure;
plot(azimuth, mag2db([Beam_d; pattern_synth_opt])','LineWidth',2);
legend('Desired','Synthesized','Location','Best');
xlabel('\theta (°)'); ylabel('Ganancia (dB)');
title('Síntesis con radios y desfase relativo optimizado');
grid on; ylim([-10 20]);

%% 7. Gráfico con las posiciones optimizadas (estilo constelación)
figure;
scatter(elementPos_opt(1,:), elementPos_opt(2,:), 60, 'filled', 'b');
hold on;
thetaC = linspace(0,2*pi,200);
plot(r1_opt*cos(thetaC), r1_opt*sin(thetaC), 'k--');  
plot(r2_opt*cos(thetaC), r2_opt*sin(thetaC), 'k--');  
axis equal; grid on;
xlabel('Eje X (m)'); ylabel('Eje Y (m)');
title('Posiciones optimizadas de los 16 elementos (dos anillos)');

