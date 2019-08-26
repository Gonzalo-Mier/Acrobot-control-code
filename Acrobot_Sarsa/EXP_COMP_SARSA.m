
% Carpeta = 'Im_S_Sutton_i_a02';

mkdir(Carpeta);
delete(strcat(Carpeta,'/*'));
rmdir(Carpeta);
mkdir(Carpeta);


% MAX_EXP = 2;


% alpha = 0.2/48;
% lambda = 0.9;
% gamma = 1;
% eps = 0;
% MaxIt=600;

% % % Parametros iniciales del problema
m1 = 1;
m2 = 1;
lc1 = 0.5;
lc2 = 0.5;
l1 = 1;
l2 = 1;
I1 = 1;
I2 = 1;
g = 9.8;
dt = 0.05;
% t_total = 1000;

% % % Inicio de parametros % % %
q1_ini = 0;
q2_ini = 0;
dq1_ini = 0;
dq2_ini = 0;
ddq1_ini = 0;
ddq2_ini = 0;
F_ini = 0;



experimentos_SARSA
experimentos_lambda_SARSA




figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3,1,2);
plot(1:MaxIt,sum(it_lS_exp)/MAX_EXP*4*dt,1:MaxIt,sum(it_S_exp)/MAX_EXP*4*dt);
title('Media de pasos por iteración en 10 iteraciones')
legend('Lambda-SARSA','SARSA')
xlabel('Iteración')
ylabel('Tiempo (s)')
grid;

subplot(3,1,1);
plot(1:MaxIt,max(it_lS_exp)*4*dt,1:MaxIt,max(it_S_exp)*4*dt);
title('Máximo de pasos por iteración en 10 iteraciones')
legend('Lambda-SARSA','SARSA')
xlabel('Iteración')
ylabel('Tiempo (s)')
grid;


subplot(3,1,3);
plot(1:MaxIt,min(it_lS_exp)*4*dt,1:MaxIt,min(it_S_exp)*4*dt);
title('Mínimo de pasos por iteración en 10 iteraciones')
legend('Lambda-SARSA','SARSA')
xlabel('Iteración')
ylabel('Tiempo (s)')
grid;


cd(Carpeta);
print('comp_SARSA','-djpeg')
cd ..


% close all