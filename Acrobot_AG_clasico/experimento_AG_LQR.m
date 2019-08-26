clear

Carpeta = 'LQR';

LIM_MAX = 1;
MAX_EXP = 10;
F_max = 10;4.5;600;
n_ind = 60;
n_it = 40;

% ind_rand = @() 7*rand(1,4)-2;
% f_eval_ind =@(x)  10.^x-10^3;
% 
% f_max_ind = @(x) x*rand(1,4)-x/2;
% ind_rand = @() f_max_ind(120*2);
% f_eval_ind =@(x)  x;
% 
f_max_ind = @(x) x*2*(rand(1,4)-0.5)+[-182.3015  -45.5932  -85.9921  -28.0590];
ind_rand = @() f_max_ind(10);
f_eval_ind =@(x)  x;




tic

% % % Inicio de parametros % % %
q1_ini = pi+.04;
q2_ini = -0.05;
dq1_ini = -0.02;
dq2_ini = 0.04;
ddq1_ini = 0;
ddq2_ini = 0;
F_ini = 0;


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
dt_controller = 0.01;
dt = dt_controller/4;
t_total = 20;


% % % % Parametros iniciales del problema
% m1 = 1.9008;
% m2 = 0.7175;
% lc1 = 1.8522e-1;
% lc2 = 6.2052e-2;
% l1 = .2;
% l2 = .2;
% I1 = 4.3399e-3;
% I2 = 5.2285e-3;
% g = 9.8;
% dt_controller = 0.01;
% dt = dt_controller/4;
% t_total = 10;



for experimento = 1:MAX_EXP
    experimento
    AG_LQR
    fitness_ex(experimento,:) = fitness_it;
    INDIV{experimento} = individuos{1};
    toc
end

figure(1);
for experimento = 1:MAX_EXP
    plot( fitness_ex(experimento,:) );
    title('Fitness por iteración')
    xlabel('Iteración')
    ylabel('Fitness')
    hold on;
end
hold off;
ylim([0,LIM_MAX]);


figure(2);
plot(sum(fitness_ex)/MAX_EXP);
title('Media de Fitness por iteración en 10 iteraciones')
xlabel('Iteración')
ylabel('Fitness')
ylim([0,LIM_MAX]);


figure(3);
plot(1:n_it,sum(fitness_ex)/MAX_EXP,1:n_it,max(fitness_ex),1:n_it,min(fitness_ex));
title('Media, máximo y mínimo de fitness por iteración en 10 iteraciones')
legend('Media','maximo','minimo')
xlabel('Iteración')
ylabel('Fitness')
ylim([0,LIM_MAX]);

mkdir(Carpeta);
cd(Carpeta);
for i = 1:3
    saveas(i,strcat('LQR_',char(string(i))),'jpeg')
end
cd ..

[~,A] = min(fitness_ex(:,end));
K = f_eval_ind(INDIV{A});
% K = [-113.8908, -9.507, -17.3951, -1.9405];

LQR_Control

grafica_res

toc