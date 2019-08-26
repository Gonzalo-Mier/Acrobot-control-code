clear

Carpeta = 'PD';
medidas_est = 0;1;
error_m = 0; 0.05;
t_total = 60;20;
dt_controller = 0.05;


LIM_MAX = 10000;
MAX_EXP = 1;2;
F_max = 1;4.5;600;
n_ind = 60;
n_it = 40;

% ind_rand = @() 7*rand(1,4)-2;
% f_eval_ind =@(x)  10.^x-10^3;
% 
% f_max_ind = @(x) x.*rand(1,3)*2-x;
% ind_rand = @() f_max_ind([1,60,5]);
% f_eval_ind =@(x)  x;
% 
% ind_rand = @() (2*rand(1,3)-1).*[2 10 4] +[1 44 2.1];
% f_eval_ind =@(x)  x;

f_max_ind = @(x) x.*rand(1,3);
ind_rand = @() f_max_ind([1,60,5]);
f_eval_ind =@(x)  x;


n = round(n_ind/3);
area = n*(n+1)/2;
prob = cumsum(n:-1:1)/area;



tic

% % % Inicio de parametros % % %
q1_ini = 0.01;
q2_ini = 0;
dq1_ini = 0;
dq2_ini = 0;
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
dt = dt_controller/4;


if medidas_est == 1
% % % Parametros iniciales estimados del problema
m1_est = 0.98;
m2_est = 1.03;
lc1_est = 0.51;
lc2_est = 0.48;
l1_est = 1.01;
l2_est = 0.96;
I1_est = 1.05;
I2_est = 1.01;
else
% % % Parametros iniciales estimados del problema
m1_est = 1;
m2_est = 1;
lc1_est = 0.5;
lc2_est = 0.5;
l1_est = 1;
l2_est = 1;
I1_est = 1;
I2_est = 1;
end
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
    AG_PD
    fitness_ex(experimento,:) = fitness_it;
    INDIV{experimento} = individuos{1};
    toc
end

% figure(1);
% for experimento = 1:MAX_EXP
%     plot( fitness_ex(experimento,:) );
%     title('Fitness por iteración')
%     xlabel('Iteración')
%     ylabel('Fitness')
%     hold on;
% end
% hold off;
% % ylim([0,LIM_MAX]);
% 
% 
% figure(2);
% plot(sum(fitness_ex)/MAX_EXP);
% title('Media de Fitness por iteración en 10 iteraciones')
% xlabel('Iteración')
% ylabel('Fitness')
% % ylim([0,LIM_MAX]);
% 
% 
% figure(3);
% plot(1:n_it,sum(fitness_ex)/MAX_EXP,1:n_it,max(fitness_ex),1:n_it,min(fitness_ex));
% title('Media, máximo y mínimo de fitness por iteración en 10 iteraciones')
% legend('Media','maximo','minimo')
% xlabel('Iteración')
% ylabel('Fitness')
% % ylim([0,LIM_MAX]);

% mkdir(Carpeta);
% cd(Carpeta);
% for i = 1:3
%     saveas(i,strcat(Carpeta,'_',char(string(i))),'jpeg')
% end
% cd ..

[~,A] = min(fitness_ex(:,end));
K = f_eval_ind(INDIV{A});

PD
grafica_PD_ES

toc