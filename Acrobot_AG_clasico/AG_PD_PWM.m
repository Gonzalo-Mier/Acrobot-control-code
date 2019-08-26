% clear; clc; close all;
% 
% tic
% n_ind = 50;
% n_it = 60;
% F_max = 100;inf;4.5;

n = round(n_ind/3);
area = n*(n+1)/2;
prob = cumsum(n:-1:1)/area;

% % % % Parametros iniciales del problema
% m1 = 1;
% m2 = 1;
% lc1 = 0.5;
% lc2 = 0.5;
% l1 = 1;
% l2 = 1;
% I1 = 1;
% I2 = 1;
% g = 9.8;
% dt_controller = 0.002;
% dt = dt_controller/4;
% t_total = 30;


% % % Inicio de parametros % % %
q1_ini = pi-rand*.1;
q2_ini = 0;
dq1_ini = 0;
dq2_ini = 0;
ddq1_ini = 0;
ddq2_ini = 0;
F_ini = 0;

F = F_ini; q1 = q1_ini; q2 = q2_ini; dq1 = dq1_ini; dq2 = dq2_ini;
F_t = [];
q1_t = [];
q2_t = [];


p_ini = [q1_ini, q2_ini, dq1_ini, dq2_ini, ddq1_ini, ddq2_ini, F_ini];
params = [m1, m2, lc1, lc2, l1, l2, I1, I2, g, dt_controller, dt, t_total, F_max];


for i=1:n_ind
    individuos{i} = ind_rand();
end


for it = 1: n_it
    it
    for i=1:n_ind
        fitness(i) = f_PD_up(f_eval_ind(individuos{i}), p_ini, params);
    end
    [~,pos] = sort(fitness);
    
    for i=1:n_ind
        ind_sort{i} = individuos{pos(i)};
    end
    
    for i=1:round(n_ind/3)
        
        padre1 = ind_sort{find(rand <= prob, 1,'first')};
        padre2 = ind_sort{find(rand <= prob, 1,'first')};
        
        
        a = rand(1,3)*2-0.5;
        hijo1 = min(padre1,padre2) + a.*(max(padre1,padre2) - min(padre1,padre2));
        ind_sort{i+round(n_ind/3)} = hijo1;
    end
    
    for i=1+2*round(n_ind/3):n_ind
        ind_sort{i} = ind_rand();
    end
    individuos = ind_sort;
    fitness(1:10)
end


toc