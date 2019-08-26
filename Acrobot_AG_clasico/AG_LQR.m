% clear;
% clc; close all;


ind_rand = @() 7*rand(1,4)-2;
f_eval_ind =@(x)  10.^x-10^3;




n = round(n_ind/3);
area = n*(n+1)/2;
prob = cumsum(n:-1:1)/area;



params = [m1, m2, lc1, lc2, l1, l2, I1, I2, g, dt_controller, dt, t_total, F_max];


for i=1:n_ind
    %     individuos{i} = (2*rand(1,4)-1)*1e3;
    individuos{i} = ind_rand();
end



for it = 1: n_it
    %     it
    
    
    p_ini = [q1_ini, q2_ini, dq1_ini, dq2_ini, ddq1_ini, ddq2_ini, F_ini];
    for i=1:n_ind
        fitness(i) = f_LQR(f_eval_ind(individuos{i}), p_ini, params);
    end
    [~,pos] = sort(fitness);
    
    for i=1:n_ind
        ind_sort{i} = individuos{pos(i)};
    end
    
    for i=1:round(n_ind/3)
        
        padre1 = ind_sort{find(rand <= prob, 1,'first')};
        padre2 = ind_sort{find(rand <= prob, 1,'first')};
        
        
        a = rand(1,4)*2-0.5;
        %         a = rand(1,4)*10-4.5;
        hijo1 = min(padre1,padre2) + a.*max((max(padre1,padre2) - min(padre1,padre2)),1);
        ind_sort{i+round(n_ind/3)} = hijo1;
    end
    
    for i=1+2*round(n_ind/3):n_ind
        %         ind_sort{i} = (2*rand(1,4)-1)*1e3;
        ind_sort{i} = ind_rand();
        
    end
    individuos = ind_sort;
    
    fitness(1:10);
    fitness_it(it) = fitness(1);
end

