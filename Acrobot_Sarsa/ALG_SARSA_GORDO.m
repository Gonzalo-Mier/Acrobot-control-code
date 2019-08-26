
clear
tic
MAX_EXP = 2;

var_1_cont = 1;

max_div = 6;



for i_alpha = 0:1/max_div:1
    var_2_cont = 1;
    for i_gamma = 0:1/max_div:1
        
        Carpeta = strcat('Im_S_Sutton_i_',char(string(var_1_cont)),'_',char(string(var_2_cont)));
        
        alpha = i_alpha/48;
        lambda = 0.9;
        gamma = i_gamma;
        eps = 0;
        MaxIt = 11;600;
        t_total = 10; 1000;
        
        EXP_COMP_SARSA;
        
        [min_t_lS_v(var_1_cont,var_2_cont), min_t_lS_p(var_1_cont,var_2_cont)] = min(min(it_lS_exp)*4*dt);
        [min_t_S_v(var_1_cont,var_2_cont), min_t_S_p(var_1_cont,var_2_cont)] = min(min(it_S_exp)*4*dt);
        
        
        
        var_2_cont = var_2_cont + 1;
        'var1:',var_1_cont,'var2:', var_2_cont
        
        toc
        
    end
    var_1_cont = var_1_cont + 1;
    
end

