
for experimento = 1:MAX_EXP
    experimento
    Sarsa_Acrobot
    it_S_exp(experimento,:) = it_t;
end

close all

for experimento = 1:MAX_EXP
    for i=1:MaxIt-9
        it_S_exp_median(experimento,i) = sum(it_S_exp(experimento,i:i+9))/10;
    end
end

figure(1);
for experimento = 1:MAX_EXP
    plot(it_S_exp(experimento,:)*4*dt);
    title('Pasos por iteración')
    xlabel('Iteración')
ylabel('Tiempo (s)')
    hold on;
end
grid;

hold off;

figure(2);

for experimento = 1:MAX_EXP
    plot(it_S_exp_median(experimento,:)*4*dt);
    title('Pasos por iteración suavizado en 10 pasos')
    xlabel('Iteración')
ylabel('Tiempo (s)')
    hold on;
end
grid;

hold off;


figure(3);
sum_it_S_exp = zeros(size(it_S_exp(experimento,:)));
S_conv = 0;
for experimento = 1:MAX_EXP 
    if it_S_exp(experimento,end) < 1e4
        sum_it_S_exp = sum_it_S_exp + it_S_exp(experimento,:);
        S_conv = S_conv + 1;
        it_S_conv_exp(S_conv,:) = it_S_exp(experimento,:);
    end
end
plot(sum_it_S_exp/S_conv*4*dt);
title('Media de pasos por iteración en 10 iteraciones')
xlabel('Iteración')
ylabel('Tiempo (s)')
grid;


figure(4);
plot(1:MaxIt,sum(it_S_conv_exp)/S_conv*4*dt,1:MaxIt,max(it_S_conv_exp)*4*dt,1:MaxIt,min(it_S_conv_exp)*4*dt);
title('Media, máximo y mínimo de pasos por iteración en 10 iteraciones')
legend('media','maximo','minimo')
xlabel('Iteración')
ylabel('Tiempo (s)')
grid;



for i = 1:4
    saveas(i,strcat(Carpeta,'/SARSA_',char(string(i))),'jpeg')
end

close all

% close all


