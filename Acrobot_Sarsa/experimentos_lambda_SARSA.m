
for experimento = 1:MAX_EXP
    experimento
    lambda_Sarsa_Acrobot
    it_lS_exp(experimento,:) = it_t;
end

close all

for experimento = 1:MAX_EXP
    for i=1:MaxIt-9
        it_lS_exp_median(experimento,i) = sum(it_lS_exp(experimento,i:i+9))/10;
    end
end

figure(1);
for experimento = 1:MAX_EXP
    plot(it_lS_exp(experimento,:)*4*dt);
    title('Pasos por iteración')
    xlabel('Iteración')
ylabel('Tiempo (s)')
    hold on;
end
grid;

hold off;

figure(2);

for experimento = 1:MAX_EXP
    plot(it_lS_exp_median(experimento,:)*4*dt);
    title('Pasos por iteración suavizado en 10 pasos')
    xlabel('Iteración')
ylabel('Tiempo (s)')
    hold on;
end
grid;

hold off;


figure(3);
sum_it_lS_exp = zeros(size(it_lS_exp(experimento,:)));
lS_conv = 0;
for experimento = 1:MAX_EXP 
    if it_lS_exp(experimento,end) < 1e4
        sum_it_lS_exp = sum_it_lS_exp + it_lS_exp(experimento,:);
        lS_conv = lS_conv + 1;
        it_lS_conv_exp(lS_conv,:) = it_lS_exp(experimento,:);
    end
end
plot(sum_it_lS_exp/lS_conv*4*dt);
title('Media de pasos por iteración en 10 iteraciones')
xlabel('Iteración')
ylabel('Tiempo (s)')
grid;


figure(4);
plot(1:MaxIt,sum(it_lS_conv_exp)/lS_conv*4*dt,1:MaxIt,max(it_lS_conv_exp)*4*dt,1:MaxIt,min(it_lS_conv_exp)*4*dt);
title('Media, máximo y mínimo de pasos por iteración en 10 iteraciones')
legend('media','maximo','minimo')
xlabel('Iteración')
ylabel('Tiempo (s)')
grid;




for i = 1:4
    saveas(i,strcat(Carpeta,'/lambda_SARSA_',char(string(i))),'jpeg')
end

close all



