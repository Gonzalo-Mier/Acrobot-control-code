% MAX_EXP = 10;

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
    plot(it_lS_exp(experimento,:));
    title('Pasos por iteración')
    xlabel('Iteración')
    ylabel('Pasos necesitados para completar la iteración')
    hold on;
end
hold off;

figure(2);

for experimento = 1:MAX_EXP
    plot(it_lS_exp_median(experimento,:));
    title('Pasos por iteración suavizado en 10 pasos')
    xlabel('Iteración')
    ylabel('Pasos necesitados para completar la iteración')
    hold on;
end
hold off;


figure(3);
plot(sum(it_lS_exp)/MAX_EXP);
title('Media de pasos por iteración en 10 iteraciones')
xlabel('Iteración')
ylabel('Pasos necesitados para completar la iteración')

figure(4);
plot(sum(it_lS_exp_median)/MAX_EXP);
title('Media de pasos por iteración en 10 iteraciones con suavizado en 10 pasos')
xlabel('Iteración')
ylabel('Pasos necesitados para completar la iteración')


figure(5);
plot(1:MaxIt,sum(it_lS_exp)/MAX_EXP,1:MaxIt,max(it_lS_exp),1:MaxIt,min(it_lS_exp));
title('Media, máximo y mínimo de pasos por iteración en 10 iteraciones')
legend('media','maximo','minimo')
xlabel('Iteración')
ylabel('Pasos necesitados para completar la iteración')



for i = 1:5
    saveas(i,strcat('Figuras/lambda_SARSA_',char(string(i))),'jpeg')
end

close all



