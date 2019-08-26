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
    title('Pasos por iteraci�n')
    xlabel('Iteraci�n')
    ylabel('Pasos necesitados para completar la iteraci�n')
    hold on;
end
hold off;

figure(2);

for experimento = 1:MAX_EXP
    plot(it_lS_exp_median(experimento,:));
    title('Pasos por iteraci�n suavizado en 10 pasos')
    xlabel('Iteraci�n')
    ylabel('Pasos necesitados para completar la iteraci�n')
    hold on;
end
hold off;


figure(3);
plot(sum(it_lS_exp)/MAX_EXP);
title('Media de pasos por iteraci�n en 10 iteraciones')
xlabel('Iteraci�n')
ylabel('Pasos necesitados para completar la iteraci�n')

figure(4);
plot(sum(it_lS_exp_median)/MAX_EXP);
title('Media de pasos por iteraci�n en 10 iteraciones con suavizado en 10 pasos')
xlabel('Iteraci�n')
ylabel('Pasos necesitados para completar la iteraci�n')


figure(5);
plot(1:MaxIt,sum(it_lS_exp)/MAX_EXP,1:MaxIt,max(it_lS_exp),1:MaxIt,min(it_lS_exp));
title('Media, m�ximo y m�nimo de pasos por iteraci�n en 10 iteraciones')
legend('media','maximo','minimo')
xlabel('Iteraci�n')
ylabel('Pasos necesitados para completar la iteraci�n')



for i = 1:5
    saveas(i,strcat('Figuras/lambda_SARSA_',char(string(i))),'jpeg')
end

close all



