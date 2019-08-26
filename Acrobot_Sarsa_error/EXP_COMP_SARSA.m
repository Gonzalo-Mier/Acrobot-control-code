clear


MAX_EXP = 1;
error_m = 0;0.05;


alpha = 0.2/48;
lambda = 0.9;
gamma = 1;
eps = 0;
MaxIt = 600;
t_total = 4000;

experimentos_SARSA
% experimentos_lambda_SARSA
% 
% 
% grafica_SARSA_ES

grafica_SARSA_EN


% figure(1);
% subplot(3,1,1);
% plot(1:MaxIt,sum(it_lS_exp)/MAX_EXP,1:MaxIt,sum(it_S_exp)/MAX_EXP);
% title('Media de pasos por iteración en 10 iteraciones')
% legend('Lambda-SARSA','SARSA')
% xlabel('Iteración')
% ylabel('Pasos necesitados para completar la iteración')
% 
% subplot(3,1,2);
% plot(1:MaxIt,max(it_lS_exp),1:MaxIt,max(it_S_exp));
% title('Máximo de pasos por iteración en 10 iteraciones')
% legend('Lambda-SARSA','SARSA')
% xlabel('Iteración')
% ylabel('Pasos necesitados para completar la iteración')
% 
% subplot(3,1,3);
% plot(1:MaxIt,min(it_lS_exp),1:MaxIt,min(it_S_exp));
% title('Mínimo de pasos por iteración en 10 iteraciones')
% legend('Lambda-SARSA','SARSA')
% xlabel('Iteración')
% ylabel('Pasos necesitados para completar la iteración')


% saveas(1,'Figuras/comp_SARSA','jpeg')


% close all