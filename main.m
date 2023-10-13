%%
% "The Capital-on-Capital Cost in Solvency II Risk Margin."
% Anna Maria Gambaro 
% October 2023
% Available at SSRN: https://ssrn.com/abstract=4418565
%%
close all;
clear;
clc;
%% Parameters from Barigou-Chen-Dhaene-2019-IME, section 5.2 and 5.3
% financial parameters
S_0 = 1;
K = 1; %0;
sig = 0.1;
r = 0.01;
mu = 0.02;
%mortality parameters: Luciano et al 2017 
%UK male inividuals 55 years old at t=0
lam_0 = 0.0087;
c = 0.0750;
eta = 0.000597;
N_0 = 1000;
%other parameters
rho= [0,0.5,0.75,1];
T = 1:40; %1:30; 
delta = 6/100; %cost of capital rate
phi = norminv(0.995);%2.58; % phi^(-1)(alpha) with alpha =99.5%
%LSMC parameters
Num_pathS = 100;
Num_pathN = 1000;
dg=2; %polynomial degree

if ~exist(['fig-K-',num2str(K)], 'dir')
   mkdir(['fig-K-',num2str(K)])
end

%% analytic BE for rho=0
[BE,price,Nsurv] = BE_indip(T,K,S_0,r,sig,N_0,lam_0,c,eta);
disp([T',BE,price,Nsurv])

figure
fig1 = plot(T',BE);
xlabel('maturity (years)')
ylabel('Best Estimate')
savefig(['fig-K-',num2str(K),'/BE_indip.fig'])
saveas(fig1,['fig-K-',num2str(K),'/BE_indip.png'])

%% Least square Monte Carlo
w = randn(T(end),Num_pathS);
S_path = stock_path(S_0, r, sig,w);

BEMC = zeros(T(end),length(rho));
ASCR = zeros(T(end),length(rho));
L_TC = zeros(T(end),length(rho));
L_TC2 = zeros(T(end),length(rho));
RM = zeros(T(end),length(rho));
RM_hat = zeros(T(end),length(rho));
RM_NTC = zeros(T(end),length(rho));

for jj = 1:length(rho)
disp(rho(jj))
BEMC_aux = zeros(T(end),Num_pathS);
ASCR_aux = zeros(T(end),Num_pathS);
SCR_NTC_aux = zeros(T(end),Num_pathS);
SCR_TC_aux = zeros(T(end),Num_pathS);
ESCR = zeros(T(end),Num_pathS);
ESCR_L = zeros(T(end),Num_pathS);
SCR_NTC = zeros(T(end),Num_pathS);
L_TC_aux = zeros(T(end),Num_pathS);
    for kk=1:Num_pathS
    z = randn(T(end),Num_pathN);
    w2 = repmat(w(:,kk),1,Num_pathN)*rho(jj)+z*sqrt(1-rho(jj)^2);
    N_path = surv_path(N_0, lam_0, T(end), Num_pathN, c, eta,w2);
    BEMC_aux(:,kk)= mean(max(repmat(S_path(:,kk),1,Num_pathN),K).*N_path,2);

    for t = 1:T(end)
    ESCR(1:t,kk) = ESCR_fun(t, r, S_path(:,kk), N_path, dg, K); 
    ASCR_aux(t,kk) = sum(ESCR(1:t,kk));
    SCR_NTC(1:t,kk) = SCR_EIOPA(t, r, S_path(:,kk), N_path, dg, K, delta);
    SCR_NTC_aux(t,kk) = sum(SCR_NTC(1:t,kk));
    [ESCR_L(1:t,kk), L_TC_aux(t,kk)] = BE_TC(t, r, S_path(:,kk), N_path,  dg, K,delta);
    SCR_TC_aux(t,kk) = sum(ESCR_L(1:t,kk));
    end

    end
ASCR(:,jj) = mean(ASCR_aux,2);   
RM_hat(:,jj) = delta*ASCR(:,jj);
RM_NTC(:,jj) = delta*mean(SCR_NTC_aux,2); 
BEMC(:,jj) = exp(-r*T').*mean(BEMC_aux,2);
RM(:,jj) = delta*mean(SCR_TC_aux,2);
L_TC(:,jj) = BEMC(:,jj)+delta*mean(SCR_TC_aux,2);
L_TC2(:,jj) = mean(L_TC_aux,2);

disp("liabilitie values time consistent")
disp([L_TC2(:,jj),L_TC(:,jj)])%verify Proposition 1
disp("RM, RM_hat, RM_NTC ")
disp([RM(:,jj), RM_hat(:,jj), RM_NTC(:,jj)])
disp("BE, L_ASCR, L_TC, L_NTC")
disp([BEMC(:,jj),BEMC(:,jj)+delta*ASCR(:,jj),L_TC(:,jj) ,BEMC(:,jj)+RM_NTC(:,jj)])

end


%% Figures
style_vec = ['-','--','-.',':'];
figure
for jj = 1:length(rho)
    plot(T',BEMC(:,jj),style_vec(jj),'DisplayName',['$\rho$ = ', num2str(rho(jj))])
    hold on
end
xlabel('maturity (years)')
ylabel('Best Estimate')
legend('interpreter','latex','Location','Best')
savefig(['fig-K-',num2str(K),'/BE.fig'])
saveas(gcf,['fig-K-',num2str(K),'/BE.png'])

%%
for jj = 1
figure
plot(T',BEMC(:,jj),'k-','DisplayName','$L^{BE}$')
hold on
plot(T',L_TC(:,jj),'-*','DisplayName','$L$')
plot(T',BEMC(:,jj)+RM_hat(:,jj),'-o','DisplayName','$\widehat{L}$')
plot(T',BEMC(:,jj)+RM_NTC(:,jj),'-s','DisplayName','$L$ NTC')
hold off
xlabel('maturity (years)')
ylabel('fair value')
legend('interpreter','latex','Location','Best')
title(['\rho = ', num2str(rho(jj))] )
%ylim([300 1100])
savefig(['fig-K-',num2str(K),'/fair-rho',num2str(rho(jj)*100),'.fig'])
saveas(gcf,['fig-K-',num2str(K),'/fair-rho',num2str(rho(jj)*100),'.png'])
end

%%
for jj = 1:length(rho)
figure
plot(T',RM(:,jj),'-*','DisplayName','$RM$')
hold on
plot(T',RM_hat(:,jj),'-o','DisplayName','$\widehat{RM}$')
plot(T',RM_NTC(:,jj),'-s','DisplayName','$RM$ NTC')
hold off
xlabel('maturity (years)')
ylabel('risk margin')
legend('interpreter','latex','Location','Best')
title(['\rho = ', num2str(rho(jj))] )
ylim([0 60])
savefig(['fig-K-',num2str(K),'/RM-rho',num2str(rho(jj)*100),'.fig'])
saveas(gcf,['fig-K-',num2str(K),'/RM-rho',num2str(rho(jj)*100),'.png'])
end

%%
for jj = 1:length(rho)
figure
plot(T',100*RM(:,jj)./BEMC(:,jj),'-*','DisplayName','$RM$')
hold on
plot(T',100*RM_hat(:,jj)./BEMC(:,jj),'-o','DisplayName','$\widehat{RM}$')
plot(T',100*RM_NTC(:,jj)./BEMC(:,jj),'-s','DisplayName','$RM$ NTC')
hold off
xlabel('maturity (years)')
ylabel('risk loading (%)')
legend('interpreter','latex','Location','Best')
title(['\rho = ', num2str(rho(jj))] )
ylim([0 40])
savefig(['fig-K-',num2str(K),'/RL-rho',num2str(rho(jj)*100),'.fig'])
saveas(gcf,['fig-K-',num2str(K),'/RL-rho',num2str(rho(jj)*100),'.png'])
end

figure
hold on
for jj = 1:length(rho)
plot(T(2:end)',100*max(RM(2:end,jj)-RM_hat(2:end,jj),0)./RM(2:end,jj),'DisplayName',['$\rho$ = ', num2str(rho(jj))])
end
hold off
xlabel('maturity (years)')
ylabel('C / RM (%)')
legend('interpreter','latex','Location','Best')
title('Cost of Capital-on-Capital relative to RM' )
savefig(['fig-K-',num2str(K),'/Cost-CoC.fig'])
saveas(gcf,['fig-K-',num2str(K),'/Cost-CoC.png'])

