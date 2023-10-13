function [value,price,Nsurv] = BE_indip(T,K,S_0,r,sig,N_0,lam_0,c,eta)
% Best Estimate of liabilities in case of indipendence 
% between stock S and force of mortality lambda
% exp(-r*T)*E^Q[ E^P[ max(S_T,K)*N(T) | S_T, N(t),lam_t] S_t] 
% = exp(-r*T)*E^Q[ max(S_T,K)*N(T) | S_t] * E^P[ N(T) | N(t),lam_t]

price = zeros(length(T),1);
for nn = 1:length(T)
    %financial price
    price(nn) = exp(-r*T(nn))*K+blsprice(S_0, K, r, T(nn), sig);   
end
%expected number of survivors 
Nsurv = N_0*SurvFunc(T',lam_0,c, eta);
value = price.*Nsurv;