function [N_path, lam_path] = surv_path(N_0, lam_0, T, Num_path, c, eta, w)

% generate random survival prob e mortality intensity
% Vigna e Luciano 2008, Barigou-Chen-Dhaene-2019

var_lam = eta^2*(exp(2*c)-1)/(2*c);
lam_path = zeros(T,Num_path);
lam_path(1,:) = lam_0*exp(c)+sqrt(var_lam)*w(1,:);
% Lambda(t) = int_t^t+1 lambda(s) ds
var_int = 2*(eta^2/(2*c^2) - eta^2*exp(c)/c^3+eta^2*exp(2*c)/(4*c^3)+3*eta^2/(4*c^3));
int_path = zeros(T,Num_path);
int_path(1,:) = lam_0*(exp(c)-1)/c + sqrt(var_int)*w(1,:);
for t=2:T
    lam_path(t,:) = lam_path(t-1,:)*exp(c)+sqrt(var_lam)*w(t,:);
    int_path(t,:) = lam_path(t-1,:)*(exp(c)-1)/c + sqrt(var_int)*w(t,:);
end
%death prob p(t) = 1- exp(-Lambda(t)) 
p_death = 1- exp(-int_path);
%survival individuals N(t) = N(t-1) - D(t) with  D(t) = binomial(N(t-1), p(t))
N_path = zeros(T,Num_path);
N_path(1,:) = N_0 - binornd(N_0*ones(1,Num_path),p_death(1,:));
for t=2:T
   N_path(t,:) = N_path(t-1,:) - binornd(N_path(t-1,:),p_death(t,:)); 
end