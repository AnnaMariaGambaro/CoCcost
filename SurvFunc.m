function res = SurvFunc(t,lam_0,a, sig)

%Vigna e Luciano 2008 eq. 3.1 and 3.6
%d lambda_t = a*lambda_t dt + sig dW_t

beta = (1-exp(a*t))/a;
alpha = sig^2*t/(2*a^2) - sig^2*exp(a*t)/a^3+sig^2*exp(2*a*t)/(4*a^3)+3*sig^2/(4*a^3);
res = exp(alpha+beta*lam_0);