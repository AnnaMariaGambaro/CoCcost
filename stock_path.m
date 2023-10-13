function [S_path, w] = stock_path(S_0, r, sigma,w)

%% risk neutral path of stock price - lognomral dynamic
dt=1;
mu = (r-0.5*sigma^2)*dt;

S_path = S_0*exp(cumsum(mu+sigma*w)); 