function SCR_BE = SCR_EIOPA(T, r, S_path, N_path,  dg,delta, K)


phi = norminv(0.995);
Num_path = size(N_path,2);
BE = zeros(T,Num_path); 

SCR_BE = zeros(T,1); 

BE(T,:) = max(S_path(T),K).*N_path(T,:);
for t = T-1:-1:1   
    base = createPolBase(S_path(t+1).*N_path(t,:)',dg);   
    coeffBE = mvregress(base,BE(t+1,:)'*exp(-r));    
    BE(t,:) = (base*coeffBE)';
    
    % Var_N = Var[N(t+1)|N(t) = mean(N(t))]
    baseN = createPolBase(N_path(t,:)',2);
    coeff2 = mvregress(baseN, N_path(t+1,:)'.^2); 
    coeff1 = mvregress(baseN, N_path(t+1,:)');
    baseN_mean = createPolBase(mean(N_path(t,:)),2); 
    Var_N = (baseN_mean*coeff2)'-(baseN_mean*coeff1)'.^2;
    N_tilde = mean(N_path(t,:)) + phi*sqrt(max(Var_N,0)); 
    
    base2 = createPolBase(S_path(t+1).*N_tilde,dg); 
    base3 = createPolBase(S_path(t+1).*mean(N_path(t,:)),dg); 
    SCR_BE(t) = exp(-r*t)*((base2*coeffBE)-(base3*coeffBE));
    BE(t,:) = BE(t,:) +delta*SCR_BE(t,:);
end
