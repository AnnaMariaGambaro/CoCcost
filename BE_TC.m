function [ESCR,L] = BE_TC(T, r, S_path, N_path,  dg, K,delta)

phi = norminv(0.995);
Num_path = size(N_path,2);
BE = zeros(T,Num_path); 
BE(T,:) = max(S_path(T),K).*N_path(T,:);
ESCR = zeros(T,1); 
for t = T-1:-1:1   
    base = createPolBase(S_path(t+1).*N_path(t,:)',dg);   
    coeffBE = mvregress(base,BE(t+1,:)'*exp(-r));    
    BE(t,:) = (base*coeffBE)';    
     coeffBE2 = mvregress(base,(BE(t+1,:)'-BE(t,:)').^2);  
     var_BE = (base*coeffBE2)';
    SCR_BE(t,:) = phi*sqrt(max(var_BE,0));
    ESCR(t) = exp(-r*t)*mean(SCR_BE(t,:));
    BE(t,:) = BE(t,:)+delta*SCR_BE(t,:);
end
L=exp(-r)*mean(BE(1,:));
