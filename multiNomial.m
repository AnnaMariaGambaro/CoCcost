function [Nmatrix Ncoef]=multiNomial(n,l,opt)
% MULTINOMIAL determines the matrix of power for a multinomial funtion as well
% as the corresponding co-efficient
% (x1 + x2 + ...+ xl)^n
% n- power of the multinomial expansion
% l- number of terms in the multinomial expresion
% opt- {'part','all'} for selectivity of the outcome
% Nmatrix - matrix of powers
% Ncoef - array of coefficients
if nargin<3
    opt='part';
end
nn=zeros((n+1)^l,l);
h=0:n;
for ii=1:l
    p=l-ii+1;
    hi=kron(ones((n+1)^(ii-1),(n+1)^(l-ii)),h);
    xx=reshape(hi,numel(hi),1);
    nn(:,p)=xx;
end
switch opt
    case 'part'
        nSum=sum(nn,2);
        Nmatrix=nn(nSum==n,:);
    case 'all'
        Nmatrix=nn;
    otherwise
        error('opt can only be ''part'' or ''all'' ')
end
Ncoef=factorial(n)*(prod(factorial(Nmatrix),2)).^(-1);