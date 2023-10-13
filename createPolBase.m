function Base = createPolBase(X,dg)
%create a polinomial base
% X = matrix MxN, X = (X_1,X_2,...X_N)
%dg=polinomial degree
% Base = polinomial base
%for instance if dg = 2 and N = 3
% Base = [ones(size(X1)), X1, X2, X3, X1.^2, X2.^2, X3.^2, X1.*X2, X1.*X3, X2.*X3]

Base = ones(size(X(:,1)));
for dg_idx =1:dg 
    Nmatrix = multiNomial(dg_idx,size(X,2),'part');
    for mm = 1:size(Nmatrix,1)
        Base_aux = prod(X.^repmat(Nmatrix(mm,:),size(X,1),1),2);
        Base = [Base,Base_aux];
    end
end
