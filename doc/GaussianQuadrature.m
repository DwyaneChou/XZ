clc
clear

ngpts = 5;

% Calculate Gauss-Laguerre Quadrature Evaluation Points and Weights
syms x
syms t
gpts_c = vpasolve(legendreP(ngpts,x) == 0); % Gaussian points
p      = 1:ngpts;
for k=1:ngpts
    j        = setdiff(p,k);
    wts_c(k) = int(prod((t - gpts_c(j)) ./ (gpts_c(k) - gpts_c(j))), t, -1, 1);
    
    disp('Gaussian quadrature points')
    disp(gpts_c(k))
    
    disp('Gaussian quadrature weights')
    disp(wts_c(k))
end
