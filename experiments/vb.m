%% setup
addpath('../WSBM_v1.2');

%% set experimental parameters
K = 3;
n = 20;
t0 = 1;

%% generate data
z = datasample(1:K,n);
z = z';                                 % make it a column vector

th = gamrnd(1,1,[K,K])

A = zeros(n,n);

for i=1:n
    for j=(i+1) : n
        A(i,j) = exprnd(th(z(i),z(j)),1);
    end
end

A = A + A';                              % make symmetric

%% test recovery


[labels,model] =  wsbm(A,3,'W_Distr','Exponential','E_Distr','none');