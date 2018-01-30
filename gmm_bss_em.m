function [model] = gmm_bss_em(X,K,varargin)
%%%%%%%% EM-algorithm for complex rank-1 Gaussian Mixture Estimation %%%%%%
% The model is:
%     p(x|z=k) = N(x;0,a_k*a_k'+sigma2_k*eye(D))
%     p(z=k) = pi_k
% where a_k \in C^D, sigma2_k>0, pi_k\in(0,1), sum(pi) = 1
%%% INPUT:
% X (D*N) : N input D-dimensional complex vectors
% K (int) : number of gaussians
% maxiter (int) : max number of iterations (default: 100)
% tol (double) : convergence tolerance (default: 1e-3)
%%% OUTPUT:
% label (1*N) \in {1..K} : label of each data point
% model {struct}:
%    .A (D*K) : each column a_k is source k's steering vector
%    .sigma2 (K*1) : noise variance for each source (>0)
%    .w (K*1) : cluster weights ( sum(w) = 1 )
% R (K*N) : posterior probabilities p(z_n=z|x_n;model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Parameters

p=inputParser();
p.KeepUnmatched=true;
p.addParameter('tol',1e-3,@isnumeric);
p.addParameter('verb',0,@isnumeric);
p.addParameter('maxiter',100,@isnumeric);

p.parse(varargin{:});
options=p.Results;

N = size(X,2);

%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize steering vectors with K random data points:
p = randperm(N);
model.A = X(:,p(1:K));
% Initialize noise variance to half the total data variance:
model.sigma2 = std(X(:))/2*ones(K,1);
% Initialize cluster weights as equal
model.w = 1/K*ones(K,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
old_R = ones(K,N);
converged = 0;
niter = 0;
if(options.verb>0);fprintf(1,'Rank-1 CGM EM Running... iter');end;
while(~converged)
    niter = niter + 1;
    if(options.verb>0);fprintf(1,' %d',niter);end;
    % Expectation
    R = expectation(X,model); % N*K : posterior probabilities p(z|x;model)
    
    % Maximization
    model = maximization(X,R); % updated model parameters
    
    % Check convergence
    if norm(old_R(:)-R(:))/norm(old_R) < options.tol || niter>=options.maxiter
        converged = true;
    end
    old_R = R;
end
if(options.verb>0);fprintf(1,'\n');end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WEIRD NORMALIZATION ??
model.A = model.A/sqrt(2);
model.sigma2 = model.sigma2/2;
model.alpha = 2*ones(1,K);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXPECTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT:
% X (D*N) : N input D-dimensional complex vectors
% model (struct) : The current parameter estimated
%    .A (D*K) : each column a_k is source k's steering vector
%    .sigma2 (K*1) : noise variance for each source (>0)
%    .w (K*1) : cluster weights ( sum(w) = 1 )
%%% OUTPUT:
% R (K*N) : posterior probabilities: p(z_n=k|x_n;model)
function R = expectation(X,model)
K = numel(model.w);
[D,N] = size(X);
logR = zeros(K,N);
get_diag = logical(eye(D));
% Compute non-normalized posterior probabilities
for k=1:K
    % Covariance matrix
    Sigma_k = model.A(:,k)*model.A(:,k)' + model.sigma2(k)*eye(D);
    % Ensures positive-definiteness :
    Sigma_k(get_diag(:)) = abs(Sigma_k(get_diag(:)));
    logR(k,:) = log(model.w(k)) + loggausspdf(X,Sigma_k); % 1*N
end
% Normalize the posterior probabilities
logR = bsxfun(@minus, logR, logsumexp(logR,1)); % K*N
R = exp(logR); % K*N
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAXIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT:
% X (D*N) : N input D-dimensional complex vectors
% R (K*N) : log-posteriors: log p(z_n=k|x_n;model)
%%% OUTPUT:
% model (struct) : Updated parameters
%    .A (D*K) : each column a_k is source k's steering vector
%    .sigma2 (K*1) : noise variance for each source (>0)
%    .w (K*1) : cluster weights ( sum(w) = 1 )
function model = maximization(X,R)
K = size(R,1);
D = size(X,1);
model.sigma2 = zeros(K,1);
model.A = zeros(D,K);
% Update cluster weights:
model.w = mean(R,2); % K*1
for k=1:K
    % Weighted data
    r_k = sqrt(R(k,:)/sum(R(k,:))); % 1*N
    X_k = bsxfun(@times,X,r_k); % D*N
    
    % Weighted sample covariance matrix
    Gamma_k = X_k*X_k'; % D*D 
    
    % Leading (eigenvector,eigenvalue) pair of Gamma_k
    [v_k,lambda2_k,~] = svd(Gamma_k); % D*1, double
    v_k = v_k(:,1);
    lambda2_k = lambda2_k(1,1);
    
    % Update noise variance
    model.sigma2(k) = (trace(Gamma_k)-lambda2_k)/(D-1);
    
    % Update steering vector
    model.A(:,k) = sqrt(lambda2_k-model.sigma2(k))*v_k;
end
% Regularize noise variance to avoid numerical problems:
model.sigma2(model.sigma2<1e-6) = 1e-6;
end

%%%%%%%%%%%% logarithm of 0-mean complex-circular Gaussian pdf %%%%%%%%%%%%
%%% INPUT:
% X (D*N) : N input D-dimensional complex vectors
% Sigma (D*D) : Covariance of the 0-mean complex-circular Gaussian
%%% OUTPUT:
% P (1*N) : probability density values for each of the N datapoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logP = loggausspdf(X, Sigma)

%%%%%%%%%%%%%%% go to real for certainty
X = [real(X); imag(X)];
Sigma = [real(Sigma) -imag(Sigma); imag(Sigma) real(Sigma)];
%%%%%%%%%%%%%%%

D = size(X,1);
[U,p]= chol(Sigma); % U is (D*D)
if p ~= 0
    error('    ERROR: Sigma is not positive definite.');
end
Q = U'\X; % Weighted data (D*N)
q = -dot(Q,Q,1);  % (1*N) Quadratic term
c = -D*log(pi)-sum(log(diag(U))); % (1*N) normalization constant
logP = (c + q)/2; % (1*N)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% logsumexp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
%   By default dim = 1 (columns).
% Written by Michael Chen (sth4nth@gmail.com).
function s = logsumexp(x, dim)
if nargin == 1,
    % Determine which dimension sum will use
    dim = find(size(x)~=1,1);
    if isempty(dim), dim = 1; end
end

% subtract the largest in each column
y = max(x,[],dim);
x = bsxfun(@minus,x,y);
s = y + log(sum(exp(x),dim));
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end
end
