function [model] = sk_bss(X,K,varargin)
%%%%%%%% Sketch-algorithm for complex rank-1 Alpha-stable Mixture Estimation %%%%%%
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

%%%%%%%%%% parsing Parameters
p=inputParser();
p.KeepUnmatched = true;

p.addParameter('sk_size', 300, @isnumeric);
p.addParameter('type_estimator', 'alpha', @(x)true);

p.parse(varargin{:});
options=p.Results;

[M,T] = size(X);

%%% unfold complex data
Xreal = [real(X); imag(X)];

%%% construct the sketch
skw = sketch_wrapper(Xreal, options.sk_size, varargin{:});
skw.set_sketch(Xreal);

switch options.type_estimator
    
    case 'alpha'
        %%% estimate model
        estimator = stable_bss_estimator(skw, varargin{:});
        estimator.estim(K,Xreal);
        alphas = estimator.mixture.alpha;
        
    case 'sk_gmm'
        estimator = gmm_bss_estimator(skw, varargin{:});
        estimator.estim(K,Xreal);
        alphas = 2*ones(1,K);
        
end

%%% store model
A = estimator.mixture.a;
model.A = A(1:M,:) + 1i*A(M+1:2*M,:);
model.sigma2 = estimator.mixture.sig;
model.w = estimator.mixture.weights;
model.alpha = alphas;
end
