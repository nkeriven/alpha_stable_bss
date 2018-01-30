function [mask,A,alpha,noise] = estimate_masks_sawada(X,K,nb_rep,verb)
%%%%%%%%%% Estimate binary masks from an input multichannel spectrogram %%%
%%%%%%%%%% using a simple Gaussian mixture model at each frequency %%%%%%%%
%%% INPUT:
% X (F*T*M): Input multichanel spectrogram
% K (int) : number of sought sources
% verb (int) : verbosity
%%% OUTPUT:
% mask (F*T*K) : F independent binary masks
% A (F*M*K) : Estimated mixing matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~exist('verb','var'))
    verb = 0;
end
if(~exist('nb_rep','var'))
    nb_rep = 1;
end
[F,T,M] = size(X);
mask = zeros(F,T,K);
A = zeros(F,M,K);
if(verb>0);fprintf(1,'Mask Estimation:\n');end;
for f=1:F
    if(verb>0);fprintf(1,'  Frequency %d/%d\n',f,F);end;
    Xf = reshape(X(f,:,:),[T,M]).'; % M*T GMM data
    best_loglike = -Inf;
    for i=1:nb_rep
        [label,model,R] = sawada_em(Xf,K,verb-1);
        if R>best_loglike
            best_model = model;
            best_label = label;
        end
    end
    model = best_model;
    label = best_label;
    % Normalize steering vectors:
    A(f,:,:) = bsxfun(@rdivide,model.A,model.A(1,:)); % M*K
    alpha(f,:) = 2*ones(1,K);
    noise(f,:) = model.sigma2;
    for k=1:K
        mask(f,:,k) = (label==k);
    end
end
end