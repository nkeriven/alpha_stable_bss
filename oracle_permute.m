function [mask_permute,A_permute,alpha_permute,noise_permute] = oracle_permute(mask,X,S_img,A,alpha,noise)
%%%%%%%%% Optimally permute masks given the true source images %%%%%%%%%%%%
%%% INPUT:
% mask (F*T*K) : input masks
% X (F*T*M): Observed multichanel spectrogram
% S_img (F*T*M*K) : true source images
% A (F*M*K) : mixing matrices (optional)
%%% OUTPUT:
% mask_permute (F*T*K) : optimally permuted masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[F,T,~,K] = size(S_img);
mask_permute = zeros(F,T,K);
if(exist('A','var'))
    A_permute = A;
else
    A_permute = NaN;
end

% Matrix whose rows are all permtations of (1:K)
P = perms(1:K); % nperms * K
nperms = size(P,1);
error_perm = zeros(nperms,1);
for f=1:F
    % Compute source estimation error for all possible permutations
    for i = 1:nperms;
        % Source estimates given current permutation
        S_mask = bsxfun(@times,X(f,:,:),...
                             reshape(mask(f,:,P(i,:)),[1,T,1,K])); %1*T*M*K
        diff = S_mask - S_img(f,:,:,:);
        error_perm(i) = norm(diff(:));
    end
    [~,ibest] = min(error_perm);
    % Select best permutation at frequency f
    mask_permute(f,:,:) = mask(f,:,P(ibest,:));
    
    % Permute mixing matrices
    if(exist('A','var'))
        A_permute(f,:,:) = A(f,:,P(ibest,:));
    end
    if(exist('alpha','var'))
        alpha_permute(f,:) = alpha(f,P(ibest,:));
    end
    if(exist('alpha','var'))
        noise_permute(f,:) = noise(f,P(ibest,:));
    end
end
end