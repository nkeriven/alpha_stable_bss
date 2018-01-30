function mask = oracle_masks(S_img)
%%%%%%%%%% Estimate binary masks from an input multichannel spectrogram %%%
%%%%%%%%%% using a simple Gaussian mixture model at each frequency %%%%%%%%
%%% INPUT:
% S_img (F*T*M*K) : true source images
%%% OUTPUT:
% mask (F*T*K) : Oracle binary masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[F,T,~,K] = size(S_img);

% Source energy in at each (f;t) bin:
E = reshape(sum(abs(S_img).^2,3),[F,T,K]);

% Compute masks:
maxE = max(E,[],3); %F*T
mask = bsxfun(@eq,E,maxE); %F*T*K
end