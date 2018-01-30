function A = get_mixing_matrices(S_img)
%%%%%%%%% Optimally permute masks given the true source images %%%%%%%%%%%%
%%% INPUT:
% S_img (F*T*M*K) : Source images
%%% OUTPUT:
% A (F*M*K) : Mixing matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[F,~,M,K] = size(S_img);
A = zeros(F,M,K);

for f=1:F
    % Divide by ref microphone
    for k=1:K
        % Ignore points where source is inactive
        active = S_img(f,:,1,k)>1e-08;
        Sd_fk = bsxfun(@rdivide,S_img(f,active,:,k),...
                                S_img(f,active,1,k)); % 1*Tactive*M*1
        A(f,:,k) = reshape(median(Sd_fk,2),[M,1]);
    end
end

end