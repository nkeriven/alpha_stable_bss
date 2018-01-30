function A = theoretical_mixing_matrices(delay,gain,F)
%%%%%%%%% Optimally permute masks given the true source images %%%%%%%%%%%%
%%% INPUT:
% delay (M*K) : Time delays in samples
% F (int) : number of positive frequency bins considered
%%% OUTPUT:
% A (F*M*K) : Mixing matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M,K] = size(delay);
A = zeros(F,M,K);

for f=1:F
    A(f,:,:) = gain.*exp(1i*(f-1)*delay*pi/F);
end

% Set first microphone to reference:
A = bsxfun(@rdivide,A,A(:,1,:)); % F*M*K

end