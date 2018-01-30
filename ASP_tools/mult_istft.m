function x = mult_istft(X,fs,time_win,overlap_percent,sig_length)
%%%%%%%%%%%%%%% Compute a spectrogram from a temporal signal %%%%%%%%%%%%%%
%%% INPUT:
% X (F*T*...) : multichannel complex sectrogram (only positive frequencies)
% fs (real) : frequency of sampling (Hz)
% time_win (real) : time window (ms)
% overlap_percent (real) : overlap between two windows (%)
% sig_length (int) : length of the original signal
%%% OUTPUT:
% x (n*...) : array whose 1st dimension contains the time-domain signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_X=size(X);
M = prod(size_X(3:end)); % Number of channels
xx = zeros(sig_length,M);

XX = reshape(X,[size_X(1:2),M]);

for m = 1:M
    xx(:,m) = istft(XX(:,:,m),fs,time_win,overlap_percent,sig_length);
end
x = reshape(xx,[sig_length,size_X(3:end)]);

end