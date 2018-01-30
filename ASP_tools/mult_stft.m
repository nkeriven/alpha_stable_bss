function X = mult_stft(x,fs,time_win,overlap_percent)
%%%%%%%%%%%%%%% Compute a spectrogram from a temporal signal %%%%%%%%%%%%%%
%%% INPUT:
% x (n*...) : array whose 1st dimension contains the time-domain signal
% fs (real) : frequency of sampling (Hz)
% time_win (real) : time window (ms)
% overlap_percent (real) : overlap between two windows (%)
%%% OUTPUT:
% X (F*T*...) : multichannel complex sectrogram (only positive frequencies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_x = size(x);
M = prod(size_x(2:end)); % Number of channels
xx = reshape(x,[size_x(1),M]);

XX1 = stft(xx(:,1),fs,time_win,overlap_percent);
FT = size(XX1);
size_XX = [FT,M];
XX = zeros(size_XX);
for m = 1:M
    XX(:,:,m) = stft(xx(:,m),fs,time_win,overlap_percent);
end
X = reshape(XX,[FT,size_x(2:end)]);

end