function x = istft(X,fs,time_win,overlap_percent,sig_length,alpha)
%%%%%%%%%%%%%% Inverse a spectrogram to get a temporal signal %%%%%%%%%%%%%
%%% INPUT:
% X (F x T) : complex spectrogram (F positive frequency bins, T time bins)
% fs (real) : frequency of sampling (Hz)
% time_win (real) : time window (ms)
% overlap_percent (real) : overlap between two windows (%)
% sig_length (int) : length of the original signal
% alpha (real,optional) : multiplicator to make stft o istft identity
%%% OUTPUT:
% x (sig_length*1) : real signal (only positive frequencies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse spectrogram parameters:
win_shift = time_win*(100-overlap_percent)/100;
F=size(X,1);
b = [X;conj(X(F-1:-1:2,:))];

ns=fs/1000*time_win;
no=fs/1000*(time_win-win_shift);
hop = ns-no;

[nfft,nframes] = size(b);

No2 = nfft/2; % nfft assumed even
a = zeros(1, nfft+(nframes-1)*hop);
xoff = 0 - No2; % output time offset = half of FFT size
for col = 1:nframes
  fftframe = b(:,col);
  xzp = ifft(fftframe);
  % xzp = real(xzp); % if signal known to be real
  x = [xzp(nfft-No2+1:nfft); xzp(1:No2)];
  if xoff<0 % FFT's "negative-time indices" are out of range
    ix = 1:xoff+nfft;
    a(ix) = a(ix) + x(1-xoff:nfft)'; % partial frames out
  else
    ix = xoff+1:xoff+nfft;
    a(ix) = a(ix) + x';  % overlap-add reconstruction
  end
  xoff = xoff + hop;
end

if(~exist('alpha','var'))
    % Compute multiplicator so that STFT o invSTFT = identity
    xrand=rand(time_win*fs*10,1);
    xrand2= istft(stft(xrand,fs,time_win,overlap_percent),fs,time_win,overlap_percent,time_win*fs*10,1);
    alpha = std(xrand)./std(xrand2);
end

x = a(1:sig_length).*alpha;

end