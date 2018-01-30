function X = stft(x,fs,time_win,overlap_percent)
%%%%%%%%%%%%%%% Compute a spectrogram from a temporal signal %%%%%%%%%%%%%%
%%% INPUT:
% x (n*1) : real temporal signal
% fs (real) : frequency of sampling (Hz)
% time_win (real) : time window (ms)
% overlap_percent (real) : overlap between two windows (%)
%%% OUTPUT:
% X (F*T) : Complex sectrogram (only positive frequencies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=x(:); % make sure it's a column

% Spectrogram parameters:
win_shift = time_win*(100-overlap_percent)/100;
ns=fs/1000*time_win; % Number of samples in each time window
noverlap=fs/1000*(time_win-win_shift);
nfft=2^nextpow2(ns);
nzp=ns/2; % Number of zero padding
zp=zeros(nzp,1);
window=hamming(ns); % Hamming window
x=[zp;x;zp]; % zero padded signal

M = length(window);
if (M<2)
   error('myspectrogram: Expect complete window, not just its length'); 
end;
if (M<2)
   error('myspectrogram: Expect complete window, not just its length'); 
end;
if length(x)<M % zero-pad to fill a window:
  x = [x;zeros(M-length(x),1)]; 
end;
Modd = mod(M,2); % 0 if M even, 1 if odd
Mo2 = (M-Modd)/2;
w = window(:); % Make sure it's a column

nhop = M-noverlap;

nx = length(x);
nframes = 1+ceil(nx/nhop);

X = zeros(nfft,nframes); % allocate output spectrogram

zp = zeros(nfft-M,1); % zero-padding for each FFT
xframe = zeros(M,1);
xoff = 0 - Mo2; % input time offset = half a frame
for m=1:nframes
%  M,Mo2,xoff,nhop
  if xoff<0
    xframe(1:xoff+M) = x(1:xoff+M); % partial input data frame
  else
    if xoff+M > nx
      xframe = [x(xoff+1:nx);zeros(xoff+M-nx,1)];
    else
      xframe = x(xoff+1:xoff+M); % input data frame
    end
  end
  xw = w .* xframe; % Apply window
  xwzp = [xw(Mo2+1:M);zp;xw(1:Mo2)];
  X(:,m) = fft(xwzp);
  xoff = xoff + nhop; % advance input offset by hop size
end

%%%%% Keep only the inside part of the spectrogram %%%%%
[nfreq,nframe]=size(X);
inside=(nzp/(ns-noverlap)+1):(nframe-nzp/(ns-noverlap));
X = X(1:(nfreq/2+1),inside);

end
