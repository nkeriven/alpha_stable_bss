function [mix,source_img,delay,gain] = random_mix(source,nmics,...
                                              max_delay,max_gain,min_delta)
%%%%% Mix monochannel sources with random gains and delays %%%%%%%%%%%%%%%%
%%% INPUT:
%  - source (L,nsrc) : input L-sample monochannel source signals
%  - nmics (int): desired nubmer of microphones
%  - max_delay (int): maximum delay with the reference microphone
%  - max_gain (dB): maximum absolute gain wrt ref. microphone in dB
%  - min_delta (int): minimum difference between any delay pair (def. 0)
%%% OUTPUT:
%  - mix (L,nmics): mixed nmics-channel signal
%  - source_img (L,nmics,nsrc): array containing the source images
%  - delay (nmics,nsrc) : array of delays
%  - gain (nmics,nrsc) : array of gains
%%% Author:
%    Antoine Deleforge (Sep. 2017)
%    antoine.deleforge@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~exist('min_delta','var') || isempty(min_delta))
    min_delta =0;
end
[L,nsrc] = size(source);

if(min_delta==0)
    delay = randi(2*max_delay+1,[nmics-1,nsrc])-max_delay-1;
    % Delays of the reference microphone are 0
    delay = [zeros(1,nsrc);delay]; 
else
    % Generate delays between -max_delay and max_delay spaced by at least
    % min_delta samples
    tau = floor(max_delay/min_delta);
    if(tau<(nmics-1)*nsrc)
        error('Cannot find enough distinct delays to perform the mix');
    end
    all_delay = [-tau:-1,1:tau]*min_delta; % Except 0
    delay = all_delay(randperm(2*tau));
    delay = reshape(delay(1:(nmics-1)*nsrc),[nmics-1,nsrc]); 
    % Delays of the reference microphone are 0
    delay = [zeros(1,nsrc);delay]; % nmics*nsrc
end

% Generate gains in [-max_gain,max_gain] (dB)
gain = 2*max_gain*(rand(nmics-1,nsrc)-0.5); % (dB) (nmics-1)*nsrc
gain = 10.^(gain/20); % Gain factor M*K
gain = [ones(1,nsrc);gain]; % Reference microphone gains are 1

% Apply gains and delays (time domain)
source_img = zeros(L,nmics,nsrc);
for m=1:nmics
    for k=1:nsrc
        % Start and end samples in the source signal:
        src_start = max(1,delay(m,k)+1);
        src_end = min(L,L+delay(m,k));
        % End sample in the image signal:
        img_start = max(1,-delay(m,k)+1);
        img_end = min(L,L-delay(m,k));
        source_img(img_start:img_end,m,k) = ...
            gain(m,k)*source(src_start:src_end,k);
    end
end

% Mix sources
mix = sum(source_img,3);

end