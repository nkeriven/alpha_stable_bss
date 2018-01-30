function source=load_sources(signal_path,fs,start_time,end_time,nmics) 
%%%%%%%% Load multiple sounds in an array with unified length and frequency
%%%%%%%% of sampling
%%% INPUT:
%  - signal_path: cell containing signal paths or string for a single
%    source
%  - fs: desired frequency of sampling (Hz)
%  - start_time: start time in seconds ([]=begining)
%  - end_time: end time in seconds ([]=end)
%  - nmics: number of channels ([]=all)
%%% OUTPUT:
%  - source (L,nmics,nsrc): output array containing the standardized
%                            source signals
%%% Author:
%    Antoine Deleforge (Sep. 2017)
%    antoine.deleforge@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isempty(nmics))
    nmics = 1;
end
if(~iscell(signal_path))
    signal_path={signal_path};
end

nsrc = numel(signal_path);
source_cell = cell(nsrc,1);
nsamples=0;

for k=1:nsrc
    [x, Fs] = audioread(signal_path{k}); % Time signal
    
    if(~isempty(start_time))
        ss=round(start_time*Fs)+1;
    else
        ss=1;       
    end
    if(~isempty(end_time))
        es=round(end_time*Fs);
    else
        es=size(x,1);       
    end    
    L = numel(resample(x(ss:es,1), fs, Fs));
    source_cell{k} = zeros(L,nmics);
    for i=1:nmics
        source_cell{k}(:,i) = resample(x(ss:es,i), fs, Fs);
    end
    nsamples = max(nsamples,L);  
end

source = zeros(L,nmics,nsrc);

for k = 1:nsrc
    Ls = size(source_cell{k},1);
    source(1:Ls,:,k) = source_cell{k};
end

end

