function [mask,A,noise,alpha,best_loglike] = estimate_masks(X,K,varargin)
%%%%%%%%%% Estimate binary masks from an input multichannel spectrogram %%%
%%%%%%%%%% using a simple Gaussian mixture model at each frequency %%%%%%%%
%%% INPUT:
% X (F*T*M): Input multichanel spectrogram
% K (int) : number of sought sources
% options ('name', 'value'):
%       'type_estimator': 'em_gmm', 'sk_gmm', 'sawada', 'alpha'
%       'nbr_replicates': number of random initializations
%       other options for the sketch method, ex: 'sk_size' (see sketchmlbox
%           for other options)
%%% OUTPUT:
% mask (F*T*K) : F independent binary masks
% A (F*M*K) : Estimated mixing matrices
% alpha (F*K) : Estimated alpha
% noise (F*K) : Estimated noise amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Antoine Deleforge
%    deleforge.antoine@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=inputParser();
p.KeepUnmatched=true;
p.addParameter('nbr_replicates',1,@isnumeric);
p.addParameter('type_estimator','em_gmm',@(x)true);
p.addParameter('verb',1,@isnumerical);

p.parse(varargin{:});
options=p.Results;

[F,T,M] = size(X);
mask = zeros(F,T,K);
A = zeros(F,M,K);
noise = zeros(F,K);
alpha = zeros(F,K);
if(options.verb>0);fprintf(1,'Mask Estimation:\n');end;
for f=1:F
    if(options.verb>0);fprintf(1,'  Frequency %d/%d\n',f,F);end;
    Xf = reshape(X(f,:,:),[T,M]).'; % M*T GMM data
    
    best_loglike(f) = -Inf;
    for it=1:options.nbr_replicates
        switch options.type_estimator
            case 'sawada'
                model = sawada_em(Xf,K);
            case 'em_gmm'
                model = gmm_bss_em(Xf,K,varargin{:});
            otherwise
                model = sk_bss(Xf,K,varargin{:});
                
        end
        emmix = mixture_gmm_bss([real(model.A);imag(model.A)],model.sigma2,model.w);
        [ll,label] = loglike(emmix.to_gmm,[real(Xf);imag(Xf)]);
        if ll>best_loglike(f)
            best_model = model;
            best_labels = label;
            best_loglike(f) = ll;
        end
    end
    
    model = best_model;
    label = best_labels;
    A(f,:,:) =  bsxfun(@rdivide,model.A,model.A(1,:));
    alpha(f,:) = model.alpha;
    noise(f,:) = model.sigma2;
    
    for k=1:K
        mask(f,:,k) = (label==k);
    end
end
end