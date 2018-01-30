%%%%%%%% class for rank-1 alpha stable mixture estimation for BSS %%%%%%
% Inherits of the mixture_estimator class
% To be used with the SketchMLbox: http://sketchml.gforge.inria.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Nicolas Keriven
%   nicolas.keriven@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef stable_bss_estimator<mixture_estimator
    
    properties
        W2;
        Wtilde;
    end
    
    methods
        %% build
        % sk_wrapper : see sketch_wrapper.m (in SketchMLbox)
        % varargin : options, can be passed as couple 'name','value'
        %   see load_estimator_options.m for available options and default
        %   values
        function self=stable_bss_estimator(sk_wrapper,varargin)
            
            self@mixture_estimator(sk_wrapper,varargin{:});
            self.p=self.d+2; % params dim: [a; sig; alpha]
            self.W2=sum((self.sk_wrapper.W).^2,1);
            self.Wtilde = [self.sk_wrapper.W(end/2+1:end,:); -self.sk_wrapper.W(1:end/2,:)];
            
            % possible constraints bounds to avoid stability problem
            self.lb=-Inf*ones(self.p,1);
            self.ub = Inf*ones(self.p,1);
            
            self.lb(end-1)=self.options.min_var; % lower bound on noise variance: avoid numerical problem
            self.lb(end)=self.options.min_alpha; % bound alpha
            
            self.ub(end)=2;
            
            if self.options.constrained
                self.lb(1:self.d)=-sqrt(2*self.sk_wrapper.rad); % lower bound on steering vector a
                
                self.ub(1:self.d)=sqrt(2*self.sk_wrapper.rad); % upper bound on steering vector a
                self.ub(end-1)=2*self.sk_wrapper.rad; % upper bound on noise variance sig
            end
            
        end
        
        
        %% Required methods
        
        function mix=construct_mixture(self,params,weights)
            a=params(1:self.d,:);
            sig=params(end-1,:);
            alpha=params(end,:);
            mix=mixture_stable_bss(a,sig,alpha,weights);
        end
        
        function [params,weights]=toparams(self,mix)
            
            weights=mix.weights;
            params=zeros(self.p,length(weights));
            params(1:self.d,:)=mix.a;
            params(end-1,:)=mix.sig;
            params(end,:)=mix.alpha';
        end
        
        function v=init_param(self,params_curr,X)
            % current parameters and eventual data not used here
            v(1:self.d)=randn(self.d,1);% blind guess for initialization
            v(1:self.d)=sqrt(self.sk_wrapper.mean_var)*v(1:self.d)/norm(v(1:self.d));
            v(self.d+1)=10^(-1);
            v(self.d+2)=1.95;
        end
        
        % cost function and gradient
        function [phi,jphi]=sketch_distrib(self,param)
            W=self.sk_wrapper.W;
            d=self.d;
            m=self.m;
            alpha=param(end);
            Wa=W'*param(1:d);
            Wta = self.Wtilde'*param(1:d);
            s=(1/2)*(Wa.^2 + Wta.^2);
            
            phi = [exp(-s.^(alpha/2)).*exp(-param(d+1)*self.W2'/2)...
                .*self.sk_wrapper.freq_weights; zeros(self.m,1)]; % the characteristic function is real
            
            if nargout>1
                phi1=phi(1:m);
                jphi=@(x)(self.intermediate_jphi(x,phi1,W,self.W2,self.Wtilde,Wa,Wta,s,alpha,m));
            end
        end
        
        function c=intermediate_jphi(self,x,phi1,W,W2,Wtilde,Wa,Wta,s,alpha,m)
            b=phi1.*x(1:m);
            u=s.^(alpha/2-1).*b;
            c=[ -(alpha/2)*(W*(Wa.*u) + Wtilde*(Wta.*u)); % a
                -W2*b/2; % sig
                -(1/2)*(log(s).*s.^(alpha/2))'*b ]; % alpha
        end
        
        
    end
end


