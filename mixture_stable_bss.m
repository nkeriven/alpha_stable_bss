%%%%%%%% class for rank-1 alpha_stable mixture %%%%%%
% Inherits of the mixture class
% To be used with the SketchMLbox: http://sketchml.gforge.inria.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Nicolas Keriven
%   nicolas.keriven@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



classdef mixture_stable_bss<mixture
    
    properties
        a; % steering vectors (d*k)
        sig; % noise levels (1*k)
        alpha; % alpha coeff (1*k)
    end
    
    methods
        %Constructor
        function self=mixture_stable_bss(a,sig,alpha,weights)
            
            self.a=a;
            self.sig=sig;
            self.alpha=alpha(:)';
            self.weights=weights(:);
            [self.d,self.k]=size(self.a);
            
        end
        
        %% Methods required by parent class
        
        function []=shift(self,s)
            % no used here
        end
        
        function []=normalize(self,m)
            % not used here
        end
        
        %% additional (for display using existing class)
        function mix=to_elliptic_stable_for_display(self)
            Sigma=zeros(self.d,self.d,self.k);
            for l=1:self.k
                a1=self.a(1:end/2,l);
                a2=self.a(end/2+1:end,l);
                A=a1*a1' + a2*a2';
                B=a1*a2' - a2*a1';
                Sigma(:,:,l) = [A B;-B A] + self.sig(l)*eye(self.d); % the sigma is technically wrong, for display here
            end
            mix=mixture_elliptic_stable(zeros(self.d,self.k),Sigma,self.alpha,self.weights);
        end
        
        
    end
end

