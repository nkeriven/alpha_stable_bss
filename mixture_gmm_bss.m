%%%%%%%% class for rank-1 gmm for BSS %%%%%%
% Inherits of the mixture class
% To be used with the SketchMLbox: http://sketchml.gforge.inria.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author's Contact:
%    Nicolas Keriven
%   nicolas.keriven@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



classdef mixture_gmm_bss<mixture
    
    properties
        a; % steering vectors (d*k)
        sig; % noise levels (1*k)
    end
    
    methods
        function self=mixture_gmm_bss(a,sig,weights)
            self.a=a;
            self.sig=sig;
            self.weights=weights(:);
            [self.d,self.k]=size(self.a);
            
        end
        
        %% Methods required by the parent class
        
        function []=shift(self,s)
            % not used here
        end
        
        function []=normalize(self,m)
            % not used here
        end
        
        %% additional (allow for the use of the already existing mixture_gmm class)
        function mix=to_gmm(self)
            Sigma=zeros(self.d,self.d,self.k);
            for l=1:self.k
                a1=self.a(1:end/2,l);
                a2=self.a(end/2+1:end,l);
                A=a1*a1' + a2*a2';
                B=a1*a2' - a2*a1';
                Sigma(:,:,l) = [A B;-B A] + self.sig(l)*eye(self.d); 
            end
            mix=mixture_gmm(zeros(self.d,self.k),Sigma,self.weights);
        end
        
    end
end

