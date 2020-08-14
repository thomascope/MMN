function [F,x2,F1,F2] = spm_Ncdf_jdwTA(x,VR1,Vx1,VR2,Vx2)
% Cumulative Distribution Function (CDF) for univariate Normal distributions: J.D.  Williams aproximation
% FORMAT F = spm_Ncdf_jdw(x,u,v)
%
% x - ordinates
% u - mean              [Defaults to 0] % TA: VR
% v - variance  (v>0)   [Defaults to 1] % TA: Vx
% F - pdf of N(u,v) at x (Lower tail probability)
%__________________________________________________________________________
%
% spm_Ncdf implements the Cumulative Distribution Function (CDF) for
% the Normal (Gaussian) family of distributions.
%
% References:
%--------------------------------------------------------------------------
% An Approximation to the Probability Integral
% J. D. Williams 
% The Annals of Mathematical Statistics, Vol. 17, No. 3. (Sep., 1946), pp.
% 363-365. 
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_Ncdf_jdw.m 4836 2012-08-10 15:55:21Z karl $


%-Approximate integral
%--------------------------------------------------------------------------
x2    = -(x - VR1)./sqrt(abs(Vx1));
F1    = sqrt(1 - exp(-(2/pi)*x2.^2))/2;
i    = x2 < 0;
F2 = F1;
F2(i) = -F1(i);
Fa    = F2 + 1/2;
Fa(x>=VR2) = 0;

x2    = (x - VR2)./sqrt(abs(Vx2));
F1    = sqrt(1 - exp(-(2/pi)*x2.^2))/2;
i    = x2 < 0;
F2 = F1;
F2(i) = -F1(i);
Fb    = F2 + 1/2;
Fb(x<VR2) = 0;
Fb = Fb.*.65;

F = Fa+Fb;



