function [A_opt, B_opt, C_opt, w_opt] = LS4p_parab(x, w0, dw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB implementation of four-parameter Least Squares sine        %%%
%%% fitting method proposed in "A Computationally Efficient           %%%
%%% Non-iterative Four-parameter Sine Fitting Method, IET Signal      %%%
%%% Processing, 2021                                                  %%%
%%%                                                                   %%%
%%% Input parameters:                                                 %%%
%%%              - x: Data samples                                    %%%
%%%              - w0: initial relative angular frequency estimate    %%%
%%%              - dw: deviation from w0, where the three-parameter   %%%
%%%                    LS cost function is evaluated                  %%%
%%%                                                                   %%%
%%% Output parameters: Optimal sine parameters of the fitting         %%%
%%%                                                                   %%%
%%% Written by: Balázs Renczes                                        %%%
%%%                                                                   %%%
%%% Last modified: April 8, 2021                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = size(x);
if (s(2) > s(1))     % We ensure that x is a column vector
    x = x.';
end
N = length(x);
n = (1:N).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A three-parameter least squares fitting is perfomred at w0, which %%%
%%% is usually the interpolated FFT estimator.                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A_ipfft, B_ipfft, C_ipfft] = LS3p(w0, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The fitting is performed with the obtained parameter set. Then    %%%
%%% the error vector and the cost function is calculated.             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_ipfft = C_ipfft + A_ipfft*cos(w0*n) + B_ipfft*sin(w0*n);
e_ipfft = x-y_ipfft;
CF_ipfft = e_ipfft.'*e_ipfft;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The procedure is repeated at w0-dw and at w0+dw                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A1, B1, C1] = LS3p(w0-dw, x);
y1 = C1 + A1*cos((w0 - dw)*n) + B1*sin((w0 - dw)*n);
e1 = x-y1;
CF1 = e1.'*e1;
[A2, B2, C2] = LS3p(w0+dw, x);
y2 = C2 + A2*cos((w0 + dw)*n) + B2*sin((w0 + dw)*n);
e2 = x-y2;
CF2 = e2.'*e2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The optimal frequency is calculated, see Equations (20) and (21)  %%%
%%% in the paper.                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = [1/2 -1 1/2; -1/2 0 1/2; 0 1 0;];
p_parab = M*[CF1; CF_ipfft; CF2;];
dw_opt = -p_parab(2)/(2*p_parab(1))*dw;
w_opt = w0 + dw_opt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The optimal parameters of A, B and C are calculated at the        %%%
%%% optimal frequency.                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A_opt, B_opt, C_opt] = LS3p(w_opt, x);
y_opt = C_opt + A_opt*cos(n*w_opt) + B_opt*sin(n*w_opt);
e_opt = x - y_opt;
CF_opt = e_opt.'*e_opt;
end

