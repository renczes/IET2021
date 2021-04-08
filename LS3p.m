function [A, B, C] = LS3p(w, x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB implementation of three Least Squares sine fitting method  %%%
%%%                                                                   %%%
%%% Supplementary file for sfit_parab.m                               %%%
%%%                                                                   %%%
%%% Input parameters:
%%%              - x: Data samples
%%%              - w: relative angular frequency                      %%%
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

D0 = [cos(w*n) sin(w*n) ones(N,1)];

p = inv*(D0'*D0) * (D0'*x);

A = p(1);
B = p(2);
C = p(3);
end
