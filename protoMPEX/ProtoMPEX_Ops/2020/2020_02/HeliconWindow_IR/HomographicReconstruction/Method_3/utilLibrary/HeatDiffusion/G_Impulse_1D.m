function [g] = G_Impulse_1D(x,t,Ns)
%==========================================================================
% 1D-TRANSIENT HEAT CONDUCTION SOLUTION WITH IMPULSE SURFACE HEAT FLUX
% Independent variables: x and t
%==========================================================================
% G_Impulse_1D is the solution of the 1D transient heat diffusion equation subject
% to an IMPULSE surface heat flux at x = 0. It can be viewed as the Green's
% function solution to the problem:

% u_t = u_xx
%
% u(x,t=0) = 0
% ux @ (x = 0) = -delta(t)
% ux @ (x = 1) = 0

% Where u is the dimensionless temperature u = (T-T0)/(q0*L/k), t is the
% dimensionless time (alpha*t/L^2) and x is the dimensionless time (x/L)
% subscripts t and xx represent first temporal and second spatial derivatives.

% x: axial location normalized to L, the axial length of slab, R^1
% t: dimensionless time value, R^1
% Ns: total number of terms in the infinite sum

if (t ==0 && x == 1)
    g  = 0;
else
    C = exp(-(([1:Ns]'*pi).^2)*t);
    E = cos([1:Ns]'*pi*x);
    g  = 1 + 2*C'*E;
end

end

