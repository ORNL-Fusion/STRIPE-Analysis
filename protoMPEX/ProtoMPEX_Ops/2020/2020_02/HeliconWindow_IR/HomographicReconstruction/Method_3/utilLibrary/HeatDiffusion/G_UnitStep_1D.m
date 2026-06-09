function [g] = G_UnitStep_1D(x,t,Ns)
%==========================================================================
% 1D-TRANSIENT HEAT CONDUCTION SOLUTION WITH UNIT STEP SURFACE HEAT FLUX
% Independent variables: x and t
%==========================================================================
% G_UnitStep_1D is the solution of the 1D transient heat diffusion equation subject
% to a UNIT STEP surface heat flux at x = 0. It can be viewed as the Green's
% function solution to the problem:

% u_t = u_xx
%
% u(x,t=0) = 0
% ux @ (x = 0) = -Heavyside(t)
% ux @ (x = 1) = 0

% where u is the dimensionless temperature u = (T-T0)/(q0*L/k), t is the
% dimensionless time (alpha*t/L^2) and x is the dimensionless time (x/L)
% subscripts t and xx represent first temporal and second spatial derivatives.

% x: axial location normalized to L, the axial length of slab, R^1
% t: dimensionless time value, R^1
% Ns: total number of terms in the infinite sum

C = 1 - exp(-(([1:Ns]*pi).^2)*t);
E = cos([1:Ns]*pi*x)./(([1:Ns]*pi).^2);
g  = t + 2*sum(C.*E);

end

