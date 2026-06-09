function [ff_IJ,A_MN,Phi_IM,Phi_JN] = FourierCosine_2D_v2(x_I,y_J,f_IJ,M,N)
% FOURIERCOSINE_2D_V2 takes a 2D array F_IJ and computes the Fourier cosine
% coefficients. The domain of X_I and Y_J are [0,1]
% In this version, we change the convention. We associate M with x, N with
% y.
%
% INPUTS:
% X_I: Column vector of "x" values with "I" elements, [0,1]
% Y_J: Column vector of "y" values with "J" elements, [0,1]
% F_IJ: 2D matrix of function to Fouerier decompose
% M: Number of terms in the p(x) eigenfunction sum
% N: Number of terms in the p(y) eigenfunction sum
% OUTPUTS:
% A_MN: 2D array of Fourier Cosine Coefficients, MxN elements
% Phi_IM: 2D array of "cos(mpix)" eigenfunctions for all "m" and all "x"
% Phi_JN: 2D array of "cos(npiy)" eigenfunctions for all "n" and all "y"
% ff_IJ : Approximated version of f_IJ

% J.F. Caneses Marin
% Created 2019_03_13

% Check size of independent variables
% -------------------------------------------------------------------------
I = length(x_I);
J = length(y_J);

dx = mean(diff(x_I));
dy = mean(diff(y_J));

% Check dimensions of inputs:
% -------------------------------------------------------------------------
if ~iscolumn(x_I)
    x_I = x_I';
end
if ~iscolumn(y_J)
    y_J = y_J';
end

if size(f_IJ) ~= [I,J]
   f_IJ = f_IJ';
end

% Define "m" and "n" numbers for the infinite sum
% -------------------------------------------------------------------------
m_M = [0:(M-1)]';
n_N = [0:(N-1)]';

% Define Kronecker delta functions
% -------------------------------------------------------------------------
d_mm = 0.5*ones(size(m_M))'; d_mm(1) = 1;
d_nn = 0.5*ones(size(n_N))'; d_nn(1) = 1;
D_MM = diag(1./d_mm);
D_NN = diag(1./d_nn);

% Define eigenfunction matrices
% -------------------------------------------------------------------------
Phi_IM = cos(pi*x_I*m_M');
Psi_IM = Phi_IM*D_MM;
Phi_JN = cos(pi*y_J*n_N');
Psi_JN = Phi_JN*D_NN;

% Compute the Fourier cosine components A_MN
% -------------------------------------------------------------------------
A_MN = Psi_IM'*f_IJ*Psi_JN*dx*dy;

% Compute the approximation of f_IJ
% -------------------------------------------------------------------------
ff_IJ = Phi_IM*A_MN*Phi_JN';
end

