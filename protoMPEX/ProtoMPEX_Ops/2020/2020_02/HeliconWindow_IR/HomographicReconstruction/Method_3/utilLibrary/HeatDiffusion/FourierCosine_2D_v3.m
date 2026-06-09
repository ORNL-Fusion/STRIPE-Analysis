function [ff_IJ,A_NM,Phi_JM,Phi_IN] = FourierCosine_2D_v3(x_J,y_I,f_IJ,M,N)
% FOURIERCOSINE_2D_V3 takes a 2D array F_IJ and computes the Fourier cosine
% coefficients. The domain of X_J and Y_I are [0,1]
% In this version, we change the dimensions of x and y so that x has the
% same elements as columns in f_IJ and y has the same number of elements as
% rows as f_IJ

% CONVENTION: Indices denote the number of elements in an object in terms
% of rows and colums, for example, f_IJ denotes an array with "I" rows and
% "J" columns. The associated coordinate variables "x" and "y" are defined
% as "x_J" and "y_I" so as to have the correct dimensions, thus x_1 and y_2
% refer to element f_21. This allows to correctly perform:
% surf(x_J,y_I,f_IJ) since size(f_IJ) = size(y*x')
%
% INPUTS:
% X_J: Column vector of "x" values with "J" elements, [0,1]
% Y_I: Column vector of "y" values with "I" elements, [0,1]
% F_IJ: 2D matrix of function to Fourier decompose
% M: Number of terms in the p(x) eigenfunction sum
% N: Number of terms in the p(y) eigenfunction sum
% OUTPUTS:
% A_NM: 2D array of Fourier Cosine Coefficients, NxM elements
% Phi_JM: 2D array of "cos(m*pi*x)" eigenfunctions for all "m" and all "x"
% Phi_IN: 2D array of "cos(n*pi*y)" eigenfunctions for all "n" and all "y"
% ff_IJ : Approximated version of f_IJ

% J.F. Caneses Marin
% Created 2019_03_13
% Modified 2020_04_07

% Check size of independent variables
% -------------------------------------------------------------------------
J = length(x_J);
I = length(y_I);

dx = x_J(2) - x_J(1);
dy = y_I(2) - y_I(1);

% Check dimensions of inputs:
% -------------------------------------------------------------------------
if ~iscolumn(x_J)
    x_J = x_J';
end
if ~iscolumn(y_I)
    y_I = y_I';
end

if size(f_IJ) ~= [I,J]
   error('size of f_IJ not consistet with x_J and y_I');
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
Phi_JM = cos(pi*x_J*m_M');
Psi_JM = Phi_JM*D_MM;
Phi_IN = cos(pi*y_I*n_N');
Psi_IN = Phi_IN*D_NN;

% Compute the Fourier cosine components A_NM
% -------------------------------------------------------------------------
A_NM = Psi_IN'*f_IJ*Psi_JM*dx*dy;

% Compute the approximation of f_IJ
% -------------------------------------------------------------------------
ff_IJ = Phi_IN*A_NM*Phi_JM';
end

